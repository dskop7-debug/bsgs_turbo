/*
 * BSGS Turbo - bsgs_engine.cpp
 * Core BSGS algorithm implementation with optimizations
 */
#include "bsgs_engine.h"
#include <thread>
#include <chrono>
#include <cstdio>
#include <cmath>
#include <random>
#include <algorithm>
#include <fstream>
#include <cstdint>

using namespace secp256k1_field;
using namespace secp256k1_constants;

static uint256_t mod_n_mul(const uint256_t& a, const uint256_t& b);

BSGSEngine::BSGSEngine()
    : baby_count_(0), total_keys_(0), all_found_(false), found_count_(0), target_found_(nullptr) {}

BSGSEngine::~BSGSEngine() {
    if (target_found_) delete[] target_found_;
}

bool BSGSEngine::init(const BSGSConfig& config, const std::vector<PointAffine>& targets) {
    config_ = config;
    targets_ = targets;
    G_ = get_G();

    if (targets_.empty()) {
        fprintf(stderr, "[!] No target public keys loaded\n");
        return false;
    }

    // Initialize found tracking
    if (target_found_) delete[] target_found_;
    target_found_ = new std::atomic<bool>[targets_.size()];
    for (size_t i = 0; i < targets_.size(); i++) {
        target_found_[i].store(false);
    }

    baby_count_ = config_.baby_steps;

    printf("[+] Initializing BSGS engine\n");
    printf("[+] Baby steps (m): %llu\n", (unsigned long long)baby_count_);
    printf("[+] Targets: %zu public keys\n", targets_.size());
    printf("[+] Threads: %d\n", config_.num_threads);
    printf("[+] Endomorphism: %s\n", config_.use_endomorphism ? "enabled (3x)" : "disabled");

    // ========== Baby Step Phase ==========
    printf("[+] Building baby step table...\n");

    auto t_start = std::chrono::high_resolution_clock::now();

    // Initialize hash table
    uint64_t effective_size = baby_count_;
    if (config_.use_endomorphism) effective_size *= 3; // 3x entries with endo
    if (!baby_table_.init(effective_size)) {
        fprintf(stderr, "[!] Failed to allocate hash table (%.2f MB needed)\n",
            (double)(effective_size * 3 / 2 * sizeof(HTEntry)) / (1024.0 * 1024.0));
        return false;
    }
    printf("[+] Hash table allocated: %.2f MB\n", baby_table_.memory_mb());

    // Precompute baby steps: i*G for i = 0..m-1
    // Use batch processing for efficiency
    // Start with G, then add G repeatedly
    const int BATCH_SIZE = 4096;
    PointJacobian* jac_batch = new PointJacobian[BATCH_SIZE];
    PointAffine* aff_batch = new PointAffine[BATCH_SIZE];
    // Store 0*G? Skip index 0 to avoid degenerate case
    // Baby step i: store point = i*G, key = i
    // We precompute in batches for batch affine conversion

    uint64_t stored = 0;
    uint64_t batch_idx = 0;

    // First point: 1*G
    jac_batch[0] = PointJacobian::from_affine(G_);
    batch_idx = 1;
    stored = 0;

    // Pre-compute 2*G for repeated addition
    PointAffine G2_affine;
    {
        PointJacobian G2 = point_double(PointJacobian::from_affine(G_));
        G2_affine = to_affine(G2);
    }

    for (uint64_t i = 1; i <= baby_count_; i++) {
        if (batch_idx == BATCH_SIZE || i == baby_count_) {
            // Convert batch to affine
            to_affine_batch(jac_batch, aff_batch, (int)batch_idx);

            // Insert into hash table
            for (uint64_t b = 0; b < batch_idx; b++) {
                uint64_t actual_index = stored + b + 1; // 1-indexed
                baby_table_.insert(aff_batch[b].x, actual_index);

                // Insert endomorphism variants
                if (config_.use_endomorphism) {
                    PointAffine endo1 = endo_apply(aff_batch[b]);
                    PointAffine endo2 = endo_apply2(aff_batch[b]);
                    // For endo variants, we use a flag in the upper bits of index
                    // to distinguish which variant matched
                    baby_table_.insert(endo1.x, actual_index | (1ULL << 46));
                    baby_table_.insert(endo2.x, actual_index | (2ULL << 46));
                }
            }

            stored += batch_idx;
            batch_idx = 0;

            // Progress
            if (stored % (baby_count_ / 10 + 1) == 0 || i == baby_count_) {
                printf("[+] Baby steps: %llu/%llu (%.0f%%)\r",
                    (unsigned long long)stored, (unsigned long long)baby_count_,
                    (double)stored / baby_count_ * 100.0);
                fflush(stdout);
            }
        }

        if (i < baby_count_) {
            // Next point: add G to current
            if (i == 1) {
                jac_batch[batch_idx] = point_double(PointJacobian::from_affine(G_));
            } else {
                jac_batch[batch_idx] = point_add_mixed(jac_batch[batch_idx > 0 ? batch_idx - 1 : 0], G_);
            }
            batch_idx++;
        }
    }
    printf("\n");

    delete[] jac_batch;
    delete[] aff_batch;

    // Compute giant step: m * G
    {
        uint256_t m_val((uint64_t)baby_count_);
        PointJacobian gj = scalar_mult(G_, m_val);
        giant_step_ = to_affine(gj);
        
        // Precompute batch of giant steps as Affine coordinates for fast subtraction lookup
        const int BATCH_SIZE = 4096;
        giant_batch_ = new PointAffine[BATCH_SIZE];
        giant_step_keys_ = new uint256_t[BATCH_SIZE];
        
        PointAffine neg_giant = giant_step_;
        neg_giant.y = mod_neg(neg_giant.y);
        
        PointJacobian* gb_jac = new PointJacobian[BATCH_SIZE];
        PointJacobian cur_gj = PointJacobian::from_affine(neg_giant);
        uint256_t cur_m = m_val;
        
        for (int i=0; i<BATCH_SIZE; i++) {
            gb_jac[i] = cur_gj;
            giant_step_keys_[i] = cur_m;
            
            cur_m = uint256_add(cur_m, m_val);
            cur_gj = point_add_mixed(cur_gj, neg_giant);
        }
        to_affine_batch(gb_jac, giant_batch_, BATCH_SIZE);
        
        // The total jump for the outer loop is BATCH_SIZE * m
        neg_batch_giant_ = giant_batch_[BATCH_SIZE - 1]; // which is BATCH_SIZE * neg_giant
        delete[] gb_jac;
    }

    // Precompute endomorphism variants of targets
    if (config_.use_endomorphism) {
        targets_endo1_.resize(targets_.size());
        targets_endo2_.resize(targets_.size());
        for (size_t i = 0; i < targets_.size(); i++) {
            targets_endo1_[i] = endo_apply(targets_[i]);
            targets_endo2_[i] = endo_apply2(targets_[i]);
        }
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(t_end - t_start).count();
    printf("[+] Baby step table built in %.2f seconds\n", elapsed);
    printf("[+] Table entries: %llu (%.2f MB)\n",
        (unsigned long long)baby_table_.count, baby_table_.memory_mb());

    return true;
}

void BSGSEngine::add_result(const uint256_t& privkey, const PointAffine& pubkey, int target_idx) {
    std::lock_guard<std::mutex> lock(results_mutex_);

    // Double-check under lock
    if (target_found_[target_idx].load()) return;
    target_found_[target_idx].store(true);
    found_count_++;

    BSGSResult r;
    r.privkey = privkey;
    r.pubkey = pubkey;
    r.target_index = target_idx;
    results_.push_back(r);

    // Print found key
    printf("[+] KEY FOUND! privkey: %s\n", privkey.to_hex().c_str());
    printf("[+] Pubkey: %s\n", pubkey_to_hex(pubkey).c_str());

    // Save to file
    if (!config_.output_file.empty()) {
        FILE* f = fopen(config_.output_file.c_str(), "a");
        if (f) {
            fprintf(f, "%s # %s\n",
                privkey.to_hex().c_str(),
                pubkey_to_hex(pubkey).c_str());
            fclose(f);
        }
    }

    // Check if all targets found
    if (found_count_ >= (int)targets_.size()) {
        all_found_.store(true);
    }
}

void BSGSEngine::check_match_affine(const uint256_t& key_checked, const uint256_t& X_match, int target_idx) {
    int64_t baby_idx = baby_table_.lookup(X_match);
    if (baby_idx >= 0) {
        uint64_t raw_idx = (uint64_t)baby_idx;
        uint64_t endo_flag = (raw_idx >> 46) & 3;
        uint64_t actual_baby = raw_idx & ((1ULL << 46) - 1);

        uint256_t baby_key = uint256_t(actual_baby);
        if (endo_flag == 1) {
            baby_key = mod_n_mul(baby_key, get_lambda2());
        } else if (endo_flag == 2) {
            baby_key = mod_n_mul(baby_key, get_lambda());
        }

        uint256_t found_key = uint256_add(key_checked, baby_key);

        PointJacobian verify_jac = scalar_mult(G_, found_key);
        PointAffine verify = to_affine(verify_jac);

        if (verify.x == targets_[target_idx].x && verify.y == targets_[target_idx].y) {
            add_result(found_key, targets_[target_idx], target_idx);
        }
    }
}

void BSGSEngine::worker_sequential(int thread_id, uint256_t start, uint256_t end) {
    uint64_t m = baby_count_;

    // For each target pubkey
    for (size_t t = 0; t < targets_.size() && !all_found_.load(); t++) {
        if (target_found_[t].load()) continue;

        // Compute: Target - start * G
        // Since baby table stores 1..m, we set current_key to start - 1 
        // to ensure we check 'start' (which corresponds to baby step 1)
        uint256_t current_key = start;
        if (!current_key.is_zero()) {
            current_key = uint256_sub(current_key, uint256_t(1ULL));
        }

        PointJacobian target_jac = PointJacobian::from_affine(targets_[t]);

        // Subtract current_key*G from target
        PointAffine startG_aff;
        {
            PointJacobian sG = scalar_mult(G_, current_key);
            startG_aff = to_affine(sG);
        }

        // point = Target - current_key*G
        PointAffine neg_startG = startG_aff;
        if (!neg_startG.infinity) {
            neg_startG.y = mod_neg(neg_startG.y);
        }

        PointJacobian base_jac = point_add_mixed(target_jac, neg_startG);
        PointAffine base_aff = to_affine(base_jac);

        if (!base_aff.infinity) {
            check_match_affine(current_key, base_aff.x, t);
        }

        const int BATCH_SIZE = 4096;
        uint256_t* dx = new uint256_t[BATCH_SIZE];
        uint256_t* inv_dx = new uint256_t[BATCH_SIZE];
        
        uint64_t steps = 0;
        
        while (current_key < end && !all_found_.load() && !target_found_[t].load()) {
            
            for (int i=0; i<BATCH_SIZE; i++) {
                dx[i] = mod_sub(giant_batch_[i].x, base_aff.x);
            }
            mod_inv_batch(dx, inv_dx, BATCH_SIZE);
            
            for (int i=0; i<BATCH_SIZE; i++) {
                if (inv_dx[i].is_zero()) continue;
                
                uint256_t key_to_check = uint256_add(current_key, giant_step_keys_[i]); 
                if (key_to_check > end) break;
                
                uint256_t dy = mod_sub(giant_batch_[i].y, base_aff.y);
                uint256_t lambda = mod_mul(dy, inv_dx[i]);
                
                uint256_t X3 = mod_sqr(lambda);
                X3 = mod_sub(X3, giant_batch_[i].x);
                X3 = mod_sub(X3, base_aff.x);
                
                int64_t baby_idx = baby_table_.lookup(X3);
                if (baby_idx >= 0) {
                    check_match_affine(key_to_check, X3, t);
                }
            }
            
            steps += BATCH_SIZE;
            total_keys_.fetch_add((uint64_t)m * BATCH_SIZE, std::memory_order_relaxed);
            
            current_key = uint256_add(current_key, giant_step_keys_[BATCH_SIZE - 1]);
            base_jac = point_add_mixed(base_jac, neg_batch_giant_);
            base_aff = to_affine(base_jac);
        }

        delete[] dx;
        delete[] inv_dx;

        if (!config_.quiet) {
            printf("[+] Thread %d: target %zu done (%llu giant steps)\n",
                thread_id, t, (unsigned long long)steps);
        }
    }
}

void BSGSEngine::worker_random(int thread_id) {
    uint64_t m = baby_count_;

    // Random number generator for this thread
    std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count() + thread_id);

    // Compute range size
    uint256_t range_size = uint256_sub(config_.range_end, config_.range_start);

    // Negate giant step for subtraction
    PointAffine neg_giant = giant_step_;
    neg_giant.y = mod_neg(neg_giant.y);
    
    // Allocate once per thread outside the infinite loop
    const int BATCH_SIZE = 4096;
    uint256_t* dx = new uint256_t[BATCH_SIZE];
    uint256_t* inv_dx = new uint256_t[BATCH_SIZE];

    while (!all_found_.load()) {
        // Pick random starting point in range
        uint256_t random_offset;
        random_offset.v[0] = rng();
        random_offset.v[1] = rng();
        random_offset.v[2] = rng();
        random_offset.v[3] = rng();

        // Reduce to range_size (simple modular reduction)
        int range_bits = range_size.highest_bit();
        if (range_bits < 255) {
            uint64_t mask_val = (range_bits >= 64) ? 0xFFFFFFFFFFFFFFFFULL : ((1ULL << (range_bits % 64)) - 1);
            for (int i = 3; i >= 0; i--) {
                if (i * 64 > range_bits) {
                    random_offset.v[i] = 0;
                } else if (i * 64 + 64 > range_bits) {
                    random_offset.v[i] &= mask_val;
                }
            }
        }
        
        // Ensure not zero if end == start somehow
        uint256_t current_key = uint256_add(config_.range_start, random_offset);
        if (!current_key.is_zero()) {
            current_key = uint256_sub(current_key, uint256_t(1ULL));
        }

        // For each target
        for (size_t t = 0; t < targets_.size() && !all_found_.load(); t++) {
            if (target_found_[t].load()) continue;

            // Compute: Target - start_key * G
            PointJacobian sG = scalar_mult(G_, current_key);
            PointAffine sG_aff = to_affine(sG);
            PointAffine neg_sG = sG_aff;
            if (!neg_sG.infinity) {
                neg_sG.y = mod_neg(neg_sG.y);
            }

            PointJacobian base_jac = point_add_mixed(PointJacobian::from_affine(targets_[t]), neg_sG);
            PointAffine base_aff = to_affine(base_jac);

            if (!base_aff.infinity) {
                check_match_affine(current_key, base_aff.x, t);
            }

            for (int i=0; i<BATCH_SIZE; i++) {
                dx[i] = mod_sub(giant_batch_[i].x, base_aff.x);
            }
            mod_inv_batch(dx, inv_dx, BATCH_SIZE);
            
            for (int i=0; i<BATCH_SIZE; i++) {
                if (inv_dx[i].is_zero()) continue;
                
                uint256_t key_to_check = uint256_add(current_key, giant_step_keys_[i]); 
                uint256_t dy = mod_sub(giant_batch_[i].y, base_aff.y);
                uint256_t lambda = mod_mul(dy, inv_dx[i]);
                
                uint256_t X3 = mod_sqr(lambda);
                X3 = mod_sub(X3, giant_batch_[i].x);
                X3 = mod_sub(X3, base_aff.x);
                
                int64_t baby_idx = baby_table_.lookup(X3);
                if (baby_idx >= 0) {
                    check_match_affine(key_to_check, X3, t);
                }
            }
            
            total_keys_.fetch_add((uint64_t)m * BATCH_SIZE, std::memory_order_relaxed);
        }
    }

    delete[] dx;
    delete[] inv_dx;
}

// Helper for mod n operations  
static uint256_t mod_n_mul(const uint256_t& a, const uint256_t& b) {
    // Simple modular multiplication mod n (group order)
    // For now use the field mul and reduce — proper implementation would
    // use dedicated mod-n reduction
    uint64_t r[8];
    secp256k1_field::mul_512(a, b, r);
    
    // Reduce mod n using simple subtraction (not as optimized as mod p)
    // This is rarely called so simplicity > speed here
    uint256_t result;
    result.v[0] = r[0]; result.v[1] = r[1]; 
    result.v[2] = r[2]; result.v[3] = r[3];
    
    // If high part is non-zero, we need full reduction
    // For baby step indices this usually fits in 256 bits
    const uint256_t& n = secp256k1_field::get_n();
    while (result >= n) {
        result = uint256_sub(result, n);
    }
    return result;
}

bool BSGSEngine::run() {
    printf("[+] Starting BSGS search...\n");
    printf("[+] Range: 0x%s -> 0x%s\n",
        config_.range_start.to_hex().c_str(),
        config_.range_end.to_hex().c_str());

    auto t_start = std::chrono::high_resolution_clock::now();
    total_keys_.store(0);

    std::vector<std::thread> threads;

    if (config_.random_mode) {
        // Random mode: each thread picks random subranges
        for (int i = 0; i < config_.num_threads; i++) {
            threads.emplace_back(&BSGSEngine::worker_random, this, i);
        }
    } else {
        // Sequential mode: divide range among threads
        uint256_t range_size = uint256_sub(config_.range_end, config_.range_start);
        
        // chunk = range_size / threads. For now, simple approximation using shr1 if power of 2
        uint256_t chunk = range_size;
        int t = config_.num_threads;
        while (t > 1) {
            chunk = chunk.shr1();
            t /= 2;
        }
        
        uint256_t thread_start = config_.range_start;
        for (int i = 0; i < config_.num_threads; i++) {
            uint256_t thread_end = (i == config_.num_threads - 1) ? config_.range_end : uint256_add(thread_start, chunk);
            
            threads.emplace_back(&BSGSEngine::worker_sequential, this, i, thread_start, thread_end);
            
            thread_start = thread_end;
        }
    }

    // Stats thread
    std::thread stats_thread;
    if (config_.stats_interval > 0) {
        stats_thread = std::thread([this, &t_start]() {
            while (!all_found_.load()) {
                std::this_thread::sleep_for(
                    std::chrono::seconds(config_.stats_interval));
                if (all_found_.load()) break;

                auto now = std::chrono::high_resolution_clock::now();
                double elapsed = std::chrono::duration<double>(now - t_start).count();
                uint64_t total = total_keys_.load();

                double speed = total / elapsed;
                const char* unit = "keys/s";
                if (speed > 1e15) { speed /= 1e15; unit = "Pkeys/s"; }
                else if (speed > 1e12) { speed /= 1e12; unit = "Tkeys/s"; }
                else if (speed > 1e9) { speed /= 1e9; unit = "Gkeys/s"; }
                else if (speed > 1e6) { speed /= 1e6; unit = "Mkeys/s"; }

                printf("[+] Speed: %.2f %s | Total: %llu | Found: %d/%zu | Time: %.1fs\n",
                    speed, unit,
                    (unsigned long long)total,
                    found_count_.load(), targets_.size(),
                    elapsed);
                fflush(stdout);
            }
        });
    }

    // Wait for all workers
    for (auto& t : threads) {
        t.join();
    }

    if (stats_thread.joinable()) {
        all_found_.store(true);
        stats_thread.join();
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(t_end - t_start).count();

    printf("[+] Search completed in %.2f seconds\n", elapsed);
    printf("[+] Total keys checked: %llu\n", (unsigned long long)total_keys_.load());
    printf("[+] Keys found: %d/%zu\n", found_count_.load(), targets_.size());

    double speed = total_keys_.load() / elapsed;
    const char* unit = "keys/s";
    if (speed > 1e15) { speed /= 1e15; unit = "Pkeys/s"; }
    else if (speed > 1e12) { speed /= 1e12; unit = "Tkeys/s"; }
    else if (speed > 1e9) { speed /= 1e9; unit = "Gkeys/s"; }
    else if (speed > 1e6) { speed /= 1e6; unit = "Mkeys/s"; }
    printf("[+] Average speed: %.2f %s\n", speed, unit);

    return found_count_.load() > 0;
}

double BSGSEngine::get_speed() const {
    return (double)total_keys_.load();
}

uint64_t BSGSEngine::get_total_keys() const {
    return total_keys_.load();
}
