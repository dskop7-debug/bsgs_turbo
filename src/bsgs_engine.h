/*
 * BSGS Turbo - bsgs_engine.h
 * Baby-step Giant-step engine with endomorphism and multi-threading
 */
#ifndef BSGS_ENGINE_H
#define BSGS_ENGINE_H

#include "uint256.h"
#include "secp256k1.h"
#include "hashtable.h"
#include <vector>
#include <string>
#include <atomic>
#include <mutex>
#include <functional>

struct BSGSConfig {
    int num_threads;        // Number of worker threads
    uint64_t baby_steps;    // Number of baby steps (m = sqrt(N))
    uint256_t range_start;  // Start of search range
    uint256_t range_end;    // End of search range
    bool random_mode;       // Random subrange selection
    int stats_interval;     // Seconds between stats output (0 = disabled)
    bool quiet;             // Suppress per-thread output
    bool use_endomorphism;  // Use secp256k1 endomorphism (3x speedup)
    std::string output_file;// Output file for found keys
};

struct BSGSResult {
    uint256_t privkey;
    PointAffine pubkey;
    int target_index;       // Which target pubkey was solved
};

class BSGSEngine {
public:
    BSGSEngine();
    ~BSGSEngine();

    // Initialize the engine and precompute baby steps
    bool init(const BSGSConfig& config, const std::vector<PointAffine>& targets);

    // Run the search
    bool run();

    // Get results
    const std::vector<BSGSResult>& get_results() const { return results_; }

    // Get current speed in keys/sec
    double get_speed() const;

    // Get total keys checked
    uint64_t get_total_keys() const;

private:
    BSGSConfig config_;
    std::vector<PointAffine> targets_;
    CompactHashTable baby_table_;

    // Precomputed values
    PointAffine G_;           // Generator
    PointAffine giant_step_;  // m * G (the giant step increment)
    uint64_t baby_count_;     // Number of baby steps stored

    PointAffine* giant_batch_ = nullptr;
    uint256_t* giant_step_keys_ = nullptr;
    PointAffine neg_batch_giant_;

    // Also store endomorphism variants of targets
    std::vector<PointAffine> targets_endo1_;  // beta*x targets
    std::vector<PointAffine> targets_endo2_;  // beta^2*x targets

    // Results and synchronization
    std::vector<BSGSResult> results_;
    std::mutex results_mutex_;
    std::atomic<uint64_t> total_keys_;
    std::atomic<bool> all_found_;
    std::atomic<int> found_count_;
    std::atomic<bool>* target_found_;

    // Worker thread function
    void worker_sequential(int thread_id, uint256_t start, uint256_t end);
    void worker_random(int thread_id);

    // Add a result
    void add_result(const uint256_t& privkey, const PointAffine& pubkey, int target_idx);

    // Reconstruct and verify match
    void check_match_affine(const uint256_t& key_checked, const uint256_t& X_match, int target_idx);

    // Check giant step point against hash table for all targets
    bool check_point(const PointAffine& point, const uint256_t& giant_key, int target_idx);
};

#endif // BSGS_ENGINE_H
