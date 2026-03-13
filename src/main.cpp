/*
 * BSGS Turbo - main.cpp
 * High-performance BSGS key search tool for secp256k1
 * Optimized version inspired by keyhunt's BSGS mode
 *
 * Usage: bsgs_turbo -f <pubkeys.txt> -b <bits> [-t <threads>] [-n <value>] [-k <factor>]
 *                   [-r <from:to>] [-R] [-s <stats_seconds>] [-q] [-e]
 */
#include "bsgs_engine.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <thread>

static void print_banner() {
    printf("====================================================\n");
    printf("  BSGS Turbo v1.0 - High Performance Key Search\n");
    printf("  secp256k1 Baby-step Giant-step with optimizations\n");
    printf("  Based on keyhunt by AlbertoBSD\n");
    printf("====================================================\n");
}

static void print_usage(const char* prog) {
    printf("\nUsage: %s [options]\n\n", prog);
    printf("Required:\n");
    printf("  -f <file>      File with target public keys (one per line)\n");
    printf("  -b <bits>      Bit range to search (e.g., 63 for puzzle #63)\n");
    printf("  OR\n");
    printf("  -r <from:to>   Explicit hex range (e.g., 4000000000000000:8000000000000000)\n");
    printf("\nOptional:\n");
    printf("  -t <threads>   Number of threads (default: CPU cores)\n");
    printf("  -n <value>     Baby step count in hex (default: 0x100000 = 2^20)\n");
    printf("  -k <factor>    K factor to multiply n (like keyhunt, default: 1)\n");
    printf("  -R             Random mode (random subranges instead of sequential)\n");
    printf("  -e             Enable endomorphism (3x speed boost)\n");
    printf("  -s <seconds>   Stats interval in seconds (0 = off, default: 10)\n");
    printf("  -q             Quiet mode (suppress per-thread output)\n");
    printf("  -o <file>      Output file for found keys (default: KEYFOUND.txt)\n");
    printf("  -h             Show this help\n");
    printf("\nExamples:\n");
    printf("  %s -f puzzle63.txt -b 63 -t 8 -e\n", prog);
    printf("  %s -f targets.txt -r 4000000000000000:8000000000000000 -t 4 -n 0x1000000 -k 16 -e\n", prog);
    printf("\n");
}

// Load public keys from file
static std::vector<PointAffine> load_pubkeys(const char* filename) {
    std::vector<PointAffine> keys;
    std::ifstream file(filename);
    if (!file.is_open()) {
        fprintf(stderr, "[!] Cannot open file: %s\n", filename);
        return keys;
    }

    std::string line;
    while (std::getline(file, line)) {
        // Trim whitespace
        size_t start = line.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) continue;
        line = line.substr(start);

        // Skip comments
        if (line[0] == '#') continue;

        // Take first token (pubkey hex), ignore rest (comments, privkey hints)
        size_t end = line.find_first_of(" \t#");
        std::string pubkey_hex = (end != std::string::npos) ? line.substr(0, end) : line;

        if (pubkey_hex.size() != 66 && pubkey_hex.size() != 130) {
            fprintf(stderr, "[!] Skipping invalid pubkey (len=%zu): %s\n",
                pubkey_hex.size(), pubkey_hex.c_str());
            continue;
        }

        PointAffine p = parse_pubkey(pubkey_hex.c_str());
        if (p.infinity) {
            fprintf(stderr, "[!] Failed to parse pubkey: %s\n", pubkey_hex.c_str());
            continue;
        }

        keys.push_back(p);
    }

    return keys;
}

int main(int argc, char* argv[]) {
    print_banner();

    // Default configuration
    BSGSConfig config;
    config.num_threads = (int)std::thread::hardware_concurrency();
    if (config.num_threads == 0) config.num_threads = 4;
    config.baby_steps = 0x100000; // 2^20 = 1M default
    config.random_mode = false;
    config.stats_interval = 10;
    config.quiet = false;
    config.use_endomorphism = false;
    config.output_file = "KEYFOUND.txt";

    const char* pubkey_file = nullptr;
    int bit_range = 0;
    bool has_range = false;
    uint256_t range_from, range_to;
    uint64_t k_factor = 1;
    uint64_t n_value = 0;
    bool n_specified = false;

    // Parse arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-f") == 0 && i + 1 < argc) {
            pubkey_file = argv[++i];
        } else if (strcmp(argv[i], "-b") == 0 && i + 1 < argc) {
            bit_range = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-r") == 0 && i + 1 < argc) {
            i++;
            // Parse range "from:to"
            char* colon = strchr(argv[i], ':');
            if (colon) {
                *colon = 0;
                range_from = uint256_t::from_hex(argv[i]);
                range_to = uint256_t::from_hex(colon + 1);
                has_range = true;
            } else {
                fprintf(stderr, "[!] Invalid range format. Use from:to\n");
                return 1;
            }
        } else if (strcmp(argv[i], "-t") == 0 && i + 1 < argc) {
            config.num_threads = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
            i++;
            n_value = strtoull(argv[i], nullptr, 0); // Supports 0x prefix
            n_specified = true;
        } else if (strcmp(argv[i], "-k") == 0 && i + 1 < argc) {
            k_factor = strtoull(argv[++i], nullptr, 0);
        } else if (strcmp(argv[i], "-R") == 0) {
            config.random_mode = true;
        } else if (strcmp(argv[i], "-e") == 0) {
            config.use_endomorphism = true;
        } else if (strcmp(argv[i], "-s") == 0 && i + 1 < argc) {
            config.stats_interval = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-q") == 0) {
            config.quiet = true;
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            config.output_file = argv[++i];
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        } else {
            fprintf(stderr, "[!] Unknown option: %s\n", argv[i]);
            print_usage(argv[0]);
            return 1;
        }
    }

    // Validate arguments
    if (!pubkey_file) {
        fprintf(stderr, "[!] No pubkey file specified (-f)\n");
        print_usage(argv[0]);
        return 1;
    }

    if (bit_range == 0 && !has_range) {
        fprintf(stderr, "[!] No search range specified (-b or -r)\n");
        print_usage(argv[0]);
        return 1;
    }

    // Set range from bit range
    if (!has_range && bit_range > 0) {
        // Range: [2^(bits-1), 2^bits)
        range_from = uint256_t(0ULL);
        range_to = uint256_t(0ULL);
        if (bit_range <= 64) {
            range_from.v[0] = 1ULL << (bit_range - 1);
            range_to.v[0] = 1ULL << bit_range;  // Might overflow for 64
            if (bit_range == 64) {
                range_to.v[0] = 0;
                range_to.v[1] = 1;
            }
        } else if (bit_range <= 128) {
            int idx = (bit_range - 1) / 64;
            int bit = (bit_range - 1) % 64;
            range_from.v[idx] = 1ULL << bit;
            idx = bit_range / 64;
            bit = bit_range % 64;
            if (bit == 0) {
                range_to.v[idx] = 1;
            } else {
                range_to.v[idx] = 1ULL << bit;
            }
        } else if (bit_range <= 192) {
            int idx = (bit_range - 1) / 64;
            int bit = (bit_range - 1) % 64;
            range_from.v[idx] = 1ULL << bit;
            idx = bit_range / 64;
            bit = bit_range % 64;
            if (bit == 0) {
                range_to.v[idx] = 1;
            } else {
                range_to.v[idx] = 1ULL << bit;
            }
        } else {
            int idx = (bit_range - 1) / 64;
            int bit = (bit_range - 1) % 64;
            range_from.v[idx] = 1ULL << bit;
            idx = bit_range / 64;
            bit = bit_range % 64;
            if (bit == 0 && idx < 4) {
                range_to.v[idx] = 1;
            } else if (idx < 4) {
                range_to.v[idx] = 1ULL << bit;
            }
        }
    }

    config.range_start = range_from;
    config.range_end = range_to;

    // Calculate baby steps
    if (n_specified) {
        config.baby_steps = n_value;
    }
    // Apply k factor (similar to keyhunt)
    // baby_steps = sqrt(n * k_factor * 2^20) roughly... simplify:
    // Actually in keyhunt, -n sets the "N" (total range per cycle) and
    // baby_steps = sqrt(N), then k amplifies the baby table size
    // Let's match keyhunt behavior: baby_steps = sqrt(n_value) * k_factor
    if (n_specified) {
        // n_value is the cycle size, baby_steps = sqrt(n_value)
        config.baby_steps = (uint64_t)sqrt((double)n_value);
        if (config.baby_steps < 1024) config.baby_steps = 1024;
    }
    config.baby_steps *= k_factor;

    // Ensure baby_steps is reasonable
    if (config.baby_steps < 1024) config.baby_steps = 1024;
    if (config.baby_steps > (1ULL << 28)) {
        printf("[!] Warning: very large baby step count (%llu), may use lots of memory\n",
            (unsigned long long)config.baby_steps);
    }

    // Load public keys
    printf("[+] Loading public keys from: %s\n", pubkey_file);
    std::vector<PointAffine> targets = load_pubkeys(pubkey_file);
    if (targets.empty()) {
        fprintf(stderr, "[!] No valid public keys loaded\n");
        return 1;
    }
    printf("[+] Loaded %zu public keys\n", targets.size());

    printf("[+] Bit range: %d\n", bit_range);
    printf("[+] Range: 0x%s -> 0x%s\n",
        config.range_start.to_hex().c_str(),
        config.range_end.to_hex().c_str());
    printf("[+] Mode: %s\n", config.random_mode ? "random" : "sequential");

    // Create and run engine
    BSGSEngine engine;

    if (!engine.init(config, targets)) {
        fprintf(stderr, "[!] Engine initialization failed\n");
        return 1;
    }

    bool found = engine.run();

    if (found) {
        printf("\n[+] SUCCESS! Keys saved to: %s\n", config.output_file.c_str());
    } else {
        printf("\n[-] No keys found in the specified range\n");
    }

    return found ? 0 : 1;
}
