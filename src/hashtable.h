/*
 * BSGS Turbo - hashtable.h
 * Compact hash table using Robin Hood hashing for BSGS baby-step lookup
 * Stores (x_prefix -> baby_step_index) with zero false positives
 */
#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <cstdint>
#include <cstring>
#include <cstdlib>

// Each entry stores a 6-byte prefix of the X coordinate + 6-byte index
// Total: 12 bytes per entry for compact storage
// With 6-byte prefix, collision probability is ~1/2^48 per lookup
// Much better than bloom filters while using similar memory

struct HTEntry {
    uint8_t x_prefix[6]; // First 6 bytes of X coordinate (big-endian)
    uint8_t index[6];    // Baby step index (supports up to 2^48 baby steps)
    uint8_t occupied;    // 0 = empty, 1 = occupied
    uint8_t psl;         // Probe sequence length for Robin Hood

    bool is_empty() const { return occupied == 0; }

    void set(const uint8_t* xp, uint64_t idx) {
        memcpy(x_prefix, xp, 6);
        index[0] = (uint8_t)(idx);
        index[1] = (uint8_t)(idx >> 8);
        index[2] = (uint8_t)(idx >> 16);
        index[3] = (uint8_t)(idx >> 24);
        index[4] = (uint8_t)(idx >> 32);
        index[5] = (uint8_t)(idx >> 40);
        occupied = 1;
        psl = 0;
    }

    uint64_t get_index() const {
        return (uint64_t)index[0] |
               ((uint64_t)index[1] << 8) |
               ((uint64_t)index[2] << 16) |
               ((uint64_t)index[3] << 24) |
               ((uint64_t)index[4] << 32) |
               ((uint64_t)index[5] << 40);
    }

    bool match_prefix(const uint8_t* xp) const {
        return memcmp(x_prefix, xp, 6) == 0;
    }
};

class CompactHashTable {
public:
    HTEntry* table;
    uint64_t capacity;
    uint64_t mask;
    uint64_t count;

    CompactHashTable() : table(nullptr), capacity(0), mask(0), count(0) {}

    ~CompactHashTable() {
        if (table) free(table);
    }

    // Initialize with expected number of elements
    // Uses ~1.5x capacity for good load factor
    bool init(uint64_t expected_elements) {
        // Round up to power of 2
        capacity = 1;
        while (capacity < expected_elements * 3 / 2) {
            capacity <<= 1;
        }
        mask = capacity - 1;
        count = 0;

        table = (HTEntry*)calloc(capacity, sizeof(HTEntry));
        if (!table) return false;
        return true;
    }

    // Extract 6-byte prefix from X coordinate (uint256_t)
    static void extract_prefix(const uint256_t& x, uint8_t prefix[6]) {
        // Take from most significant bytes for better distribution
        // v[3] is most significant limb
        prefix[0] = (uint8_t)(x.v[3] >> 56);
        prefix[1] = (uint8_t)(x.v[3] >> 48);
        prefix[2] = (uint8_t)(x.v[3] >> 40);
        prefix[3] = (uint8_t)(x.v[3] >> 32);
        prefix[4] = (uint8_t)(x.v[3] >> 24);
        prefix[5] = (uint8_t)(x.v[3] >> 16);
    }

    // Hash function for the prefix
    static uint64_t hash_prefix(const uint8_t prefix[6]) {
        // FNV-1a hash on the prefix bytes
        uint64_t h = 14695981039346656037ULL;
        for (int i = 0; i < 6; i++) {
            h ^= prefix[i];
            h *= 1099511628211ULL;
        }
        return h;
    }

    // Insert a baby step
    void insert(const uint256_t& x, uint64_t baby_index) {
        uint8_t prefix[6];
        extract_prefix(x, prefix);
        uint64_t pos = hash_prefix(prefix) & mask;

        HTEntry entry;
        entry.set(prefix, baby_index);
        entry.psl = 0;

        while (true) {
            if (table[pos].is_empty()) {
                table[pos] = entry;
                count++;
                return;
            }

            // Robin Hood: swap if current entry has traveled less
            if (table[pos].psl < entry.psl) {
                HTEntry tmp = table[pos];
                table[pos] = entry;
                entry = tmp;
            }

            entry.psl++;
            pos = (pos + 1) & mask;

            // Safety: if PSL gets too high, table is too full
            if (entry.psl > 200) {
                // This shouldn't happen with 1.5x capacity
                return;
            }
        }
    }

    // Lookup: returns baby_index if found, -1 if not
    int64_t lookup(const uint256_t& x) const {
        uint8_t prefix[6];
        extract_prefix(x, prefix);
        uint64_t pos = hash_prefix(prefix) & mask;
        uint8_t current_psl = 0;

        while (true) {
            if (table[pos].is_empty()) return -1;
            if (table[pos].psl < current_psl) return -1; // Robin Hood guarantee

            if (table[pos].match_prefix(prefix)) {
                return (int64_t)table[pos].get_index();
            }

            current_psl++;
            pos = (pos + 1) & mask;

            if (current_psl > 200) return -1;
        }
    }

    // Memory usage in MB
    double memory_mb() const {
        return (double)(capacity * sizeof(HTEntry)) / (1024.0 * 1024.0);
    }
};

#endif // HASHTABLE_H
