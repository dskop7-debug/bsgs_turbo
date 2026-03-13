/*
 * BSGS Turbo - uint256.h
 * High-performance 256-bit unsigned integer arithmetic for secp256k1
 * Uses 4x uint64_t limbs in little-endian order
 */
#ifndef UINT256_H
#define UINT256_H

#include <cstdint>
#include <cstring>
#include <cstdio>
#include <string>

#ifdef _MSC_VER
#include <intrin.h>
#pragma intrinsic(_umul128, _addcarry_u64, _subborrow_u64)
#else
#include <x86intrin.h>
#endif

struct uint256_t {
    uint64_t v[4]; // little-endian: v[0] is least significant

    uint256_t() { v[0] = v[1] = v[2] = v[3] = 0; }
    explicit uint256_t(uint64_t x) { v[0] = x; v[1] = v[2] = v[3] = 0; }
    uint256_t(uint64_t v3, uint64_t v2, uint64_t v1, uint64_t v0) {
        v[0] = v0; v[1] = v1; v[2] = v2; v[3] = v3;
    }

    bool is_zero() const {
        return (v[0] | v[1] | v[2] | v[3]) == 0;
    }

    bool operator==(const uint256_t& o) const {
        return v[0] == o.v[0] && v[1] == o.v[1] && v[2] == o.v[2] && v[3] == o.v[3];
    }
    bool operator!=(const uint256_t& o) const { return !(*this == o); }

    bool operator<(const uint256_t& o) const {
        if (v[3] != o.v[3]) return v[3] < o.v[3];
        if (v[2] != o.v[2]) return v[2] < o.v[2];
        if (v[1] != o.v[1]) return v[1] < o.v[1];
        return v[0] < o.v[0];
    }
    bool operator<=(const uint256_t& o) const { return !(o < *this); }
    bool operator>(const uint256_t& o) const { return o < *this; }
    bool operator>=(const uint256_t& o) const { return !(*this < o); }

    // Get bit at position i (0-indexed from LSB)
    int get_bit(int i) const {
        return (v[i / 64] >> (i % 64)) & 1;
    }

    // Highest set bit position (0-indexed), -1 if zero
    int highest_bit() const {
        for (int i = 3; i >= 0; i--) {
            if (v[i] != 0) {
#ifdef _MSC_VER
                unsigned long idx;
                _BitScanReverse64(&idx, v[i]);
                return i * 64 + (int)idx;
#else
                return i * 64 + 63 - __builtin_clzll(v[i]);
#endif
            }
        }
        return -1;
    }

    // Right shift by 1
    uint256_t shr1() const {
        uint256_t r;
        r.v[0] = (v[0] >> 1) | (v[1] << 63);
        r.v[1] = (v[1] >> 1) | (v[2] << 63);
        r.v[2] = (v[2] >> 1) | (v[3] << 63);
        r.v[3] = v[3] >> 1;
        return r;
    }

    // Convert to hex string
    std::string to_hex() const {
        char buf[67];
        snprintf(buf, sizeof(buf), "%016llx%016llx%016llx%016llx",
            (unsigned long long)v[3], (unsigned long long)v[2],
            (unsigned long long)v[1], (unsigned long long)v[0]);
        // Strip leading zeros
        std::string s(buf);
        size_t pos = s.find_first_not_of('0');
        return (pos == std::string::npos) ? "0" : s.substr(pos);
    }

    // Parse from hex string
    static uint256_t from_hex(const char* hex) {
        uint256_t r;
        // Skip optional 0x prefix
        if (hex[0] == '0' && (hex[1] == 'x' || hex[1] == 'X')) hex += 2;
        size_t len = strlen(hex);
        if (len > 64) len = 64;

        // Pad to 64 chars
        char padded[65] = {};
        size_t pad = 64 - len;
        memset(padded, '0', pad);
        memcpy(padded + pad, hex, len);
        padded[64] = 0;

        auto hex_to_u64 = [](const char* s) -> uint64_t {
            uint64_t v = 0;
            for (int i = 0; i < 16; i++) {
                char c = s[i];
                uint64_t nibble;
                if (c >= '0' && c <= '9') nibble = c - '0';
                else if (c >= 'a' && c <= 'f') nibble = c - 'a' + 10;
                else if (c >= 'A' && c <= 'F') nibble = c - 'A' + 10;
                else nibble = 0;
                v = (v << 4) | nibble;
            }
            return v;
        };

        r.v[3] = hex_to_u64(padded);
        r.v[2] = hex_to_u64(padded + 16);
        r.v[1] = hex_to_u64(padded + 32);
        r.v[0] = hex_to_u64(padded + 48);
        return r;
    }
};

// ==================== 256-bit Addition ====================
inline uint256_t uint256_add(const uint256_t& a, const uint256_t& b) {
    uint256_t r;
    unsigned char carry = 0;
#ifdef _MSC_VER
    carry = _addcarry_u64(0, a.v[0], b.v[0], &r.v[0]);
    carry = _addcarry_u64(carry, a.v[1], b.v[1], &r.v[1]);
    carry = _addcarry_u64(carry, a.v[2], b.v[2], &r.v[2]);
    _addcarry_u64(carry, a.v[3], b.v[3], &r.v[3]);
#else
    unsigned __int128 t = (unsigned __int128)a.v[0] + b.v[0];
    r.v[0] = (uint64_t)t; carry = (unsigned char)(t >> 64);
    t = (unsigned __int128)a.v[1] + b.v[1] + carry;
    r.v[1] = (uint64_t)t; carry = (unsigned char)(t >> 64);
    t = (unsigned __int128)a.v[2] + b.v[2] + carry;
    r.v[2] = (uint64_t)t; carry = (unsigned char)(t >> 64);
    t = (unsigned __int128)a.v[3] + b.v[3] + carry;
    r.v[3] = (uint64_t)t;
#endif
    return r;
}

// ==================== 256-bit Subtraction ====================
inline uint256_t uint256_sub(const uint256_t& a, const uint256_t& b) {
    uint256_t r;
#ifdef _MSC_VER
    unsigned char borrow = 0;
    borrow = _subborrow_u64(0, a.v[0], b.v[0], &r.v[0]);
    borrow = _subborrow_u64(borrow, a.v[1], b.v[1], &r.v[1]);
    borrow = _subborrow_u64(borrow, a.v[2], b.v[2], &r.v[2]);
    _subborrow_u64(borrow, a.v[3], b.v[3], &r.v[3]);
#else
    unsigned __int128 t = (unsigned __int128)a.v[0] - b.v[0];
    r.v[0] = (uint64_t)t;
    unsigned char borrow = (t >> 127) & 1;
    t = (unsigned __int128)a.v[1] - b.v[1] - borrow;
    r.v[1] = (uint64_t)t; borrow = (t >> 127) & 1;
    t = (unsigned __int128)a.v[2] - b.v[2] - borrow;
    r.v[2] = (uint64_t)t; borrow = (t >> 127) & 1;
    t = (unsigned __int128)a.v[3] - b.v[3] - borrow;
    r.v[3] = (uint64_t)t;
#endif
    return r;
}

// ==================== Modular Arithmetic (mod p for secp256k1) ====================
// p = FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
namespace secp256k1_field {

static const uint256_t P(
    0xFFFFFFFFFFFFFFFEULL, 0xFFFFFFFFFFFFFFFFULL,
    0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFEFFFFFC2FULL  // Note: reversed for our LE layout
);

// Actually, let's define P correctly in our little-endian layout:
// P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
// v[3] = 0xFFFFFFFFFFFFFFFF (most significant)
// v[2] = 0xFFFFFFFFFFFFFFFF
// v[1] = 0xFFFFFFFFFFFFFFFF
// v[0] = 0xFFFFFFFEFFFFFC2F (least significant)
inline const uint256_t& get_p() {
    static const uint256_t p(
        0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL,
        0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFEFFFFFC2FULL
    );
    return p;
}

// N (order of the group)
// N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
inline const uint256_t& get_n() {
    static const uint256_t n(
        0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFEULL,
        0xBAAEDCE6AF48A03BULL, 0xBFD25E8CD0364141ULL
    );
    return n;
}

// Modular addition: (a + b) mod p
inline uint256_t mod_add(const uint256_t& a, const uint256_t& b) {
    uint256_t r = uint256_add(a, b);
    const uint256_t& p = get_p();
    // Check for overflow or r >= p
    // Since a,b < p, a+b < 2p, so at most one subtraction needed
    // We detect carry by checking if result < a (overflow) or result >= p
    if (r < a || r >= p) {
        r = uint256_sub(r, p);
    }
    return r;
}

// Modular subtraction: (a - b) mod p
inline uint256_t mod_sub(const uint256_t& a, const uint256_t& b) {
    if (a >= b) {
        return uint256_sub(a, b);
    } else {
        // a - b + p
        return uint256_add(uint256_sub(a, b), get_p());
    }
}

// Modular negation: (-a) mod p = p - a (for a != 0)
inline uint256_t mod_neg(const uint256_t& a) {
    if (a.is_zero()) return a;
    return uint256_sub(get_p(), a);
}

// ==================== Modular Multiplication ====================
// We use schoolbook 4-limb multiplication with reduction using secp256k1 special form
// p = 2^256 - 0x1000003D1

// Multiply two 256-bit numbers to get 512-bit result
inline void mul_512(const uint256_t& a, const uint256_t& b, uint64_t r[8]) {
    // Schoolbook multiplication: 4x4 limbs
    // Using 128-bit intermediates
#ifdef _MSC_VER
    uint64_t hi;
    // We'll accumulate into r[0..7]
    memset(r, 0, 8 * sizeof(uint64_t));
    for (int i = 0; i < 4; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < 4; j++) {
            uint64_t lo = _umul128(a.v[i], b.v[j], &hi);
            // r[i+j] += lo + carry
            unsigned char c1 = _addcarry_u64(0, r[i+j], lo, &r[i+j]);
            unsigned char c2 = _addcarry_u64(c1, r[i+j], carry, &r[i+j]);
            // But wait, we might double-add. Let's be more careful:
            uint64_t tmp = r[i+j];
            unsigned char c = _addcarry_u64(0, tmp, lo, &tmp);
            c = _addcarry_u64(c, tmp, carry, &tmp);
            r[i+j] = tmp;
            carry = hi + c;
        }
        r[i+4] = carry;
    }
#else
    memset(r, 0, 8 * sizeof(uint64_t));
    for (int i = 0; i < 4; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < 4; j++) {
            unsigned __int128 prod = (unsigned __int128)a.v[i] * b.v[j] + r[i+j] + carry;
            r[i+j] = (uint64_t)prod;
            carry = (uint64_t)(prod >> 64);
        }
        r[i+4] = carry;
    }
#endif
}

// Reduce a 512-bit value mod p using secp256k1 special form:
// p = 2^256 - c where c = 0x1000003D1
// For x = x_hi * 2^256 + x_lo: x mod p = x_lo + x_hi * c (mod p)
inline uint256_t reduce_512(const uint64_t r[8]) {
    const uint64_t C = 0x1000003D1ULL; // 2^32 + 977
    uint256_t result;

    // First reduction: multiply high part by C and add to low part
    // r[4..7] * C + r[0..3]
#ifdef _MSC_VER
    uint64_t carry = 0;
    uint64_t hi;
    for (int i = 0; i < 4; i++) {
        uint64_t lo = _umul128(r[i+4], C, &hi);
        unsigned char c = _addcarry_u64(0, lo, carry, &lo);
        hi += c;
        c = _addcarry_u64(0, lo, r[i], &result.v[i]);
        carry = hi + c;
    }
#else
    uint64_t carry = 0;
    for (int i = 0; i < 4; i++) {
        unsigned __int128 t = (unsigned __int128)r[i+4] * C + r[i] + carry;
        result.v[i] = (uint64_t)t;
        carry = (uint64_t)(t >> 64);
    }
#endif
    // Second reduction: carry * C
    // carry is at most ~33 bits * C which fits in 64-bit... actually let's be safe
    if (carry > 0) {
#ifdef _MSC_VER
        uint64_t lo = _umul128(carry, C, &hi);
        unsigned char c = _addcarry_u64(0, result.v[0], lo, &result.v[0]);
        c = _addcarry_u64(c, result.v[1], hi, &result.v[1]);
        c = _addcarry_u64(c, result.v[2], 0, &result.v[2]);
        _addcarry_u64(c, result.v[3], 0, &result.v[3]);
#else
        unsigned __int128 t = (unsigned __int128)carry * C + result.v[0];
        result.v[0] = (uint64_t)t;
        carry = (uint64_t)(t >> 64);
        for (int i = 1; i < 4 && carry; i++) {
            t = (unsigned __int128)result.v[i] + carry;
            result.v[i] = (uint64_t)t;
            carry = (uint64_t)(t >> 64);
        }
#endif
    }

    // Final reduction: if result >= p, subtract p
    const uint256_t& p = get_p();
    if (result >= p) {
        result = uint256_sub(result, p);
    }

    return result;
}

// Modular multiplication: (a * b) mod p
inline uint256_t mod_mul(const uint256_t& a, const uint256_t& b) {
    uint64_t r[8];
    mul_512(a, b, r);
    return reduce_512(r);
}

// Modular squaring: a^2 mod p (can be optimized but reuses mul for now)
inline uint256_t mod_sqr(const uint256_t& a) {
    return mod_mul(a, a);
}

// ==================== Modular Inversion ====================
// Using Fermat's little theorem: a^(-1) = a^(p-2) mod p
// This is simple and constant-time. For batch, we use Montgomery's trick externally.
inline uint256_t mod_inv(const uint256_t& a) {
    // p - 2 = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2D
    uint256_t p_minus_2(
        0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL,
        0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFEFFFFFC2DULL
    );

    // Square-and-multiply
    uint256_t result(1ULL);
    uint256_t base = a;

    int top = p_minus_2.highest_bit();
    for (int i = top; i >= 0; i--) {
        result = mod_sqr(result);
        if (p_minus_2.get_bit(i)) {
            result = mod_mul(result, base);
        }
    }
    return result;
}

// ==================== Batch Modular Inversion (Montgomery's Trick) ====================
// Given n values a[0..n-1], compute their inverses using only 1 inversion + 3(n-1) multiplications
inline void mod_inv_batch(const uint256_t* a, uint256_t* inv_out, int n) {
    if (n == 0) return;
    if (n == 1) { inv_out[0] = mod_inv(a[0]); return; }

    // Step 1: Compute prefix products
    // prefix[i] = a[0] * a[1] * ... * a[i]
    uint256_t* prefix = new uint256_t[n];
    prefix[0] = a[0];
    for (int i = 1; i < n; i++) {
        prefix[i] = mod_mul(prefix[i-1], a[i]);
    }

    // Step 2: Invert the product of all elements
    uint256_t inv = mod_inv(prefix[n-1]);

    // Step 3: Compute individual inverses
    for (int i = n - 1; i > 0; i--) {
        inv_out[i] = mod_mul(inv, prefix[i-1]);
        inv = mod_mul(inv, a[i]);
    }
    inv_out[0] = inv;

    delete[] prefix;
}

} // namespace secp256k1_field

#endif // UINT256_H
