#include "uint256.h"
#include <cstdio>
#include <cassert>

using namespace secp256k1_field;

void test_addition() {
    uint256_t a(0, 0, 0, 1);
    uint256_t b(0, 0, 0, 2);
    uint256_t c = uint256_add(a, b);
    assert(c.v[0] == 3);
    assert(c.v[1] == 0);
    
    // Test carry
    uint256_t d(0, 0, 0, 0xFFFFFFFFFFFFFFFFULL);
    uint256_t e = uint256_add(d, a);
    assert(e.v[0] == 0);
    assert(e.v[1] == 1);
    printf("Addition tests passed.\n");
}

void test_modular_arithmetic() {
    // p = FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
    uint256_t a = uint256_sub(get_p(), uint256_t(5ULL));
    uint256_t b = uint256_t(10ULL);
    
    uint256_t c = mod_add(a, b); // -5 + 10 = 5 mod p
    assert(c.v[0] == 5 && c.v[1] == 0 && c.v[2] == 0 && c.v[3] == 0);
    
    uint256_t d = mod_sub(uint256_t(5ULL), uint256_t(10ULL)); // -5 mod p = p - 5
    assert(d == a);
    
    printf("Modular add/sub tests passed.\n");
}

void test_modular_multiplication() {
    // 2 * 3 = 6
    uint256_t a = uint256_t(2ULL);
    uint256_t b = uint256_t(3ULL);
    uint256_t c = mod_mul(a, b);
    assert(c.v[0] == 6 && c.v[1] == 0);
    
    // (p-1) * (p-1) = 1 mod p
    uint256_t p_minus_1 = uint256_sub(get_p(), uint256_t(1ULL));
    uint256_t d = mod_mul(p_minus_1, p_minus_1);
    assert(d.v[0] == 1 && d.v[1] == 0 && d.v[2] == 0 && d.v[3] == 0);
    
    printf("Modular multiplication tests passed.\n");
}

void test_modular_inversion() {
    uint256_t a = uint256_t(5ULL);
    uint256_t a_inv = mod_inv(a);
    
    uint256_t check = mod_mul(a, a_inv);
    assert(check.v[0] == 1 && check.v[1] == 0);
    
    printf("Modular inversion tests passed.\n");
}

int main() {
    test_addition();
    test_modular_arithmetic();
    test_modular_multiplication();
    test_modular_inversion();
    printf("All math tests passed!\n");
    return 0;
}
