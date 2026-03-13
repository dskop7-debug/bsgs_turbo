#include "secp256k1.h"
#include <cstdio>
#include <cassert>

using namespace secp256k1_field;
using namespace secp256k1_constants;

void test_point_double() {
    PointAffine G = get_G();
    PointJacobian G_jac = PointJacobian::from_affine(G);
    
    PointJacobian G2_jac = point_double(G_jac);
    PointAffine G2 = to_affine(G2_jac);
    
    // 2G x-coordinate is expected to be c6047f9441ed7d6d3045406e95c07cd85c778e4b8cef3ca7abac09b95c709ee5
    std::string expected_x = "c6047f9441ed7d6d3045406e95c07cd85c778e4b8cef3ca7abac09b95c709ee5";
    assert(G2.x.to_hex() == expected_x);
    
    printf("Point doubling test passed.\n");
}

void test_point_addition() {
    PointAffine G = get_G();
    PointJacobian G_jac = PointJacobian::from_affine(G);
    
    PointJacobian G2_jac = point_double(G_jac);
    PointJacobian G3_jac = point_add_mixed(G2_jac, G);
    PointAffine G3 = to_affine(G3_jac);
    
    // 3G x-coordinate is expected to be f9308a019258c31049344f85f89d5229b531c845836f99b08601f113bce036f9
    std::string expected_x = "f9308a019258c31049344f85f89d5229b531c845836f99b08601f113bce036f9";
    assert(G3.x.to_hex() == expected_x);
    
    printf("Point addition test passed.\n");
}

void test_scalar_mult() {
    // 5 * G
    PointJacobian G5_jac = scalar_mult(get_G(), uint256_t(5ULL));
    PointAffine G5 = to_affine(G5_jac);
    
    // 5G x = 2f8bde4d1a07209355b4a7250a5c5128e88b84bddc619ab7cba8d569b240efe4
    std::string expected_x = "2f8bde4d1a07209355b4a7250a5c5128e88b84bddc619ab7cba8d569b240efe4";
    if (G5.x.to_hex() != expected_x) {
        printf("FAILED! Expected: %s\nActual  : %s\n", expected_x.c_str(), G5.x.to_hex().c_str());
    }
    assert(G5.x.to_hex() == expected_x);
    
    printf("Scalar multiplication test passed.\n");
}

void test_endomorphism() {
    PointAffine G = get_G();
    PointAffine endo1 = endo_apply(G);
    PointAffine endo2 = endo_apply2(G);
    
    PointJacobian lam_G = scalar_mult(G, get_lambda());
    PointAffine lam_G_aff = to_affine(lam_G);
    
    assert(endo1.x == lam_G_aff.x && endo1.y == lam_G_aff.y);
    printf("Endomorphism test passed.\n");
}

int main() {
    test_point_double();
    test_point_addition();
    test_scalar_mult();
    test_endomorphism();
    printf("All EC tests passed!\n");
    return 0;
}
