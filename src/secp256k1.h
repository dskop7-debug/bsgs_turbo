/*
 * BSGS Turbo - secp256k1.h
 * Elliptic curve operations on secp256k1 using Jacobian coordinates
 * Supports batch point operations and endomorphism
 */
#ifndef SECP256K1_H
#define SECP256K1_H

#include "uint256.h"
#include <vector>

using namespace secp256k1_field;

// ==================== Point Types ====================

struct PointAffine {
    uint256_t x, y;
    bool infinity;

    PointAffine() : infinity(true) {}
    PointAffine(const uint256_t& x_, const uint256_t& y_) : x(x_), y(y_), infinity(false) {}
};

struct PointJacobian {
    uint256_t X, Y, Z;

    PointJacobian() : X(0ULL), Y(1ULL), Z(0ULL) {} // Point at infinity: Z=0
    PointJacobian(const uint256_t& x, const uint256_t& y, const uint256_t& z)
        : X(x), Y(y), Z(z) {}

    bool is_infinity() const { return Z.is_zero(); }

    static PointJacobian from_affine(const PointAffine& p) {
        if (p.infinity) return PointJacobian();
        return PointJacobian(p.x, p.y, uint256_t(1ULL));
    }
};

// ==================== secp256k1 Constants ====================
namespace secp256k1_constants {

// Generator point G
inline const PointAffine& get_G() {
    static const PointAffine G(
        uint256_t::from_hex("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798"),
        uint256_t::from_hex("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8")
    );
    return G;
}

// Endomorphism constants
// beta: cube root of unity in Fp, beta^3 = 1 mod p
// lambda: cube root of unity in Fn, lambda^3 = 1 mod n
// If P = (x, y), then lambda*P = (beta*x, y)
inline const uint256_t& get_beta() {
    static const uint256_t beta = uint256_t::from_hex(
        "7AE96A2B657C07106E64479EAC3434E99CF0497512F58995C1396C28719501EE"
    );
    return beta;
}

inline const uint256_t& get_lambda() {
    static const uint256_t lambda = uint256_t::from_hex(
        "5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72"
    );
    return lambda;
}

// Second endomorphism values (beta^2, lambda^2)
inline const uint256_t& get_beta2() {
    static const uint256_t beta2 = uint256_t::from_hex(
        "851695D49A83F8EF919BB86153CBCB16630FB68AED0A766A3EC693D68E6AFA40"
    );
    return beta2;
}

inline const uint256_t& get_lambda2() {
    static const uint256_t lambda2 = uint256_t::from_hex(
        "AC9C52B33FA3CF1F5AD9E3FD77ED9BA4A880B9FC8EC739C2E0CFC810B51283CF"
    );
    return lambda2;
}

} // namespace secp256k1_constants

// ==================== EC Point Operations ====================

// Point doubling in Jacobian: 2P
// Cost: 1S + 5M (S=squaring, M=multiplication)
inline PointJacobian point_double(const PointJacobian& p) {
    if (p.is_infinity()) return p;

    // For secp256k1, a=0 in y^2 = x^3 + 7
    // S = 4*X*Y^2
    // M = 3*X^2 (since a=0)
    // X' = M^2 - 2*S
    // Y' = M*(S - X') - 8*Y^4
    // Z' = 2*Y*Z

    uint256_t YY = mod_sqr(p.Y);         // Y^2
    uint256_t YYYY = mod_sqr(YY);         // Y^4
    uint256_t XX = mod_sqr(p.X);          // X^2
    uint256_t S = mod_mul(p.X, YY);       // X*Y^2
    S = mod_add(S, S);
    S = mod_add(S, S);                    // S = 4*X*Y^2

    uint256_t M = mod_add(XX, mod_add(XX, XX));  // M = 3*X^2

    uint256_t X3 = mod_sqr(M);
    X3 = mod_sub(X3, mod_add(S, S));      // X' = M^2 - 2*S

    uint256_t Y3 = mod_sub(S, X3);
    Y3 = mod_mul(M, Y3);
    uint256_t YYYY8 = mod_add(YYYY, YYYY);
    YYYY8 = mod_add(YYYY8, YYYY8);
    YYYY8 = mod_add(YYYY8, YYYY8);        // 8*Y^4
    Y3 = mod_sub(Y3, YYYY8);             // Y' = M*(S - X') - 8*Y^4

    uint256_t Z3 = mod_mul(p.Y, p.Z);
    Z3 = mod_add(Z3, Z3);                // Z' = 2*Y*Z

    return PointJacobian(X3, Y3, Z3);
}

// Point addition in Jacobian: P + Q (where Q is in affine)
// This is faster than general Jacobian addition: 7M + 4S vs 12M + 4S
inline PointJacobian point_add_mixed(const PointJacobian& p, const PointAffine& q) {
    if (q.infinity) return p;
    if (p.is_infinity()) return PointJacobian::from_affine(q);

    // Using "madd-2008-g" formulas
    uint256_t Z1Z1 = mod_sqr(p.Z);           // Z1^2
    uint256_t U2 = mod_mul(q.x, Z1Z1);       // U2 = X2*Z1^2
    uint256_t S2 = mod_mul(q.y, mod_mul(p.Z, Z1Z1)); // S2 = Y2*Z1^3

    uint256_t H = mod_sub(U2, p.X);          // H = U2 - X1
    uint256_t HH = mod_sqr(H);              // H^2
    uint256_t I = mod_add(HH, HH);
    I = mod_add(I, I);                       // I = 4*H^2
    uint256_t J = mod_mul(H, I);             // J = H*I
    uint256_t rr = mod_sub(S2, p.Y);
    rr = mod_add(rr, rr);                   // r = 2*(S2 - Y1)
    uint256_t V = mod_mul(p.X, I);           // V = X1*I

    uint256_t X3 = mod_sqr(rr);
    X3 = mod_sub(X3, J);
    X3 = mod_sub(X3, mod_add(V, V));         // X3 = r^2 - J - 2*V

    uint256_t Y3 = mod_sub(V, X3);
    Y3 = mod_mul(rr, Y3);
    uint256_t t = mod_mul(p.Y, J);
    t = mod_add(t, t);
    Y3 = mod_sub(Y3, t);                    // Y3 = r*(V - X3) - 2*Y1*J

    uint256_t Z3 = mod_add(p.Z, H);
    Z3 = mod_sqr(Z3);
    Z3 = mod_sub(Z3, Z1Z1);
    Z3 = mod_sub(Z3, HH);                   // Z3 = (Z1+H)^2 - Z1^2 - H^2

    return PointJacobian(X3, Y3, Z3);
}

// General Jacobian point addition
inline PointJacobian point_add(const PointJacobian& p1, const PointJacobian& p2) {
    if (p1.is_infinity()) return p2;
    if (p2.is_infinity()) return p1;

    uint256_t Z1Z1 = mod_sqr(p1.Z);
    uint256_t Z2Z2 = mod_sqr(p2.Z);
    uint256_t U1 = mod_mul(p1.X, Z2Z2);
    uint256_t U2 = mod_mul(p2.X, Z1Z1);
    uint256_t S1 = mod_mul(p1.Y, mod_mul(p2.Z, Z2Z2));
    uint256_t S2 = mod_mul(p2.Y, mod_mul(p1.Z, Z1Z1));

    uint256_t H = mod_sub(U2, U1);
    uint256_t R = mod_sub(S2, S1);

    if (H.is_zero()) {
        if (R.is_zero()) {
            return point_double(p1); // P == Q
        }
        return PointJacobian(); // P == -Q, return infinity
    }

    uint256_t HH = mod_sqr(H);
    uint256_t HHH = mod_mul(H, HH);
    uint256_t V = mod_mul(U1, HH);

    uint256_t X3 = mod_sqr(R);
    X3 = mod_sub(X3, HHH);
    X3 = mod_sub(X3, mod_add(V, V));

    uint256_t Y3 = mod_sub(V, X3);
    Y3 = mod_mul(R, Y3);
    Y3 = mod_sub(Y3, mod_mul(S1, HHH));

    uint256_t Z3 = mod_mul(p1.Z, p2.Z);
    Z3 = mod_mul(Z3, H);

    return PointJacobian(X3, Y3, Z3);
}

// Convert Jacobian to Affine
inline PointAffine to_affine(const PointJacobian& p) {
    if (p.is_infinity()) return PointAffine();

    uint256_t z_inv = mod_inv(p.Z);
    uint256_t z_inv2 = mod_sqr(z_inv);
    uint256_t z_inv3 = mod_mul(z_inv2, z_inv);

    return PointAffine(mod_mul(p.X, z_inv2), mod_mul(p.Y, z_inv3));
}

// Batch convert Jacobian to Affine (using batch inversion)
inline void to_affine_batch(const PointJacobian* points, PointAffine* out, int n) {
    if (n == 0) return;

    // Collect Z values
    uint256_t* z_vals = new uint256_t[n];
    uint256_t* z_invs = new uint256_t[n];

    int valid_count = 0;
    int* valid_indices = new int[n];

    for (int i = 0; i < n; i++) {
        if (!points[i].is_infinity()) {
            z_vals[valid_count] = points[i].Z;
            valid_indices[valid_count] = i;
            valid_count++;
        } else {
            out[i] = PointAffine();
        }
    }

    if (valid_count > 0) {
        mod_inv_batch(z_vals, z_invs, valid_count);

        for (int k = 0; k < valid_count; k++) {
            int i = valid_indices[k];
            uint256_t zi2 = mod_sqr(z_invs[k]);
            uint256_t zi3 = mod_mul(zi2, z_invs[k]);
            out[i] = PointAffine(
                mod_mul(points[i].X, zi2),
                mod_mul(points[i].Y, zi3)
            );
        }
    }

    delete[] z_vals;
    delete[] z_invs;
    delete[] valid_indices;
}

// Scalar multiplication: k * P
inline PointJacobian scalar_mult(const PointAffine& p, const uint256_t& k) {
    if (k.is_zero() || p.infinity) return PointJacobian();

    PointJacobian result;
    int top = k.highest_bit();

    for (int i = top; i >= 0; i--) {
        result = point_double(result);
        if (k.get_bit(i)) {
            result = point_add_mixed(result, p);
        }
    }
    return result;
}

// ==================== Endomorphism ====================
// For secp256k1: given P = (x, y), we have lambda*P = (beta*x, y)
// This means: if we're looking for k such that k*G = Target,
// we can also check (k*lambda mod n) and (k*lambda^2 mod n) for free

// Apply endomorphism: returns (beta*x, y))
inline PointAffine endo_apply(const PointAffine& p) {
    if (p.infinity) return p;
    return PointAffine(mod_mul(secp256k1_constants::get_beta(), p.x), p.y);
}

// Apply second endomorphism: returns (beta^2 * x, y)
inline PointAffine endo_apply2(const PointAffine& p) {
    if (p.infinity) return p;
    return PointAffine(mod_mul(secp256k1_constants::get_beta2(), p.x), p.y);
}

// ==================== Pubkey Parsing ====================

// Parse compressed or uncompressed public key from hex string
inline PointAffine parse_pubkey(const char* hex) {
    size_t len = strlen(hex);

    if (len == 66) {
        // Compressed: 02/03 + 32-byte x
        uint8_t prefix = 0;
        if (hex[0] == '0' && (hex[1] == '2' || hex[1] == '3')) {
            prefix = hex[1] - '0';
        } else {
            return PointAffine(); // invalid
        }

        uint256_t x = uint256_t::from_hex(hex + 2);

        // Compute y from x: y^2 = x^3 + 7
        uint256_t y2 = mod_add(mod_mul(mod_sqr(x), x), uint256_t(7ULL));

        // Square root: y = y2^((p+1)/4) since p ≡ 3 (mod 4)
        // (p+1)/4 = 0x3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFFFF0C
        uint256_t exp = uint256_t::from_hex(
            "3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFFFF0C"
        );
        uint256_t y(1ULL);
        uint256_t base = y2;
        int top = exp.highest_bit();
        for (int i = top; i >= 0; i--) {
            y = mod_sqr(y);
            if (exp.get_bit(i)) {
                y = mod_mul(y, base);
            }
        }

        // Check parity and negate if needed
        if ((y.v[0] & 1) != (prefix & 1)) {
            y = mod_neg(y);
        }

        // Verify: y^2 should equal y2
        if (mod_sqr(y) != y2) {
            return PointAffine(); // invalid point
        }

        return PointAffine(x, y);

    } else if (len == 130) {
        // Uncompressed: 04 + 32-byte x + 32-byte y
        if (hex[0] != '0' || hex[1] != '4') return PointAffine();
        uint256_t x = uint256_t::from_hex(hex + 2);
        uint256_t y = uint256_t::from_hex(hex + 66);
        return PointAffine(x, y);
    }

    return PointAffine(); // invalid format
}

// Format pubkey as compressed hex
inline std::string pubkey_to_hex(const PointAffine& p) {
    if (p.infinity) return "(infinity)";
    char buf[67];
    snprintf(buf, sizeof(buf), "%02x%s",
        (p.y.v[0] & 1) ? 3 : 2,
        p.x.to_hex().c_str());
    // Pad x to 64 chars
    std::string xhex = p.x.to_hex();
    while (xhex.size() < 64) xhex = "0" + xhex;
    return std::string(1, (p.y.v[0] & 1) ? '0' : '0') +
           std::string(1, (p.y.v[0] & 1) ? '3' : '2') + xhex;
}


#endif // SECP256K1_H
