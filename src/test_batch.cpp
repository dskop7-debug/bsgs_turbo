#include "secp256k1.h"
#include "hashtable.h"
#include <cstdio>
#include <cstdint>
#include <vector>

using namespace secp256k1_field;
using namespace secp256k1_constants;

static uint256_t mod_n_mul(const uint256_t& a, const uint256_t& b) {
    uint64_t r[8];
    secp256k1_field::mul_512(a, b, r);
    uint256_t result;
    result.v[0] = r[0]; result.v[1] = r[1]; 
    result.v[2] = r[2]; result.v[3] = r[3];
    const uint256_t& n = secp256k1_field::get_n();
    while (result >= n) {
        result = uint256_sub(result, n);
    }
    return result;
}

void test_endo_math() {
    PointAffine G = secp256k1_constants::get_G();
    PointJacobian GJ = PointJacobian::from_affine(G);
    
    CompactHashTable ht;
    uint64_t m = 1000;
    ht.init(m * 3);
    
    // Baby steps
    PointJacobian current = GJ;
    for (uint64_t i = 1; i <= m; i++) {
        PointAffine a = to_affine(current);
        ht.insert(a.x, i); // baby index i
        
        PointAffine endo1 = endo_apply(a);
        PointAffine endo2 = endo_apply2(a);
        
        ht.insert(endo1.x, i | (1ULL << 46));
        ht.insert(endo2.x, i | (2ULL << 46));
        
        current = point_add(current, GJ); // i+1 G
    }
    
    uint256_t target_key = uint256_t(2500ULL);
    PointAffine Target = to_affine(scalar_mult(G, target_key));
    PointAffine TargetE1 = endo_apply(Target);
    PointAffine TargetE2 = endo_apply2(Target);
    
    // Find TargetE1
    uint256_t start = uint256_t(100ULL);
    uint256_t current_key = start;
    if (!current_key.is_zero()) current_key = uint256_sub(current_key, uint256_t(1ULL));
    
    // In sequential mode, target_jac is Target. Not TargetE1!
    // The giant step loop uses Target.
    // The baby step table stores endo1.x which is (lambda * baby * G).x
    // BUT, wait!!
    // If the baby table stores `lambda * baby * G`
    // And we check `Target - giant * G` against the table...
    // We are asking: does `Target - giant * G == lambda * baby * G` ?
    // Which means `Target == giant * G + lambda * baby * G`!
    // BUT what we ACTUALLY want is to find a key such that `key * G == Target`
    // If the target is just some point, and we found it matches `lambda * baby`, that assumes Target is the point we reached.
    // Wait, endomorphism in secp256k1 means: lambda * (x, y) = (beta * x, y).
    // It's a property of the curve.
    // If a public key P = k * G, then lambda * P = (lambda * k) * G.
    // If we want to find `k` such that `k * G = Target`, 
    // keyhunt BSGS does this: 
    // It stores (baby * G).x, AND (lambda * baby * G).x, AND (lambda^2 * baby * G).x
    // Then it checks `Target - giant * G` against the table.
    // If `Target - giant * G` matches `(lambda * baby * G).x`, it means:
    // `Target - giant * G = +/-(lambda * baby * G)`
    // `Target = giant * G +/- lambda * baby * G`
    // So the private key is `giant +/- lambda * baby`!
    // Let's test EXACTLY this!
    
    PointJacobian target_jac = PointJacobian::from_affine(Target); // Using real Target!
    PointAffine startG_aff = to_affine(scalar_mult(G, current_key));
    PointAffine neg_startG = startG_aff;
    neg_startG.y = mod_neg(neg_startG.y);
    
    PointJacobian point_jac = point_add_mixed(target_jac, neg_startG);
    
    PointAffine giant_step = to_affine(scalar_mult(G, uint256_t(m)));
    PointAffine neg_giant = giant_step;
    neg_giant.y = mod_neg(neg_giant.y);

    bool found = false;
    for(int step = 0; step < 10; step++) {
        PointAffine point = to_affine(point_jac);
        int64_t baby_idx = ht.lookup(point.x);
        if (baby_idx >= 0) {
            uint64_t raw_idx = (uint64_t)baby_idx;
            uint64_t endo_flag = (raw_idx >> 46) & 3;
            uint64_t actual_baby = raw_idx & ((1ULL << 46) - 1);
            
            // Wait, we inserted endo1.x meaning (lambda * actual_baby * G).x
            // And point = TargetE1 - giant_step
            // So TargetE1 - giant_step = lambda * actual_baby * G
            // TargetE1 = giant_step + lambda * actual_baby * G
            // But TargetE1 is lambda * Target
            // Therefore: lambda * Target = giant_step + lambda * actual_baby * G
            // Wait, our loop:
            // target_jac = TargetE1
            // point_jac = TargetE1 - startG
            // step 1: point = TargetE1 - startG - m*G
            // If point matches (lambda * actual_baby * G), then:
            // TargetE1 - startG - step*m*G = lambda * actual_baby * G
            // TargetE1 = startG + step*m*G + lambda * actual_baby * G
            // Wait, TargetE1 is lambda * Target for endo1...
            // Or is it?
            // Actually, in bsgs_engine.cpp, worker_sequential loops over:
            // point_jac = Target - startG - step*m*G
            // If point.x == endo1.x, then:
            // +/- (Target - startG - step*m*G) = lambda * actual_baby * G
            // So Target = startG + step*m*G +/- lambda * actual_baby * G
            // AHA! The target itself doesn't have lambda applied to it in the search loop!
            // The search loop just checks if (Target - Giant) matches Hash(lambda * baby)!
            
            uint256_t baby_key = uint256_t(actual_baby);
            if (endo_flag == 1) {
                baby_key = mod_n_mul(baby_key, get_lambda2());
            } else if (endo_flag == 2) {
                baby_key = mod_n_mul(baby_key, get_lambda());
            }
            
            // To test this exactly, I'll print the recovered key.
            uint256_t recovered = uint256_add(current_key, baby_key);
            
            printf("Found ENDO flag %llu at step %d! Recovered key: %llu, target: %llu\n", (unsigned long long)endo_flag, step, (unsigned long long)recovered.v[0], (unsigned long long)target_key.v[0]);
            
            if (recovered.v[0] == target_key.v[0]) {
                printf("SUCCESS!\n");
            }
            found = true;
            break;
        }
        
        point_jac = point_add_mixed(point_jac, neg_giant);
        current_key = uint256_add(current_key, uint256_t(m));
    }
    
    if(!found) {
        printf("Key not found in simulation!\n");
    }
}

int main() {
    test_endo_math();
    return 0;
}
