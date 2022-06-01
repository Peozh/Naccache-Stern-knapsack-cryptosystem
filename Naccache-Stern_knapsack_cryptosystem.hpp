#include <iostream>
#include <vector>
#include <random>

/*  
    Naccache–Stern knapsack cryptosystem
    https://en.wikipedia.org/wiki/Naccache%E2%80%93Stern_knapsack_cryptosystem
    https://www.di.ens.fr/~stern/data/St63.pdf

    Pohlig–Hellman algorithm
    https://en.wikipedia.org/wiki/Pohlig%E2%80%93Hellman_algorithm
    https://mathphysics.tistory.com/1148
*/

// N : unit length, bit count of unit data.
template <std::size_t N = 8ULL>
class KnapsackCryptosystem {
public:
    uint64_t largePrime_p;
    std::vector<uint64_t> coprimes_vis;

    KnapsackCryptosystem() {
        init();
    }

    uint64_t encrypt(uint8_t message) {
        uint64_t activeBit = 1;
        uint64_t encrypted_c = 1;
        for (std::size_t idx = 0; idx < N; ++idx) {
            if (activeBit & message) {
                encrypted_c *= coprimes_vis[idx];
                encrypted_c %= largePrime_p;
            }
            activeBit <<= 1;
        }
        return encrypted_c;
    }

    uint8_t decrypt(uint64_t encrypted_c) {
        uint64_t data = 1;
        uint64_t activeBit = 1ULL << 63;
        for (std::size_t idx = 0; idx < 64; ++idx) {
            data *= data;
            data %= largePrime_p;
            if (activeBit & secretKey_s) {
                data *= encrypted_c;
                data %= largePrime_p;
            }
            activeBit >>= 1;
        }
        
        uint8_t message = 0;
        uint8_t activeBit2 = 1;
        for (std::size_t idx = 0; idx < N; ++idx) {
            if (data % primes[idx] == 0) {
                message |= activeBit2;
            }
            activeBit2 <<= 1;
        }
        return message;
    }
    
    // Pohlig-Hellman 
    // https://www.youtube.com/watch?v=B0p0jbCGvWk
    // p_0 = v_0^s (mod p)
    // suppose we know p_0 = 2
    uint64_t sFinder() {
        // find factors of p-1. factors = { (p_i)^e_i } = { factor_i }
        // p-1 = 3037000492 = 2^2 * 1543 * 492061
        std::vector<uint64_t> factors {4ULL, 1543ULL, 492061ULL};
        std::vector<uint64_t> ris {0ULL, 0ULL, 0ULL};
        std::vector<uint64_t> exponents {0ULL, 0ULL, 0ULL};
        auto b = coprimes_vis[0];
        std::uint64_t a = 2ULL;
        // a = b^s (mod p). find s
        for (std::size_t i = 0; i < factors.size(); ++i) {
            /*  
                (q stands for quotient, r stands for remainder)
                s = (factor_i)q_i + r_i
                a = b^s (mod p)
                a^((p-1)/factor_i) = b^(r_i*(p-1)/factor_i) (mod p)
                
                lhs = a^((p-1)/factor_i) (mod p)
                rhs = b^(r_i*(p-1)/factor_i) (mod p)
            */
            std::uint64_t exponent_i = (largePrime_p-1)/factors[i];
            exponents[i] = exponent_i;
            std::uint64_t lhs = modular_exponentiation(a, exponent_i, largePrime_p);
            /*
                R_i = { 0 ~ factor_i-1 }, r_i ∈ R_i
                s = r_i (mod factor_i)
            */
            for (std::uint64_t r_i = 0; r_i < factors[i]; ++r_i) {
                std::uint64_t rhs = modular_exponentiation(b, (r_i*exponent_i), largePrime_p);
                if (lhs == rhs) {
                    ris[i] = r_i;
                    break;
                }
            }
        }
        /*
            for all r_i, s = r_i (mod factor_i) is true.

            Chinese remainder theorem : s exists uniquly in range [0, p-1).
                s = (M_i*exponent_i + m_i*factor_i) * r_i (mod factor_i) --- (M_i, m_i is some coefficient)
            Bézout's identity
             : when gcd(exponent_i, factor_i) = 1 (coprime condition),
             : M_i, m_i pair that satisfying below exist.
                M_i*exponent_i + m_i*factor_i = 1

            s = Sigma{ (M_i*exponent_i) * r_i }

            s = (M_i*exponent_i) * r_i (mod factor_i)
            M_i*exponent_i = 1 (mod factor_i)
            M_i = 1/exponent_i (mod factor_i)

            s = Sigma{ (1/exponent_i*exponent_i) * r_i }
        */
        
        for (auto ri : ris) {
            std::cout << "ri : " << ri << std::endl;
        }
        uint64_t deduced_s = 0ULL;
        for (size_t i = 0; i < ris.size(); ++i) {
           auto M_i = modular_multiplicative_inverse(factors[i], exponents[i]);
           uint64_t sigma_i = M_i;
           sigma_i *= exponents[i];
           sigma_i %= largePrime_p-1;
           sigma_i *= ris[i];
           sigma_i %= largePrime_p-1;

           deduced_s += sigma_i;
           deduced_s %= largePrime_p-1;
        }
        // s = deduced_s (mod p-1)
        return deduced_s;
    }


private:
    
    void init() {
        setLargePrime();
        setSecretKey();
        setCoprimes();
    }

    bool isPrime(uint64_t num) {
        uint64_t bound = sqrtl(num);
        for (uint64_t factor = 3; factor <= bound; factor += 2) {
            if (num % factor == 0) {
                return false;
            }
        }
        std::cout << "found prime = " << num << std::endl;
        return true;
    }

    // p > p_0 * ... * p_(N-1)
    void setLargePrime() {
        primes = { 2,  3,  5,  7, 11, 13, 17, 19, // 9,699,690
                  23, 29, 31, 37, 41, 43, 47}; // 614,889,782,588,491,410
        largePrime_p = 9700247ULL; // max 18,446,744,073,709,551,616(2^64) - 59 with uint64_t
        largePrime_p = 18'446'744'073'709'551'557ULL; 
        // 18,446,744,073,709,551,616 : 2^64 (uint64_t)
        //  9 223 372 036 854 775 808 : 2^63 (int64_t)
        //  9 223 372 035 854 775 787 : biggest prime in 2^63 (int64_t)
        //              3 037 000 499 : sqrt(2^63) (int64_t)
        //              3 037 000 493 : biggest prime in sqrt(2^63) (int64_t)

        uint64_t primeUpperbound = sqrtl(9'223'372'035'854'775'808ULL);
        std::size_t maximumPrimeGap = sqrtl(primeUpperbound); // Oppermann's conjecture : (n^2-n, n^2)
        std::cout << "maximum prime gap = " << maximumPrimeGap << std::endl;
        for (std::size_t idx = 0; idx < maximumPrimeGap; idx += 2) { // density ~= 1/log(n) = 4.58%
            if (isPrime(primeUpperbound - idx)) {
                largePrime_p = primeUpperbound - idx;
                break;
            }
        }

        largePrime_p = 3'037'000'493ULL; // p-1 = 2^2 * 1543 * 492061
    }

    /* 
        gcd(s, p-1) = 1
        it guarantees existance of MMI (s^-1 modulo p-1)
    */
    void setSecretKey() {
        //secretKey_s = 5642069;
        std::random_device seed;
        std::mt19937_64 randomEngine(seed());
        std::uniform_int_distribution<uint64_t> distribution(2ULL, largePrime_p-2);
        secretKey_s = distribution(randomEngine);
        while(gcd(secretKey_s, largePrime_p - 1) != 1) { // Probability of coprimality : 60.8%
            std::cout << "retrying s ... : " << secretKey_s << " -> ";
            secretKey_s = distribution(randomEngine);
            std::cout << secretKey_s << std::endl;
        }
        std::cout << "confirmed s = " << secretKey_s << std::endl;
    }

    // Greatest Common Devisor
    uint64_t gcd(uint64_t x, uint64_t y) {
        while (x != 0) {
            y = y % x;
            std::swap(x, y);
        }
        return y;
    }

    // Set public keys v_i (= encrypted p_i) as coprime with p
    void setCoprimes() {
        //coprimes_vis = {8567078, 5509479, 2006538, 4340987, 8643477, 6404090, 1424105, 7671241};
        coprimes_vis.resize(N);
        for (std::size_t idx = 0; idx < coprimes_vis.size(); ++idx) {
            coprimes_vis[idx] = s_th_root_modulo_p(primes[idx]);
        }
    }

    /*  
        (v_i)^s = p_i mod p
        v_i = s-th_root(p_i) mod p
        v_i = p_i ^ (s^-1) mod p
        (Fermat's little Theorem, FlT) applied... guaranteed by (p is prime)
        v_i = p_i ^ (s^-1 mod p-1) mod p
        (Modular Multiplicative Inverse, MMI) applied... guaranteed by (gcd(s, p-1) = 1)
        v_i = p_i ^ (MMI(p-1, s)) mod p
    */
    uint64_t s_th_root_modulo_p(uint64_t p_i) {
        uint64_t v_i = 1;
        uint64_t exponent = modular_multiplicative_inverse(largePrime_p-1, secretKey_s);
        uint64_t activeBit = 1ULL << 63;
        for (std::size_t idx = 0; idx < 64; ++idx) {
            v_i *= v_i; // can overflow if (largePrime_p > sqrt(2^63) > 3037000493)
            v_i %= largePrime_p;
            if (activeBit & exponent) {
                v_i *= p_i;
                v_i %= largePrime_p;
            }
            activeBit >>= 1;
        }
        return v_i;
    }

    /*  
        Modular Multiplicative Inverse w\ Extended Euclidean Algorithm
        a*x = 1 mod modulus
        return a = x^-1 mod modulus
    */
    uint64_t modular_multiplicative_inverse(int64_t modulus, int64_t x) {
        int64_t k[2] = {1, 0};
        int64_t mmi[2] = {0, 1};
        int64_t right[2] = {modulus, x};
        while (right[1] != 0) {
            int64_t quotient = right[0] / right[1];
            k[0] -= quotient * k[1];
            mmi[0] -= quotient * mmi[1];
            right[0] -= quotient * right[1];
            std::swap(k[0], k[1]);
            std::swap(mmi[0], mmi[1]);
            std::swap(right[0], right[1]);
        }
        return (mmi[0] + modulus) % modulus;
    }

    uint64_t modular_exponentiation(uint64_t base, uint64_t exponent, uint64_t modulus) {
        uint64_t result = 1ULL;
        uint64_t activeBit = 1ULL << 63;
        for (size_t shift = 0; shift < 64; ++shift) {
            result *= result;
            result %= modulus;
            if (activeBit & exponent) {
                result *= base;
                result %= modulus;
            }
            activeBit >>= 1;
        }
        return result;
    }

private:
    std::vector<uint64_t> primes;
    uint64_t secretKey_s;
};