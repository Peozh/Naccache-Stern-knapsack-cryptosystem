#pragma once
#include <iostream>
#include <vector>
#include <random>
#include "BigNumberUnsigned64_16.hpp"

/*  
    Naccache–Stern knapsack cryptosystem
    https://en.wikipedia.org/wiki/Naccache%E2%80%93Stern_knapsack_cryptosystem
    https://www.di.ens.fr/~stern/data/St63.pdf

    Pohlig–Hellman algorithm
    https://en.wikipedia.org/wiki/Pohlig%E2%80%93Hellman_algorithm
    https://mathphysics.tistory.com/1148
*/

/*
    N : unit's length, bit count of unit data.
    (N, Bits, primeMaker) = (16, 256, 59), (32, 512, 229), (64, 1024, 331)
*/
template <size_t N = 64ULL, size_t Bits = N*16ULL>
class NaccacheSternKnapsackCryptosystem {
public:
    static constexpr size_t Bit = Bits;
    BigNumber<Bits> largePrime_p;
    std::vector<BigNumber<Bits>> coprimes_vis;
    BigNumber<Bits> primeMaker;

    NaccacheSternKnapsackCryptosystem() {
    }

    void setPrimeMaker(uint64_t pm) {
        primeMaker = BigNumber<Bits>(pm, 0);
    }

    void init() {
        setPrimes(N);
        //setLargePrime(); // it takes too much time(isPrime function is the reason). changed to set (largePrime_p = mult + primeMaker) in setPrimes()
        setSecretKey();
        setCoprimes();
    }

    // init with text file
    void init (std::ifstream& ifs) {
        setPrimes(N);
        for (auto& node : secretKey_s.container) {
            ifs >> node.num;
        }
        coprimes_vis.resize(N);
        for (auto& vi : coprimes_vis) {
            for (auto& node : vi.container) {
                ifs >> node.num;
            }
        }
    }

    //todo : save it with binary or text file
    void saveStates(std::ofstream& ofs) {
        for (auto& node : secretKey_s.container) {
            ofs << node.num;
        }
        ofs << std::endl;
        for (auto& vi : coprimes_vis) {
            for (auto& node : vi.container) {
                ofs << node.num;
            }
            ofs << std::endl;
        }
        ofs << std::endl;
    }

    //todo : secret key regeneration
    


    BigNumber<Bits> encrypt(uint64_t message) {
        auto zero = BigNumber<Bits>(0);
        BigNumber<Bits> activeBit = BigNumber<Bits>(1);
        BigNumber<Bits> encrypted_c = BigNumber<Bits>(1);
        auto msg = BigNumber<Bits>(message, 0);
        
        for (std::size_t idx = 0; idx < N; ++idx) {
            if ((activeBit & msg) > zero) {
                encrypted_c *= coprimes_vis[idx];
                encrypted_c %= largePrime_p;
            }
            activeBit <<= 1;
        }
        return encrypted_c;
    }

    // gcd(p_i, c^s) == p_i or 1 ?
    uint64_t decrypt(BigNumber<Bits> encrypted_c) {
        auto zero = BigNumber<Bits>(0);
        BigNumber<Bits> data = BigNumber<Bits>(1);
        BigNumber<Bits> activeBit = BigNumber<Bits>(1ULL, Bits/2 - 1);
        for (std::size_t idx = Bits/2; idx < Bits; ++idx) { // compute data = c^s
            data *= data;
            data %= largePrime_p;
            if ((activeBit & secretKey_s) > zero) {
                data *= encrypted_c;
                data %= largePrime_p;
            }
            activeBit >>= 1;
        }

        uint64_t message = 0;
        uint64_t activeBit2 = 1;
        for (std::size_t idx = 0; idx < N; ++idx) {
            if ((data % primes[idx]) == zero) {
                message |= activeBit2;
            }
            activeBit2 <<= 1;
        }
        return message;
    }
    
private:

    void setPrimes(size_t testN) {
        primes.reserve(testN);
        primes.push_back(BigNumber<Bits>(2ULL, 0));
        BigNumber<Bits> testingNumber = BigNumber<Bits>(3ULL, 0);
        BigNumber<Bits> step = BigNumber<Bits>(2ULL, 0);

        BigNumber<Bits> mult = BigNumber<Bits>(1ULL, 0);
        while (primes.size() < testN) {
            if (isPrime(testingNumber)) {
                primes.push_back(testingNumber);
            }
            testingNumber += step;
        }
        for (auto prime : primes) {
            mult *= prime;
        }

        std::cout << "mult = " << mult << std::endl;
        largePrime_p = mult + primeMaker;
        std::cout << "largePrime_p = " << largePrime_p << std::endl;
    }

    
    bool isPrime(BigNumber<Bits> num) {
        BigNumber<Bits> bound = BigNumber<Bits>::sqrt(num);
        if ((num % BigNumber<Bits>(2)) == BigNumber<Bits>(0)) return false;
        for (BigNumber<Bits> factor = BigNumber<Bits>(3); factor <= bound; factor += BigNumber<Bits>(2)) {
            if ((num % factor) == BigNumber<Bits>(0)) {
                return false;
            }
        }
        
        return true;
    }

    // p > p_0 * ... * p_(N-1)
    void setLargePrime() {
        BigNumber<Bits> primeUpperbound = BigNumber<Bits>::sqrt(BigNumber<Bits>(1ULL, Bits/2-4)) - BigNumber<Bits>(1ULL, 0);
        BigNumber<Bits> maximumPrimeGap = BigNumber<Bits>::sqrt(primeUpperbound); // Oppermann's conjecture : (n^2-n, n^2)

        for (BigNumber<Bits> idx = BigNumber<Bits>(0); idx < maximumPrimeGap; idx += BigNumber<Bits>(2)) { // density ~= 1/log(mult) = 0.789% ~ 4.58%
            if (isPrime(primeUpperbound - idx)) {
                largePrime_p = primeUpperbound - idx;
                break;
            }
        }

        std::cout << "largePrime_p = " << largePrime_p << std::endl;
    }

    /* 
        gcd(s, p-1) = 1
        it guarantees existance of MMI (s^-1 modulo p-1)
    */
    void setSecretKey() {
        std::random_device seed;
        std::mt19937_64 randomEngine(seed());

        BigNumber<Bits> upperBound = largePrime_p - BigNumber<Bits>(4);

        for (size_t idx = 0; idx < upperBound.container.size(); ++idx) {
            std::uniform_int_distribution<uint64_t> distribution(0ULL, upperBound.container[idx].num);
            secretKey_s.container[idx].num = distribution(randomEngine);
        }
        secretKey_s += BigNumber<Bits>(2);
        while(gcd(secretKey_s, largePrime_p - BigNumber<Bits>(1)) != BigNumber<Bits>(1)) { // Probability of coprimality : 60.8%
            std::cout << "retrying s ... : " << secretKey_s << " -> ";
            bool isSmallS = false;
            for (size_t idx = 0; idx < upperBound.container.size(); ++idx) {
                std::uniform_int_distribution<uint64_t> distribution(0ULL, isSmallS ? 65535ULL : upperBound.container[idx].num);
                secretKey_s.container[idx].num = distribution(randomEngine);
                if (!isSmallS) isSmallS = upperBound.container[idx].num > secretKey_s.container[idx].num;
            }
            secretKey_s += BigNumber<Bits>(2);

            std::cout << secretKey_s << std::endl;
        }
        std::cout << "confirmed s = " << secretKey_s << std::endl;
    }

    // Greatest Common Devisor
    BigNumber<Bits> gcd(BigNumber<Bits> x, BigNumber<Bits> y) {
        while (x != BigNumber<Bits>(0)) {
            y %= x;
            std::swap(x, y);
        }
        return y;
    }

    // Set public keys v_i (= encrypted p_i with s)
    void setCoprimes() {
        coprimes_vis.resize(N);
        auto exponent = modular_multiplicative_inverse(largePrime_p - BigNumber<Bits>(1), secretKey_s);
        for (std::size_t idx = 0; idx < coprimes_vis.size(); ++idx) {
            /*  
                (v_i)^s = p_i mod p
                v_i = s-th_root(p_i) mod p
                v_i = p_i ^ (s^-1) mod p
                (Fermat's little Theorem, FlT) applied... guaranteed by (p is prime)
                v_i = p_i ^ (s^-1 mod p-1) mod p
                (Modular Multiplicative Inverse, MMI) applied... guaranteed by (gcd(s, p-1) = 1) : element's minimum gap = gcd(s, p-1)
                v_i = p_i ^ (MMI(p-1, s)) mod p
            */
            coprimes_vis[idx] = modular_exponentiation(primes[idx], exponent, largePrime_p);
        }
    }

    /*  
        Modular Multiplicative Inverse w\ Extended Euclidean Algorithm
        a*x = 1 mod modulus
        return a = x^-1 mod modulus
    */
    BigNumber<Bits> modular_multiplicative_inverse(BigNumber<Bits> modulus, BigNumber<Bits> x) {
        BigNumber<Bits> k[2] = {BigNumber<Bits>(1), BigNumber<Bits>(0)};
        BigNumber<Bits> mmi[2] = {BigNumber<Bits>(0), BigNumber<Bits>(1)};
        BigNumber<Bits> right[2] = {modulus, x};
        while (right[1] != BigNumber<Bits>(0)) {
            BigNumber<Bits> quotient = right[0] / right[1];
            k[0] -= quotient * k[1];
            mmi[0] -= quotient * mmi[1];
            right[0] -= quotient * right[1];
            std::swap(k[0], k[1]);
            std::swap(mmi[0], mmi[1]);
            std::swap(right[0], right[1]);
        }
        return (mmi[0] + modulus) % modulus;
    }

    // result = base^exponent (mod modulus)
    BigNumber<Bits> modular_exponentiation(BigNumber<Bits> base, BigNumber<Bits> exponent, BigNumber<Bits> modulus) { // O(Bits*{ O(*=) + O(%=) })
        BigNumber<Bits> result = BigNumber<Bits>(1ULL);
        BigNumber<Bits> activeBit = BigNumber<Bits>(1ULL, Bits-1);
        
        for (size_t shift = 0; shift < Bits; ++shift) {
            result *= result;
            result %= modulus;
            if ((activeBit & exponent) != BigNumber<Bits>(0)) {
                result *= base;
                result %= modulus;
            }
            activeBit >>= 1;
        }
        return result;
    }

private:
    std::vector<BigNumber<Bits>> primes;
    BigNumber<Bits> secretKey_s;
};