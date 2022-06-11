#include <fstream>
#include "Naccache-Stern_knapsack_cryptosystem.hpp"
#include "Naccache-Stern_knapsack_cryptosystem_BigNum.hpp"

int main() {
    // Naccache-Stern_knapsack_cryptosystem.hpp  with  N = 8, uint64_t
    if (false) {
        KnapsackCryptosystem<8> kc;

        uint8_t message = 202; // 1100 1010
        std::cout << "message : " << static_cast<int>(message) << std::endl;

        uint64_t encrypted_c = kc.encrypt(message);
        std::cout << "encrypted message : " << encrypted_c << std::endl;

        uint8_t decrypted_message = kc.decrypt(encrypted_c);
        std::cout << "decrypted message : " << static_cast<int>(decrypted_message) << std::endl;
        
        auto deduced_s = kc.sFinder();
        std::cout << "deduced s = " << deduced_s << std::endl;
    }

    /*
        Naccache-Stern_knapsack_cryptosystem_BigNum.hpp 
        (N, Bits, primeMaker) = (64, 1024, 331)
        text file init
    */
    if (false) { 
        NaccacheSternKnapsackCryptosystem<64> kc;
        kc.setPrimeMaker(331ULL);

        std::ifstream ifs("init.txt");
        kc.init(ifs);
        ifs.close();

        uint64_t message = 0ULL-1ULL;
        std::cout << "message : " << static_cast<uint64_t>(message) << std::endl;

        BigNumber<NaccacheSternKnapsackCryptosystem<64>::Bit> encrypted_c = kc.encrypt(message);
        std::cout << "encrypted message : " << encrypted_c << std::endl;

        uint64_t decrypted_message = kc.decrypt(encrypted_c);
        std::cout << "decrypted message : " << static_cast<uint64_t>(decrypted_message) << std::endl;
    }

    /*
        Naccache-Stern_knapsack_cryptosystem_BigNum.hpp
        (N, Bits, primeMaker) = (32, 512, 229)
    */
    if (false) { 
        NaccacheSternKnapsackCryptosystem<32> kc;
        kc.setPrimeMaker(229ULL);
        kc.init();

        uint32_t message = 0UL-1UL;
        std::cout << "message : " << static_cast<uint64_t>(message) << std::endl;

        BigNumber<NaccacheSternKnapsackCryptosystem<32>::Bit> encrypted_c = kc.encrypt(message);
        std::cout << "encrypted message : " << encrypted_c << std::endl;

        uint64_t decrypted_message = kc.decrypt(encrypted_c);
        std::cout << "decrypted message : " << static_cast<uint64_t>(decrypted_message) << std::endl;
    }

    /*
        Naccache-Stern_knapsack_cryptosystem_BigNum.hpp
        (N, Bits, primeMaker) = (16, 256, 59)
    */
    if (false) { 
        NaccacheSternKnapsackCryptosystem<16> kc;
        kc.setPrimeMaker(59ULL);
        kc.init();

        uint16_t message = -1;
        std::cout << "message : " << static_cast<uint64_t>(message) << std::endl;

        BigNumber<NaccacheSternKnapsackCryptosystem<16>::Bit> encrypted_c = kc.encrypt(message);
        std::cout << "encrypted message : " << encrypted_c << std::endl;

        uint64_t decrypted_message = kc.decrypt(encrypted_c);
        std::cout << "decrypted message : " << static_cast<uint64_t>(decrypted_message) << std::endl;
    }
    return 0;
}