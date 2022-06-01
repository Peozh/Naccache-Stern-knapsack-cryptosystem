#include "Naccache-Stern_knapsack_cryptosystem.hpp"


int main() {
    KnapsackCryptosystem<8> kc;

    for (auto v : kc.coprimes_vis) {
        std::cout << "vi : " << v << std::endl;
    }

    uint8_t message = 202; // 1100 1010
    std::cout << "message : " << static_cast<int>(message) << std::endl;
    uint64_t encrypted_c = kc.encrypt(message);
    uint8_t decrypted_message = kc.decrypt(encrypted_c);

    std::cout << "encrypted message : " << encrypted_c << std::endl;
    std::cout << "decrypted message : " << static_cast<int>(decrypted_message) << std::endl;
    
    auto deduced_s = kc.sFinder();
    std::cout << "deduced s = " << deduced_s << std::endl;
    return 0;
}