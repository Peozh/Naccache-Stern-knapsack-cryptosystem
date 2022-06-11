#pragma once
#include <utility>
#include <iostream>
#include <math.h>
#include <deque>
#include <vector>
#include <complex>
#include <algorithm>
#include <iomanip>

/*
    19 digit accurate
    PI = 3.141592653589793238'46L;
    we want ... 3.141592653589793238'46 264338327950288419716939937510...
    it's value .. 3.141592653589793238'51 280895940618620443274267017841339111328125
*/
constexpr long double PI = 3.141592653589793238'46L; 

/* 
    arr.size() = 2^k

    error of PI is ruining FFT. error is multiplied repeatedly.
    to solve this...
     - should use long double PI (significant figures : 19 digit)
     - use small valued argument (has smaller significant figures) : use only lower 16-bits of node.value
*/
void FFT_inplace(std::vector<std::complex<long double>>& arr, bool isInv) {
    for (int i = 0; i < arr.size(); ++i) { // reverse bit expression : to find swap position
        int j = 0;
        for (int r = 1, l = arr.size() >> 1; l > 0; r <<= 1, l >>= 1) {
            j |= i & l ? r : 0; 
        }
        if (i < j) swap(arr[i], arr[j]);
    }
    for (int i = 2; i <= arr.size(); i <<= 1) { // i = group size : 2 4 8 ...
        long double theta = 2.0L * (PI / i); 
        theta *= isInv ? -1.0L : 1.0L;
        for (int group_st = 0; group_st < arr.size(); group_st += i) { // group_st = offset of group
            std::complex<long double> w { 1.0L, 0.0L };
            for (int pair_st = 0; pair_st < i/2; ++pair_st) {
                w = { cos(theta * pair_st), sin(theta * pair_st) };
                std::complex<long double> right = arr[group_st + pair_st + i/2] * w; // odd * w
                arr[group_st + pair_st + i/2] = arr[group_st + pair_st] - right; // odd
                arr[group_st + pair_st] += right; // even
            }
        }
    }
    if (isInv) {
        for (int i = 0; i < arr.size(); ++i) {
            arr[i] /= arr.size();
        }
    }
};

// actually, it is 128-bit
struct int65_t {
    int64_t sign = 1LL;
    uint64_t value = 0ULL;
};

struct invCarry {
    int64_t sign = 1LL;
    uint64_t value = 0ULL;
    size_t radix_point = 0ULL; // 0~15
};


/* 
    Node
        : &=, |=, <<=, >>=, +=, -=

    sign is saved separately.
    use only lower 16-bits of node.value to save some value.
    upperbits is moved into carry.
*/
class Node {
public:
    int64_t sign = 1LL;
    uint64_t num = 0ULL;
    int65_t carry = { 1LL, 0ULL }; // sign, value
    invCarry invCarry = { 1LL, 0ULL, 0 }; // sign, value, radix_point
    Node() = default;
    ~Node() = default;
    Node(const Node& other) = default;
    Node(Node&& other) = default;
    Node& operator=(const Node& other) = default;
    Node& operator=(Node&& other) = default;

    Node(int64_t value) {
        if (value < 0LL) {
            sign = -1LL;
        }
        num = static_cast<uint64_t>(sign*value);
        updateCarry();
    }
    Node(int64_t sign, uint64_t num)
     : sign(sign), num(num)
    {}

    // |shift| < 16, carry updated, should manage carry
    Node& operator<<=(const size_t shift) { 
        num <<= shift;
        updateCarry();

        return *this;
    }
    // |shift| < 16, carry updated, should manage carry
    Node operator<<(const size_t shift) {
        Node result = *this;
        result <<= shift;

        return result;
    }

    // |shift| < 16, invCarry updated, should manage invCarry
    Node& operator>>=(const size_t shift) { 
        invCarry.sign = sign;
        invCarry.value = num & ((1ULL << shift) - 1);
        invCarry.radix_point = shift;
        num >>= shift;

        return *this;
    }
    // |shift| < 16, invCarry updated, should manage invCarry
    Node operator>>(const size_t shift) {
        Node result = *this;
        result >>= shift;

        return result;
    }

    // carry updated, should manage carry
    Node& operator+=(const Node& rhs) { // 0 인 경우 기존 부호 유지
        auto new_sign = num >= rhs.num ? sign : rhs.sign;
        num = new_sign*(sign*(int64_t)(num) + rhs.sign*(int64_t)(rhs.num));
        sign = new_sign;
        updateCarry();

        return *this;
    }
    // carry updated, should manage carry
    Node operator+(const Node& rhs) {
        Node result = *this;
        result += rhs;

        return result;
    }

    // carry updated, should manage carry
    Node& operator-=(const Node& rhs) {
        Node negated_rhs = rhs;
        negated_rhs.sign *= -1LL;
        *this += negated_rhs;

        return *this;
    }
    // carry updated, should manage carry
    Node operator-(const Node& rhs) {
        Node result = *this;
        result -= rhs;

        return result;
    }

    // no carry, don't care sign
    Node& operator|=(const Node& rhs) {
        num |= rhs.num;

        return *this;
    }
    Node operator|(const Node& rhs) {
        Node result = *this;
        result.num |= rhs.num;

        return result;
    }

    // no carry, don't care sign
    Node& operator&=(const Node& rhs) {
        num &= rhs.num;

        return *this;
    }
    Node operator&(const Node& rhs) {
        Node result = *this;
        result.num &= rhs.num;

        return result;
    }

    void updateCarry() {
        carry.value += num >> 16;
        carry.sign = sign;
        num &= 0xffffULL;
    }

    // high's carry updated, thus low to high propagation needed
    // sign can be different after operator- 
    // if (high.num == 0, carry.value == 0) carry.sign is dominant
    // we need sign for carry
    void applyCarry(Node& high) { 
        auto new_sign = (high.num <= carry.value) ? carry.sign : high.sign;
        high.num = new_sign * (high.sign * (int64_t)high.num + carry.sign * (int64_t)carry.value);
        high.sign = new_sign;
        carry.value = 0ULL;
        carry.sign = 1LL;
        high.updateCarry();
    }

    // for >> shift (radix point > 0)
    // for / remainder (radix point = 0)
    // applyCarry should be preceded. 
    // same sign in normal case. otherwise, invCarry is dominant.
    void applyInvCarry(Node& low) {
        auto new_sign = invCarry.sign;
        low.num = new_sign * (low.sign*(int64_t)(low.num) + invCarry.sign*(int64_t)(invCarry.value << (16 - invCarry.radix_point)));
        low.sign = new_sign;
        invCarry.sign = 1LL;
        invCarry.value = 0ULL;
        invCarry.radix_point = 0;
    }

    static std::complex<long double> toComplex(const Node& node) {
        return std::complex<long double> { (long double)node.sign * (int64_t)(node.num), 0.0L };
    }
};

/*
    Bits >= 64

    BigNumber
        : &=, |=, <<=, >>=, +=, -=, *=, >, <, ==, /=, %=, sqrt
        (inner implementation)
        : +    : { += }
        : >=   : { > }
        : *=   : { FFT }
        : /=   : { compare , -= , bit shift }
        : %=   : { compare , -= , bit shift } or { /=, *=, -= }
        : sqrt : { bits hift, -= }

    it can only comapre abs values
*/
template<size_t Bits = 64ULL>
class BigNumber {
private:
    
public:
    static constexpr size_t Count = Bits / 16ULL;
    std::deque<Node> container;

    BigNumber() 
     : container(Count, Node{ 1LL, 0ULL })
    {
    }

    BigNumber(int64_t num) 
     : container(Count, Node{ 1LL, 0ULL })
    {
        if (num < 0) {
            container[Count-4].sign = -1LL;
            container[Count-3].sign = -1LL;
            container[Count-2].sign = -1LL;
            container[Count-1].sign = -1LL;
            container[Count-4].num = ((uint64_t)(num * -1LL) >> 48) & 0xffffU;
            container[Count-3].num = ((uint64_t)(num * -1LL) >> 32) & 0xffffU;
            container[Count-2].num = ((uint64_t)(num * -1LL) >> 16) & 0xffffU;
            container[Count-1].num = (uint64_t)(num * -1LL) & 0xffffU;
            signCorrection();
        } else {
            container[Count-4].num = (num >> 48) & 0xffffU;
            container[Count-3].num = (num >> 32) & 0xffffU;
            container[Count-2].num = (num >> 16) & 0xffffU;
            container[Count-1].num = num & 0xffffU;
        }
    }

    BigNumber(uint64_t num, size_t shift)
     : container(Count, Node{ 1LL, 0ULL })
    {
        container[Count-4].num = (num >> 48) & 0xffffU;
        container[Count-3].num = (num >> 32) & 0xffffU;
        container[Count-2].num = (num >> 16) & 0xffffU;
        container[Count-1].num = num & 0xffffU;
        *this <<= shift;
    }

    ~BigNumber() = default;
    BigNumber(const BigNumber<Bits>& other) = default;
    BigNumber(BigNumber<Bits>&& other) = default;
    BigNumber& operator=(const BigNumber<Bits>& other) = default;
    BigNumber& operator=(BigNumber<Bits>&& other) = default;

    BigNumber<Bits>& operator<<=(size_t shift) { // O(Bits)
        size_t nodeShift = shift / 16ULL;
        shift %= 16ULL;
        int64_t sign = container[0].sign;
        for (size_t count = 0; count < nodeShift; ++count) {
            container.pop_front();
            container.push_back(Node { sign, 0ULL });
        }
        for (auto& node : container) {
            node <<= shift;
        }
        for (size_t idx = container.size()-1; idx >= 1; --idx) {
            container[idx].applyCarry(container[idx-1]);
        }
        if (container[0].carry.value != 0ULL) {
            std::cout << "overflowed : carry.value exists in maxNode" << std::endl;
        }
        //signCorrection();
        flowCheck();
        return *this;
    }
    BigNumber<Bits> operator<<(const size_t shift) {
        BigNumber<Bits> result = *this;
        result <<= shift;
        return result;
    }

    BigNumber<Bits>& operator>>=(size_t shift) { // O(Bits)
        size_t nodeShift = shift / 16ULL;
        shift %= 16ULL;
        int64_t sign = container[0].sign;
        for (size_t count = 0; count < nodeShift; ++count) {
            container.pop_back();
            container.push_front(Node { sign, 0ULL });
        }
        for (auto& node : container) {
            node >>= shift;
        }
        for (size_t idx = 0; idx < container.size()-1; ++idx) {
            container[idx].applyInvCarry(container[idx+1]);
        }
        //signCorrection();
        return *this;
    }
    BigNumber<Bits> operator>>(const size_t shift) {
        BigNumber<Bits> result = *this;
        result >>= shift;
        return result;
    }

    BigNumber<Bits>& operator|=(const BigNumber<Bits>& rhs) { // O(Bits)
        for (size_t idx = 0; idx < container.size(); ++idx) {
            container[idx] |= rhs.container[idx];
        }
        return *this;
    }
    BigNumber<Bits> operator|(const BigNumber<Bits>& rhs) {
        BigNumber<Bits> result = *this;
        result |= rhs;
        return result;
    }
    
    BigNumber<Bits>& operator&=(const BigNumber<Bits>& rhs) { // O(Bits)
        for (size_t idx = 0; idx < container.size(); ++idx) {
            container[idx] &= rhs.container[idx];
        }
        return *this;
    }
    BigNumber<Bits> operator&(const BigNumber<Bits>& rhs) {
        BigNumber<Bits> result = *this;
        result &= rhs;
        return result;
    }

    BigNumber<Bits>& operator+=(const BigNumber<Bits>& rhs) { // O(Bits)
        for (size_t idx = 0; idx < container.size(); ++idx) {
            container[idx] += rhs.container[idx];
        }
        for (size_t idx = container.size() - 1; idx >= 1; --idx) {
            container[idx].applyCarry(container[idx-1]);
        }
        signCorrection();
        flowCheck();
        return *this;
    }
    BigNumber<Bits> operator+(const BigNumber<Bits>& rhs) {
        BigNumber<Bits> result = *this;
        result += rhs;
        return result;
    }

    BigNumber<Bits>& operator-=(const BigNumber<Bits>& rhs) { // O(Bits)
        for (size_t idx = 0; idx < container.size(); ++idx) {
            container[idx] -= rhs.container[idx];
        }
        for (size_t idx = container.size() - 1; idx >= 1; --idx) {
            container[idx].applyCarry(container[idx-1]); // 왼쪽으로 쓸어올리기
        }
        signCorrection(); // 오른쪽으로 쓸어내리기
        return *this;
    }
    BigNumber<Bits> operator-(const BigNumber<Bits>& rhs) {
        BigNumber<Bits> result = *this;
        result -= rhs;
        return result;
    }

    BigNumber<Bits>& operator*=(const BigNumber<Bits>& rhs) { // O(Bits*log(Bits))
        std::vector<std::complex<long double>> A(2 * container.size());
        std::transform(container.begin(), container.end(), A.begin() + container.size(), Node::toComplex);
        std::vector<std::complex<long double>> B(2 * container.size());
        std::transform(rhs.container.begin(), rhs.container.end(), B.begin() + container.size(), Node::toComplex);

        FFT_inplace(A, false);
        FFT_inplace(B, false);

        std::vector<std::complex<long double>> C(2 * container.size());
        for (int i = 0; i < C.size(); ++i) {
            C[i] = A[i] * B[i];
        }

        FFT_inplace(C, true);

        // std::cout.precision(3);
        // std::cout << std::scientific;
        // std::cout << "complex C<" << Bits * 2 << "> : \t\t";
        // for (auto c : C) {
        //     std::cout.width(10);
        //     std::cout << c.real() << ' ';
        // }
        // std::cout << std::endl;

        for (size_t idx = 0; idx < container.size(); ++idx) {
            int64_t sign = C[idx + container.size()-1].real() > 0.0L ? 1LL : -1LL;
            container[idx] = Node { sign, (uint64_t)(round(abs(C[idx + container.size()-1].real()))) };
            //container[idx].updateCarry();
        }
        for (size_t idx = container.size() - 1; idx >= 1; --idx) {
            //container[idx].applyCarry(container[idx-1]); // 왼쪽으로 쓸어올리기
        }
        signCorrection(); // 오른쪽으로 쓸어내리기
        container[container.size() - 1].updateCarry();
        for (size_t idx = container.size() - 1; idx >= 1; --idx) {
            container[idx].applyCarry(container[idx-1]); // 왼쪽으로 쓸어올리기
        }
        flowCheck();
        return *this;
    }
    BigNumber<Bits> operator*(const BigNumber<Bits>& rhs) {
        BigNumber<Bits> result = *this;
        result *= rhs;

        return result;
    }

    bool operator>(const BigNumber<Bits>& rhs) const { // O(Bits)
        for (size_t idx = 0; idx < container.size(); ++idx) {
            if (container[idx].num > rhs.container[idx].num) return true;
            else if (container[idx].num < rhs.container[idx].num) return false;
        }
        return false;
    }
    bool operator>=(const BigNumber<Bits>& rhs) const {
        for (size_t idx = 0; idx < container.size(); ++idx) {
            if (container[idx].num > rhs.container[idx].num) return true;
            else if (container[idx].num < rhs.container[idx].num) return false;
        }
        return true;
    }
    bool operator<(const BigNumber<Bits>& rhs) const {
        return !(*this >= rhs);
    }
    bool operator<=(const BigNumber<Bits>& rhs) const {
        return !(*this > rhs);
    }
    bool operator==(const BigNumber<Bits>& rhs) const {
        for (size_t idx = 0; idx < container.size(); ++idx) {
            if (container[idx].num != rhs.container[idx].num) return false;
        }
        return true;
    }
    bool operator!=(const BigNumber<Bits>& rhs) const {
        return !(*this == rhs);
    }

    BigNumber<Bits>& operator/=(BigNumber<Bits> rhs) { // O(Bits^2)
        BigNumber<Bits> result;
        int64_t new_sign = this->container[0].sign * rhs.container[0].sign;
        *this *= BigNumber<Bits>(this->container[0].sign);
        rhs *= BigNumber<Bits>(rhs.container[0].sign);

        if (rhs == BigNumber<Bits>()) {
            std::cout << "divided by 0 rhs : BigNumber<Bits>::operator/=" << std::endl;
        } else {
            if (*this >= rhs) {
                size_t maxShift = 0ULL;
                while (*this >= (rhs << 1)) { // *this >= rhs
                    rhs <<= 1;
                    ++maxShift;
                }
                *this -= rhs;
                result += BigNumber<Bits>(1ULL, maxShift);
                
                for (size_t shift = 1; shift <= maxShift; ++shift) {
                    rhs >>= 1;
                    if (*this >= rhs) {
                        *this -= rhs;
                        result += BigNumber<Bits>(1ULL, maxShift-shift);
                    }
                }
            }
        }
        *this = result;
        *this *= BigNumber<Bits>(new_sign);
        return *this;
    }
    BigNumber<Bits> operator/(BigNumber<Bits> rhs) {
        BigNumber<Bits> result = *this;
        result /= rhs;
        return result;
    }
    
    BigNumber<Bits>& operator%=(BigNumber<Bits> rhs) { // O(Bits^2)
        if (rhs == BigNumber<Bits>()) {
            std::cout << "divided by 0 rhs : BigNumber<Bits>::operator%=" << std::endl;
        } else {
            if (*this >= rhs) {
                size_t maxShift = 0;
                while (*this >= (rhs << 1ULL)) { // *this >= rhs
                    rhs <<= 1ULL;
                    ++maxShift;
                }
                *this -= rhs;
                
                for (size_t shift = 1; shift <= maxShift; ++shift) {
                    rhs >>= 1ULL;
                    if (*this >= rhs) {
                        *this -= rhs;
                    }
                }
            }
        }
        return *this;
    }
    BigNumber<Bits> operator%(BigNumber<Bits> rhs) {
        BigNumber<Bits> result = *this;
        result %= rhs;
        return result;
    }

    void signCorrection() { // O(Bits)
        int firstNonZeroIdx = -1;
        for (auto node : container) { // find valid node
            ++firstNonZeroIdx;
            if (node.num != 0) break;
        }
        for (size_t idx = 0; idx < firstNonZeroIdx; ++idx) { // make (front non-valid(=0) bits's sign) to (valid bit's sign)
            container[idx].sign = container[firstNonZeroIdx].sign;
        }
        for (size_t idx = firstNonZeroIdx; idx < container.size()-1; ++idx) { // if sign is flipped, generate invCarry
            if (container[idx].sign * container[idx+1].sign < 0LL) {
                container[idx].invCarry = { container[idx].sign, 1ULL, 0 };
                container[idx] -= Node { container[idx].sign, 1ULL };
                container[idx].applyInvCarry(container[idx+1]);
            }
        }
    }

    void print(bool numOnly) {
        for (const auto& node : container) {
            if (node.sign == -1) std::cout << " -";
            std::cout << " " << node.num;
        }
        std::cout << std::endl;
        if (numOnly == true) return;
        for (const auto& node : container) {
            std::cout << " " << node.carry.value;
        }
        std::cout << std::endl;
        for (const auto& node : container) {
            std::cout << " " << node.invCarry.value;
        }
        std::cout << std::endl;
    }

    void flowCheck() {
        if (container[0].carry.value != 0ULL) {
            if (container[0].carry.sign > 0LL) std::cout << "overflow" << std::endl;
            else std::cout << "underflow" << std::endl;
        }
    }

    static BigNumber<Bits> sqrt(const BigNumber<Bits>& rhs) { // O(Bits)
        BigNumber<Bits> result;
        BigNumber<Bits> tool;
        BigNumber<Bits> target;
        for (auto node : rhs.container) {
            if (node.num  > 0) {
                int a = 0;
            }
            for (size_t shift = 2; shift <= 16; shift += 2) {
                if (shift == 16) {
                    int b = 0;
                }
                target <<= 2;
                target.container.back().num |= ((node.num >> (16 - shift)) & 0b11ULL);

                tool.container.back().num += tool.container.back().num & 0b1ULL;
                tool <<= 1;
                tool.container.back().num |= 0b1ULL;
                
                result <<= 1;
                result.container.back().num |= 0b1ULL;
                
                if (target < tool) {
                    tool.container.back().num &= 0xfffeULL;
                    result.container.back().num &= 0xfffeULL;
                } else {
                    target -= tool;
                }
            }
        }

        return result;
    }
    
    template <size_t Bits2>
    friend std::ostream& operator<<(std::ostream& os, const BigNumber<Bits2>& bignum);
};

template <size_t Bits>
std::ostream& operator<<(std::ostream& os, const BigNumber<Bits>& bignum)
{
    for (const auto& node : bignum.container) {
        if (node.sign == -1) os << " -";
        os << " " << node.num;
    }
    os << " = " << bignum.container[bignum.container.size()-1].num
     + 0x10000ULL * bignum.container[bignum.container.size()-2].num
     + 0x100000000ULL * bignum.container[bignum.container.size()-3].num
     + 0x1000000000000ULL * bignum.container[bignum.container.size()-4].num;
    return os;
};

