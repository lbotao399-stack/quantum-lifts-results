// Exact four-flavour exterior coefficients, using sparse matchings of the native quartic.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#include <vector>

using Wide = __int128_t;
using Entry = std::array<int,4>;

int parity(int mask) { return __builtin_popcount(static_cast<unsigned>(mask)) % 2; }

std::string decimal(Wide n) {
    if (n == 0) return "0";
    bool negative = n < 0;
    if (negative) n = -n;
    std::string text;
    while (n) {
        text.push_back(static_cast<char>('0' + n % 10));
        n /= 10;
    }
    if (negative) text.push_back('-');
    std::reverse(text.begin(), text.end());
    return text;
}

int factorial(int n) {
    int result = 1;
    for (int i=2; i<=n; ++i) result *= i;
    return result;
}

struct Complement {
    int mask, coefficient;
};

struct Tensor {
    std::array<std::vector<Entry>,14> edges;
    std::array<std::array<int,3>,7> symplectic;
    int components = 0;

    void read() {
        for (auto &p : symplectic) std::cin >> p[0] >> p[1] >> p[2];
        int count;
        std::cin >> count;
        for (int row=0; row<count; ++row) {
            Entry ix;
            int value;
            for (int &i : ix) std::cin >> i;
            std::cin >> value;
            assert(std::is_sorted(ix.begin(),ix.end()) && ix.front()>=0 && ix.back()<14);
            do {
                edges[ix[0]].push_back({1<<ix[1],1<<ix[2],1<<ix[3],value});
                ++components;
            } while (std::next_permutation(ix.begin(),ix.end()));
        }
        assert(std::cin && count == 44 && components == 840);
    }

    std::vector<Complement> complements(int vertices) const {
        std::vector<Complement> out;
        for (int selection=0; selection<128; ++selection) {
            if (__builtin_popcount(static_cast<unsigned>(selection)) != vertices/2) continue;
            int colour_mask=0, coefficient=1;
            std::vector<int> order;
            for (int p=0; p<7; ++p) {
                auto [i,j,b]=symplectic[p];
                if (selection & (1<<p)) {
                    colour_mask |= (1<<i)|(1<<j);
                } else {
                    coefficient *= b;
                    order.push_back(i);
                    order.push_back(j);
                }
            }
            for (int i=0; i<14; ++i) if (colour_mask & (1<<i)) order.push_back(i);
            for (int i=0; i<14; ++i) for (int j=0; j<i; ++j)
                if (order[j]>order[i]) coefficient = -coefficient;
            out.push_back({colour_mask,coefficient});
        }
        return out;
    }

    // Fix the order in the first flavour. In the other three flavours the
    // ordered determinant signs arise by moving each selected bit to the front.
    int64_t four_column_minor(int a, int b, int c, int d) const {
        if (!a) return 1;
        int i=__builtin_ctz(static_cast<unsigned>(a));
        int tail=a&(a-1);
        int64_t result=0;
        for (auto [j,k,l,value] : edges[i]) {
            if ((b&j) && (c&k) && (d&l)) {
                int sign=parity(((j-1)&b)^((k-1)&c)^((l-1)&d)) ? -1 : 1;
                result += sign*value*four_column_minor(tail,b^j,c^k,d^l);
            }
        }
        return result;
    }

    void run(int vertices) const {
        auto choices=complements(vertices);
        Wide numerator=0;
        int64_t tested=0, nonzero=0;
        for (int i=0; i<static_cast<int>(choices.size()); ++i) {
            for (int j=i; j<static_cast<int>(choices.size()); ++j)
            for (int k=j; k<static_cast<int>(choices.size()); ++k)
            for (int l=k; l<static_cast<int>(choices.size()); ++l) {
                Entry label={i,j,k,l};
                int copies=24;
                for (int index : std::set<int>(label.begin(),label.end()))
                    copies/=factorial(static_cast<int>(std::count(label.begin(),label.end(),index)));
                int64_t minor=four_column_minor(choices[i].mask,choices[j].mask,choices[k].mask,choices[l].mask);
                tested += copies;
                if (!minor) continue;
                nonzero += copies;
                Wide summand=static_cast<Wide>(copies)*minor;
                for (int v : label) summand*=choices[v].coefficient;
                numerator+=summand;
            }
        }
        const int power=7-vertices/2;
        // Vertices are even, so their exponential sign is positive. Ordering
        // by the first flavour already includes the divided power 1 / vertices!.
        for (int i=0; i<4; ++i) numerator*=factorial(power);
        const int denom=1<<vertices; // input quartic is 2 times the polarized tensor
        assert(numerator % denom == 0);
        std::cout << vertices << " " << power << " " << tested << " " << nonzero
                  << " " << decimal(numerator/denom) << "\n";
    }
};

int main() {
    Tensor t;
    t.read();
    for (int m : {2,4,6}) t.run(m);
}
