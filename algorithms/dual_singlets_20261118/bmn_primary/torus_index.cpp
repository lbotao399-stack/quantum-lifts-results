// Refined BMN Euler and primary indices by exact finite torus integration.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using Residue = std::uint32_t;
using Wide = std::uint64_t;
using Triple = std::array<int, 3>;
using Weight = std::array<int, 4>;

struct Field {
    Residue p;
    Residue multiply(Residue x, Residue y) const { return Wide(x) * y % p; }
    Residue add(Residue x, Residue y) const { return x + y >= p ? x + y - p : x + y; }
    Residue subtract(Residue x, Residue y) const { return x >= y ? x - y : x + p - y; }
    Residue power(Residue x, unsigned n) const {
        Residue value = 1;
        while (n) {
            if (n & 1) value = multiply(value, x);
            x = multiply(x, x);
            n >>= 1;
        }
        return value;
    }
};

struct RootSystem {
    std::vector<Weight> positive;
    int rank, weyl_order;
    bool signed_permutations;
};

RootSystem roots(const std::string &name) {
    if (name == "G2") return {{{1,0},{0,1},{1,1},{1,-1},{2,1},{1,2}}, 2,12,false};
    if (name == "A2") return {{{1,-1},{2,1},{1,2}}, 2,6,false};
    const std::map<std::string,int> orders{
        {"C2",8},{"B3",48},{"C3",48},{"B4",384},{"C4",384},{"D4",192},{"F4",1152}};
    if (!orders.count(name)) throw std::runtime_error("unsupported root system");
    RootSystem g{{},name[1]-'0',orders.at(name),true};
    const int scale = name=="F4" ? 2 : 1;
    for (int i=0; i<g.rank; ++i) {
        if (name[0]!='D') {
            Weight root{};
            root[i] = name[0]=='C' ? 2 : scale;
            g.positive.push_back(root);
        }
        for (int j=i+1; j<g.rank; ++j) for (int sign : {-1,1}) {
            Weight root{};
            root[i]=scale;
            root[j]=sign*scale;
            g.positive.push_back(root);
        }
    }
    if (name=="F4")
        for (int a : {-1,1}) for (int b : {-1,1}) for (int c : {-1,1})
            g.positive.push_back({1,a,b,c});
    return g;
}

std::map<std::vector<int>, Wide> torus(const RootSystem &group, int grid) {
    std::map<std::vector<int>, Wide> counts;
    auto contribute = [&](const Weight &point, int count) {
        std::vector<int> powers;
        for (const auto &r : group.positive) {
            int k=0;
            for (int i=0; i<group.rank; ++i) k+=r[i]*point[i];
            k = ((k % grid) + grid) % grid;
            if (k == 0) return; // The Weyl denominator vanishes.
            powers.push_back(k);
            powers.push_back(grid-k);
        }
        std::sort(powers.begin(), powers.end());
        counts[powers]+=count;
    };
    if (!group.signed_permutations) {
        for (int i=0;i<grid;++i) for (int j=0;j<grid;++j) contribute({i,j,0,0},1);
    } else {
        int factorial=1;
        for (int j=2;j<=group.rank;++j) factorial*=j;
        Weight point{};
        auto visit = [&](auto &&self, int axis, int first, int orbit) -> void {
            if (axis==group.rank) { contribute(point,orbit); return; }
            // Equal absolute coordinates lie on a root hyperplane in every listed group.
            for (int k=first;k<=grid/2;++k) {
                point[axis]=k;
                self(self,axis+1,k+1,orbit*(k==0 || k==grid/2 ? 1 : 2));
            }
        };
        visit(visit,0,0,factorial);
    }
    return counts;
}

struct Flavour {
    int side, degree;
    std::vector<Triple> weights;
    std::vector<int> indices;

    int index(int a, int b, int c) const { return (a*side+b)*side+c; }
    explicit Flavour(int maximum, int total) : side(maximum+3), degree(total) {
        for (int a = 0; a < side; ++a)
            for (int b = 0; b < side && a+b <= degree; ++b)
                for (int c = 0; c < side && a+b+c <= degree; ++c) {
                    weights.push_back({a,b,c});
                    indices.push_back(index(a,b,c));
                }
    }

    void convolve(std::vector<Residue> &values, std::vector<Residue> &out, Triple step,
                  const std::vector<Residue> &coeff, const Field &f) const {
        const int stride = index(step[0], step[1], step[2]);
        for (size_t i = 0; i < weights.size(); ++i) {
            int limit = int(coeff.size())-1;
            for (int j = 0; j < 3; ++j) if (step[j]) limit = std::min(limit, weights[i][j]);
            const int k = indices[i];
            // The supported primes are below 10^9. Reduce after each eight products.
            Wide sum = 0;
            for (int n = 0; n <= limit; ++n) {
                sum += Wide(coeff[n])*values[k-n*stride];
                if ((n & 7) == 7) sum %= f.p;
            }
            out[k] = sum % f.p;
        }
        values.swap(out);
    }
};

int main(int argc, char **argv) {
    if (argc != 6 && argc != 7) throw std::runtime_error("usage: torus_index group degree max_first_part prime grid [ordinary]");
    const bool ordinary = argc==7;
    if (ordinary && std::string(argv[6])!="ordinary") throw std::runtime_error("unknown index mode");
    const auto group = roots(argv[1]);
    const int degree = std::stoi(argv[2]), maximum = std::stoi(argv[3]), grid = std::stoi(argv[5]);
    const Field f{static_cast<Residue>(std::stoul(argv[4]))};
    int width = 0;
    for (int i = 0; i < group.rank; ++i) {
        int sum = 0, largest = 0;
        for (auto r : group.positive) {
            sum += std::abs(r[i]);
            largest = std::max(largest,std::abs(r[i]));
        }
        width = std::max(width, degree*largest+sum);
    }
    if (degree < 1 || maximum < 2 || maximum > degree || grid <= width ||
        (grid & (grid-1)) || f.p >= 1000000000 || (f.p-1)%grid)
        throw std::runtime_error("invalid integration grid, prime or degree");
    Residue root = 1;
    for (Residue x = 2; x < f.p && f.power(root,grid/2) == 1; ++x)
        root = f.power(x,(f.p-1)/grid);
    if (f.power(root,grid) != 1 || f.power(root,grid/2) == 1)
        throw std::runtime_error("no root of the requested order");
    std::vector<Residue> unity(grid,1);
    for (int k=1; k<grid; ++k) unity[k] = f.multiply(unity[k-1],root);
    Flavour flavour(maximum, degree);
    const int side = flavour.side, size = side*side*side;
    std::vector<Residue> integral(size), values(size), work(size);
    const auto points = torus(group,grid);
    for (const auto &[exponents, count] : points) {
        Residue weyl = count % f.p;
        std::vector<Residue> eigenvalues(group.rank,1);
        for (int k : exponents) {
            weyl = f.multiply(weyl, f.subtract(1,unity[k]));
            eigenvalues.push_back(unity[k]);
        }
        std::vector<Residue> exterior{1};
        for (Residue v : eigenvalues) {
            exterior.push_back(0);
            for (int j = int(exterior.size())-1; j>0; --j)
                exterior[j] = f.subtract(exterior[j], f.multiply(v,exterior[j-1]));
        }
        std::vector<Residue> symmetric(side), axis(side);
        symmetric[0] = axis[0] = 1;
        for (int n=1; n<side; ++n) {
            for (int k=1; k<=n && k<int(exterior.size()); ++k)
                symmetric[n] = f.subtract(symmetric[n], f.multiply(exterior[k],symmetric[n-k]));
            axis[n] = ordinary ? symmetric[n] : f.add(axis[n-1],symmetric[n]);
        }
        for (size_t i=0; i<flavour.weights.size(); ++i) {
            auto [a,b,c] = flavour.weights[i];
            values[flavour.indices[i]] = f.multiply(axis[a],f.multiply(axis[b],axis[c]));
        }
        for (Triple step : {Triple{1,1,0},Triple{0,1,1},Triple{1,0,1}})
            flavour.convolve(values,work,step,exterior,f);
        flavour.convolve(values,work,{1,1,1},symmetric,f);
        for (int k : flavour.indices) integral[k] = f.add(integral[k],f.multiply(weyl,values[k]));
    }
    Residue scale = f.power(f.multiply(group.weyl_order,f.power(grid,group.rank)),f.p-2);
    const int minimum = ordinary ? 0 : 2;
    for (int a=minimum; a<=maximum; ++a) for (int b=minimum; b<=a; ++b)
        for (int c=minimum; c<=b && a+b+c<=degree; ++c) {
            Residue total=0;
            Triple perm{0,1,2};
            do {
                const int i=a+2-perm[0], j=b+1-perm[1], k=c-perm[2];
                if (i<0 || j<0 || k<0) continue;
                const int sign = (perm[0]<perm[1])+(perm[0]<perm[2])+(perm[1]<perm[2]);
                const Residue v = integral[flavour.index(i,j,k)];
                total = sign%2 ? f.subtract(total,v) : f.add(total,v);
            } while (std::next_permutation(perm.begin(),perm.end()));
            if ((a+b+c)%2) total = f.subtract(0,total);
            std::cout << a << ' ' << b << ' ' << c << ' ' << f.multiply(scale,total) << '\n';
        }
}
