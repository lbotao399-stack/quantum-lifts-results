// Exact polynomial and exterior arithmetic for three-colour trace products.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

struct Monomial {
    std::array<unsigned char, 16> bytes{};
    bool operator==(const Monomial& other) const { return bytes == other.bytes; }
    int size() const { return bytes[15]; }
};

struct Hash {
    std::size_t operator()(const Monomial& m) const {
        std::uint64_t first, second;
        std::memcpy(&first, m.bytes.data(), 8);
        std::memcpy(&second, m.bytes.data()+8, 8);
        first ^= second + 0x9e3779b97f4a7c15ULL + (first << 6) + (first >> 2);
        return first ^ (first >> 29);
    }
};

using Integer = std::int64_t;
using Polynomial = std::unordered_map<Monomial,Integer,Hash>;
using Matrix = std::array<Polynomial,9>;
static std::array<int,128> odd;
static std::unordered_map<char,Matrix> letters;
static std::unordered_map<std::string,Matrix> matrices;
static std::unordered_map<std::string,Polynomial> traces;
static std::size_t cached_entries = 0;
static int degree_limit, root_limit;

static Integer times(Integer x, Integer y) {
    Integer value;
    if (__builtin_mul_overflow(x,y,&value)) throw std::overflow_error("coefficient product");
    return value;
}

static void add(Polynomial& out, const Monomial& m, Integer value) {
    if (!value) return;
    auto insertion = out.try_emplace(m,value);
    if (insertion.second) return;
    Integer& current = insertion.first->second;
    if (__builtin_add_overflow(current,value,&current)) throw std::overflow_error("coefficient sum");
    if (!current) out.erase(insertion.first);
}

static int multiply(const Monomial& left, const Monomial& right, Monomial& result) {
    const int n = left.size(), k = right.size(), length = n+k;
    if (length > degree_limit) return 0;
    if (root_limit < 4) {
        int roots = 0;
        for (int i=n-1; i>=0 && left.bytes[i]>=96; --i) ++roots;
        for (int j=k-1; j>=0 && right.bytes[j]>=96; --j) ++roots;
        if (roots > root_limit) return 0;
    }
    int i=0, j=0, right_parity=0, inversions=0;
    for (int s=0; s<length; ++s) {
        if (j==k || (i<n && left.bytes[i]<=right.bytes[j])) {
            const auto a=left.bytes[i++];
            if (j<k && a==right.bytes[j] && odd[a]) return 0;
            inversions ^= right_parity & odd[a];
            result.bytes[s]=a;
        } else {
            const auto a=right.bytes[j++];
            right_parity ^= odd[a];
            result.bytes[s]=a;
        }
    }
    result.bytes[15]=length;
    return inversions ? -1 : 1;
}

static void product(Polynomial& out, const Polynomial& left, const Polynomial& right) {
    for (const auto& [m,x]:left) for (const auto& [n,y]:right) {
        Monomial combined;
        const int sign=multiply(m,n,combined);
        if (sign) add(out,combined,times(sign,times(x,y)));
    }
}

static const Matrix& word_matrix(const std::string& word) {
    if (word.size()==1) return letters.at(word[0]);
    auto it=matrices.find(word);
    if (it!=matrices.end()) return it->second;
    const Matrix& head=word_matrix(word.substr(0,word.size()-1));
    const Matrix& tail=letters.at(word.back());
    Matrix result;
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) {
        for (int k=0;k<3;++k) product(result[3*i+j],head[3*i+k],tail[3*k+j]);
        cached_entries+=result[3*i+j].size();
    }
    return matrices.emplace(word,std::move(result)).first->second;
}

static Polynomial trace(const std::string& w) {
    const auto mid=w.size()/2;
    const Matrix& head=word_matrix(w.substr(0,mid));
    const Matrix& tail=word_matrix(w.substr(mid));
    Polynomial result;
    for (int i=0;i<3;++i) for (int j=0;j<3;++j)
        product(result,head[3*i+j],tail[3*j+i]);
    return result;
}

static const Polynomial& short_trace(const std::string& word) {
    auto it=traces.find(word);
    if (it!=traces.end()) return it->second;
    return traces.emplace(word,trace(word)).first->second;
}

static std::string decode(const std::string& text) {
    if (text.size()%2) throw std::invalid_argument("odd hexadecimal length");
    std::string result;
    for (std::size_t i=0;i<text.size();i+=2)
        result.push_back(static_cast<char>(std::stoul(text.substr(i,2),nullptr,16)));
    return result;
}

static void read_letters(std::istream& in) {
    for (int& k:odd) in>>k;
    int count;
    in>>count;
    for (int t=0;t<count;++t) {
        int id;
        in>>id;
        Matrix letter;
        for (auto& polynomial:letter) {
            int terms;
            in>>terms;
            for (int i=0;i<terms;++i) {
                int size;
                Integer c;
                in>>c>>size;
                if (size<0 || size>2) throw std::runtime_error("invalid letter degree");
                Monomial m;
                m.bytes[15]=size;
                for (int j=0;j<size;++j) {
                    int value;
                    in>>value;
                    if (value<0 || value>=128) throw std::runtime_error("bad generator");
                    m.bytes[j]=value;
                }
                add(polynomial,m,c);
            }
        }
        letters.emplace(static_cast<char>(id),std::move(letter));
    }
}

int main(int argc, char** argv) {
    if (argc!=7) return 1;
    std::ifstream in(argv[1]);
    in.exceptions(std::ios::badbit | std::ios::failbit);
    degree_limit=std::stoi(argv[3]);
    root_limit=std::stoi(argv[4]);
    int first=std::stoi(argv[5]), stop=std::stoi(argv[6]), group_count;
    if (degree_limit<1 || degree_limit>15 || root_limit<0 || root_limit>4) return 1;
    read_letters(in);
    in>>group_count;
    if (first<0 || stop>group_count || first>=stop) return 1;
    Polynomial total;
    for (int group=0;group<stop;++group) {
        std::string longest;
        int rows;
        in>>longest>>rows;
        Polynomial remaining;
        for (int row=0;row<rows;++row) {
            Integer coefficient;
            int n;
            in>>coefficient>>n;
            Polynomial spectator;
            if (group>=first) spectator[Monomial{}]=coefficient;
            for (int i=0;i<n;++i) {
                std::string word;
                in>>word;
                if (spectator.empty()) continue;
                Polynomial next;
                product(next,spectator,short_trace(decode(word)));
                spectator.swap(next);
            }
            for (const auto& [m,c]:spectator) add(remaining,m,c);
        }
        if (!remaining.empty()) product(total,trace(decode(longest)),remaining);
        if (cached_entries>500000) { matrices.clear(); cached_entries=0; }
        if (group>=first && (group-first+1)%10000==0)
            std::cerr<<group+1<<"/"<<stop<<" groups, "<<total.size()<<" monomials\n";
    }
    std::ofstream out(argv[2],std::ios::binary);
    out.exceptions(std::ios::badbit | std::ios::failbit);
    // Each monomial is 15 variable slots and one length byte, followed by a signed little-endian integer.
    for (const auto& [m,c]:total) {
        out.write(reinterpret_cast<const char*>(m.bytes.data()),16);
        std::uint64_t bits=static_cast<std::uint64_t>(c);
        for (int i=0;i<8;++i) out.put(static_cast<char>((bits>>(8*i))&255));
    }
}
