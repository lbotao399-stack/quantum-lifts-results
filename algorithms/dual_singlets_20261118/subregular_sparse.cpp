// Construct the Taylor differential and solve its coefficient equations modulo a prime.
#include <algorithm>
#include <array>
#include <cstdint>
#include <deque>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using Word = std::string;
using Row = std::vector<std::pair<int, int>>;
struct Rule {
    std::array<int, 3> factors{};
    int length = 2;
    int coefficient = 0;
};
using Rules = std::array<std::vector<Rule>, 128>;
constexpr char digit[] = "0123456789abcdef";

Word decode(const std::string& hex) {
    Word w;
    for (size_t i = 0; i < hex.size(); i += 2)
        w.push_back(static_cast<char>(std::stoi(hex.substr(i, 2), nullptr, 16)));
    return w;
}
std::string encode(const Word& w) {
    std::string hex;
    for (unsigned char c : w) { hex += digit[c >> 4]; hex += digit[c & 15]; }
    return hex;
}
int mod(std::int64_t n, int prime) {
    n %= prime;
    return static_cast<int>(n < 0 ? n + prime : n);
}
int inverse(int a, int prime) {
    int result = 1;
    for (int power = prime-2; power; power >>= 1) {
        if (power & 1) result = mod(std::int64_t(result)*a, prime);
        a = mod(std::int64_t(a)*a, prime);
    }
    return result;
}

std::unordered_map<Word, int> variation(const Word& parent, const Rules& rules, const std::array<int,128>& odd) {
    std::unordered_map<Word,int> image;
    int prefix = 1;
    for (size_t pos = 0; pos < parent.size(); ++pos) {
        int old = parent[pos];
        for (auto rule : rules[old]) {
            int sign = prefix;
            for (size_t h = 0; h < parent.size() && sign; ++h) {
                int x = parent[h];
                if (h == pos || !odd[x]) continue;
                for (int b=0; b<rule.length; ++b) {
                    int k = rule.factors[b];
                    if (!odd[k]) continue;
                    if (k == x) { sign = 0; break; }
                    if ((h < pos && x > k) || (h > pos && x < k)) sign = -sign;
                }
            }
            if (!sign) continue;
            Word target = parent;
            target.erase(pos, 1);
            for (int b=0; b<rule.length; ++b) {
                int k = rule.factors[b];
                target.insert(std::lower_bound(target.begin(), target.end(), k), static_cast<char>(k));
            }
            image[target] += sign*rule.coefficient;
        }
        if (odd[old]) prefix = -prefix;
    }
    return image;
}

// Each bit in a higher level represents one nonempty 64-bit word below it.
// During elimination the next pivot is the least active column. A zero
// coefficient may keep its bit until extraction, which avoids scanning tails.
class ActiveColumns {
    std::vector<std::vector<std::uint64_t>> levels;
public:
    explicit ActiveColumns(size_t columns) {
        size_t words=std::max<size_t>(1,(columns+63)/64);
        levels.emplace_back(words,0);
        while (words>1) {
            words=(words+63)/64;
            levels.emplace_back(words,0);
        }
    }
    bool empty() const { return !levels.back()[0]; }
    void add(size_t column) {
        for (auto& level : levels) {
            const auto bit=std::uint64_t{1}<<(column%64);
            auto& word=level[column/64];
            const bool was_empty=!word;
            word|=bit;
            if (!was_empty) break;
            column/=64;
        }
    }
    int pop_least() {
        size_t position=0;
        for (auto it=levels.rbegin(); it!=levels.rend(); ++it)
            position=64*position+__builtin_ctzll((*it)[position]);
        const int column=static_cast<int>(position);
        for (auto& level : levels) {
            auto& word=level[position/64];
            word&=~(std::uint64_t{1}<<(position%64));
            if (word) break;
            position/=64;
        }
        return column;
    }
};

// A basis of the alternating quotient for exchanging the two even Taylor directions.
// At a fixed word the colour-coordinate exterior sign decides whether it survives.
int alternating_representative(Word& w, const std::array<int,128>& odd) {
    Word swapped = w;
    int sign = 1;
    for (size_t i=0; i<w.size(); ++i) {
        const int v = w[i];
        const char x = static_cast<char>((v & ~12) | ((v & 4)<<1) | ((v & 8)>>1));
        size_t pos=i;
        while (pos && swapped[pos-1] > x) {
            if (odd[x] && odd[swapped[pos-1]]) sign=-sign;
            swapped[pos]=swapped[pos-1];
            --pos;
        }
        swapped[pos]=x;
    }
    if (swapped == w) return sign == -1 ? 1 : 0;
    if (swapped > w) return 1;
    w = std::move(swapped);
    return -sign;
}

int main(int argc, char** argv) {
    const bool fixed_support = argc >= 7 && std::string(argv[argc-2]) == "--support";
    std::vector<Word> supplied;
    if (fixed_support) {
        std::ifstream previous(argv[argc-1]);
        int rows, columns, singles, rank;
        if (!(previous >> rows >> columns >> singles >> rank) || rows < 1 || columns < 1)
            throw std::runtime_error("invalid support header");
        std::string hex;
        int residue;
        while (previous >> hex >> residue) {
            if (hex.empty() || hex.size()%2 || !residue)
                throw std::runtime_error("invalid support entry");
            const Word w = decode(hex);
            if (std::any_of(w.begin(), w.end(), [](unsigned char k) { return k >= 128; })
                    || !std::is_sorted(w.begin(),w.end()))
                throw std::runtime_error("unexpected Taylor generator in support");
            supplied.push_back(w);
        }
        std::sort(supplied.begin(),supplied.end());
        if (!previous.eof() || supplied.empty() || supplied.size() > static_cast<size_t>(columns)
                || std::adjacent_find(supplied.begin(),supplied.end()) != supplied.end())
            throw std::runtime_error("invalid previous solution");
        argc -= 2;
    }
    const bool spin_odd = argc>5 && std::string(argv[argc-1])=="spin-odd";
    if (spin_odd) --argc;
    const bool mixed = argc>=6 && argc<=8 && std::string(argv[5])=="mixed";
    if (argc!=5 && !mixed)
        throw std::runtime_error("expected INPUT OUTPUT PRIME INCIDENCE_ROUNDS [mixed [MIN_ROOTS [MAX_ROOTS]]] [spin-odd] [--support PREVIOUS]");
    const int min_roots = argc>=7 ? std::stoi(argv[6]) : 0;
    const int max_roots = argc==8 ? std::stoi(argv[7]) : 4;
    if (min_roots<0 || min_roots>3 || max_roots<min_roots || max_roots>4)
        throw std::runtime_error("invalid range of off-diagonal ghost degrees");
    // f and g occupy the last two groups of sixteen Taylor coefficients.
    auto roots = [](const Word& w) { return w.end()-std::lower_bound(w.begin(),w.end(),6*16); };
    int prime=std::stoi(argv[3]), rounds=std::stoi(argv[4]);
    if (fixed_support && rounds != 1)
        throw std::runtime_error("a fixed support is supplied for one round");
    std::ifstream in(argv[1]);
    std::array<int,128> odd{};
    for (int& v : odd) in >> v;
    Rules rules;
    std::array<std::vector<int>,128*128> backwards;
    std::unordered_map<int,std::vector<int>> cubic;
    int number;
    in >> number;
    for (int h=0; h<number; ++h) {
        int k;
        Rule rule;
        in >> k;
        if (mixed) in >> rule.length;
        if (rule.length!=2 && rule.length!=3) throw std::runtime_error("invalid Taylor rule length");
        int encoded=0;
        for (int b=0; b<rule.length; ++b) {
            int& x=rule.factors[b];
            in >> x;
            if (x<0 || x>=128 || (b && x<rule.factors[b-1])) throw std::runtime_error("invalid Taylor factors");
            encoded=128*encoded+x;
        }
        in >> rule.coefficient;
        rules.at(k).push_back(rule);
        if (rule.length==2) backwards.at(encoded).push_back(k);
        else cubic[encoded].push_back(k);
    }
    int rhs_count;
    in >> rhs_count;
    std::vector<int> rhs;
    std::unordered_map<Word,int> row_id;
    std::vector<Word> targets;
    auto append = [&](const Word& w) {
        auto it = row_id.emplace(w, static_cast<int>(targets.size()));
        if (it.second) { targets.push_back(w); rhs.push_back(0); }
        return it.first->second;
    };
    for (int t=0; t<rhs_count; ++t) {
        std::string hex; std::int64_t n; in >> hex >> n;
        Word w=decode(hex);
        if (argc==8 && roots(w)>max_roots) throw std::runtime_error("right hand side outside quotient");
        const int sign = spin_odd ? alternating_representative(w,odd) : 1;
        if (sign) {
            int i=append(w);
            rhs[i]=mod(rhs[i]+sign*(n%prime),prime);
        }
    }
    if (!in) throw std::runtime_error("incomplete input");
    std::vector<Word> parents;
    std::unordered_set<Word> seen;
    std::vector<Row> columns;
    size_t first=0;
    for (int layer=0; layer<rounds && first<targets.size(); ++layer) {
        size_t stop=targets.size();
        std::vector<Word> fresh = std::move(supplied);
        for (size_t t=first; !fixed_support && t<stop; ++t) {
            const Word& w=targets[t];
            auto contract = [&](const std::vector<int>& keys, size_t i, size_t j, size_t kpos) {
                if (keys.empty()) return;
                Word rest=w;
                if (kpos<w.size()) rest.erase(kpos,1);
                rest.erase(j,1); rest.erase(i,1);
                for (int k : keys) {
                    if (odd[k] && std::binary_search(rest.begin(),rest.end(),k)) continue;
                    Word p=rest; p.insert(std::lower_bound(p.begin(),p.end(),k),static_cast<char>(k));
                    if (min_roots && roots(p)<min_roots) continue;
                    if (spin_odd && !alternating_representative(p,odd)) continue;
                    if (seen.insert(p).second) fresh.push_back(std::move(p));
                }
            };
            for (size_t i=0; i<w.size(); ++i) for (size_t j=i+1; j<w.size(); ++j) {
                const int pair=128*int(w[i])+int(w[j]);
                contract(backwards[pair],i,j,w.size());
                if (!cubic.empty()) for (size_t k=j+1; k<w.size(); ++k) {
                    auto it=cubic.find(128*pair+int(w[k]));
                    if (it!=cubic.end()) contract(it->second,i,j,k);
                }
            }
        }
        std::sort(fresh.begin(), fresh.end());
        for (const Word& p : fresh) {
            auto image=variation(p,rules,odd);
            if (spin_odd) {
                std::unordered_map<Word,int> projected;
                for (const auto& term : image) {
                    Word w=term.first;
                    int sign=alternating_representative(w,odd);
                    if (sign) projected[w]+=sign*term.second;
                }
                image=std::move(projected);
            }
            std::vector<std::pair<Word,int>> ordered(image.begin(),image.end());
            std::sort(ordered.begin(),ordered.end());
            Row column;
            for (const auto& item : ordered)
                if (item.second && (argc!=8 || roots(item.first)<=max_roots))
                    column.emplace_back(append(item.first),item.second);
            if (!column.empty()) { parents.push_back(p); columns.push_back(std::move(column)); }
        }
        first=stop;
    }
    const int m=static_cast<int>(rhs.size()), n=static_cast<int>(parents.size());
    std::vector<Row> rows(m);
    std::vector<int> degree(n);
    std::deque<int> single;
    for (int j=0; j<n; ++j) {
        degree[j]=static_cast<int>(columns[j].size());
        if (degree[j]==1) single.push_back(j);
        for (auto [i,c] : columns[j]) rows[i].emplace_back(j,c);
    }
    std::vector<bool> alive(m,true);
    std::vector<std::pair<int,int>> peeled;
    while (!single.empty()) {
        int j=single.front(); single.pop_front();
        if (degree[j]!=1) continue;
        for (auto [i,c] : columns[j]) if (alive[i]) {
            alive[i]=false; peeled.emplace_back(i,j);
            for (auto [other,value] : rows[i]) if (--degree[other]==1) single.push_back(other);
            break;
        }
    }
    std::vector<int> order(n), indices(n), core;
    std::iota(order.begin(),order.end(),0);
    std::stable_sort(order.begin(),order.end(),[&](int x,int y){return degree[x]<degree[y];});
    for (int i=0; i<n; ++i) indices[order[i]]=i;
    for (int i=0; i<m; ++i) if (alive[i]) core.push_back(i);
    std::stable_sort(core.begin(),core.end(),[&](int x,int y){return rows[x].size()<rows[y].size();});
    std::cerr << m << " equations " << n << " parents " << peeled.size() << " singleton pivots\n";
    std::vector<Row> pivots(n);
    std::vector<int> values(n), solution(n), accumulator(n);
    ActiveColumns active(n);
    int rank=0;
    for (int i : core) {
        for (auto [column,c] : rows[i]) {
            int j=indices[column];
            accumulator[j]=mod(c,prime);
            active.add(j);
        }
        int y=rhs[i];
        int new_pivot=-1;
        while (!active.empty()) {
            int j=active.pop_least(), factor=accumulator[j];
            accumulator[j]=0;
            if (!factor) continue;
            if (pivots[j].empty()) {
                const int inv=inverse(factor,prime);
                Row& normalized=pivots[j];
                normalized.emplace_back(j,1);
                while (!active.empty()) {
                    int k=active.pop_least(), c=accumulator[k];
                    accumulator[k]=0;
                    if (c) normalized.emplace_back(k,mod(std::int64_t(c)*inv,prime));
                }
                values[j]=mod(std::int64_t(y)*inv,prime);
                new_pivot=j;
                ++rank;
                break;
            }
            const Row& pivot=pivots[j];
            for (size_t pos=1; pos<pivot.size(); ++pos) {
                auto [k,c]=pivot[pos];
                if (!accumulator[k]) active.add(k);
                accumulator[k]=mod(accumulator[k]-std::int64_t(factor)*c,prime);
            }
            y=mod(y-std::int64_t(factor)*values[j],prime);
        }
        if (new_pivot<0 && y) { std::cerr << "inconsistent in this parent span\n"; return 2; }
    }
    for (int j=n-1; j>=0; --j) if (!pivots[j].empty()) {
        int y=values[j];
        for (auto [k,c] : pivots[j]) if (k!=j)
            y=mod(y-std::int64_t(c)*solution[order[k]],prime);
        solution[order[j]]=y;
    }
    for (auto it=peeled.rbegin(); it!=peeled.rend(); ++it) {
        auto [i,j]=*it;
        int y=rhs[i], pivot=0;
        for (auto [k,c] : rows[i]) {
            if (k==j) pivot=c;
            else y=mod(y-std::int64_t(c)*solution[k],prime);
        }
        solution[j]=mod(std::int64_t(y)*inverse(mod(pivot,prime),prime),prime);
    }
    for (int i=0; i<m; ++i) {
        int y=rhs[i];
        for (auto [j,c] : rows[i]) y=mod(y-std::int64_t(c)*solution[j],prime);
        if (y) throw std::runtime_error("modular solution failed replay");
    }
    std::ofstream out(argv[2]);
    out << m << " " << n << " " << peeled.size() << " " << rank << "\n";
    for (int j=0; j<n; ++j) if (solution[j]) out << encode(parents[j]) << " " << solution[j] << "\n";
    if (!out) throw std::runtime_error("could not save solution");
}
