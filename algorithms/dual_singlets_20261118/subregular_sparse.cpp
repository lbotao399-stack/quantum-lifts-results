// Construct the quadratic Taylor matrix and solve its coefficient equations modulo a prime.
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
struct Rule { int i, j, coefficient; };
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
                for (int k : {rule.i, rule.j}) {
                    if (!odd[k]) continue;
                    if (k == x) { sign = 0; break; }
                    if ((h < pos && x > k) || (h > pos && x < k)) sign = -sign;
                }
            }
            if (!sign) continue;
            Word target = parent;
            target.erase(pos, 1);
            target.insert(std::lower_bound(target.begin(), target.end(), rule.i), static_cast<char>(rule.i));
            target.insert(std::lower_bound(target.begin(), target.end(), rule.j), static_cast<char>(rule.j));
            image[target] += sign*rule.coefficient;
        }
        if (odd[old]) prefix = -prefix;
    }
    return image;
}

// Subtract factor times a normalized pivot. Both supports are sorted.
void eliminate(Row& row, const Row& pivot, int factor, int prime) {
    Row out;
    out.reserve(row.size() + pivot.size());
    size_t i=0, j=0;
    while (i<row.size() || j<pivot.size()) {
        if (j==pivot.size() || (i<row.size() && row[i].first<pivot[j].first)) {
            out.push_back(row[i++]);
        } else {
            int k = pivot[j].first;
            std::int64_t n = -std::int64_t(factor)*pivot[j++].second;
            if (i<row.size() && row[i].first==k) n += row[i++].second;
            int value=mod(n, prime);
            if (value) out.emplace_back(k, value);
        }
    }
    row.swap(out);
}

int main(int argc, char** argv) {
    if (argc != 5) throw std::runtime_error("expected INPUT OUTPUT PRIME INCIDENCE_ROUNDS");
    int prime=std::stoi(argv[3]), rounds=std::stoi(argv[4]);
    std::ifstream in(argv[1]);
    std::array<int,128> odd{};
    for (int& v : odd) in >> v;
    Rules rules;
    std::array<std::vector<int>,128*128> backwards;
    int number;
    in >> number;
    for (int h=0; h<number; ++h) {
        int k,i,j,c; in >> k >> i >> j >> c;
        rules.at(k).push_back({i,j,c});
        backwards.at(128*i+j).push_back(k);
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
        int i=append(decode(hex)); rhs[i]=mod(n,prime);
    }
    if (!in) throw std::runtime_error("incomplete input");
    std::vector<Word> parents;
    std::unordered_set<Word> seen;
    std::vector<Row> columns;
    size_t first=0;
    for (int layer=0; layer<rounds && first<targets.size(); ++layer) {
        size_t stop=targets.size();
        std::vector<Word> fresh;
        for (size_t t=first; t<stop; ++t) {
            const Word& w=targets[t];
            for (size_t i=0; i<w.size(); ++i) for (size_t j=i+1; j<w.size(); ++j) {
                const auto& keys=backwards[128*int(w[i])+int(w[j])];
                if (keys.empty()) continue;
                Word rest=w; rest.erase(j,1); rest.erase(i,1);
                for (int k : keys) {
                    if (odd[k] && std::binary_search(rest.begin(),rest.end(),k)) continue;
                    Word p=rest; p.insert(std::lower_bound(p.begin(),p.end(),k),static_cast<char>(k));
                    if (seen.insert(p).second) fresh.push_back(std::move(p));
                }
            }
        }
        std::sort(fresh.begin(), fresh.end());
        for (const Word& p : fresh) {
            auto image=variation(p,rules,odd);
            std::vector<std::pair<Word,int>> ordered(image.begin(),image.end());
            std::sort(ordered.begin(),ordered.end());
            Row column;
            for (const auto& item : ordered)
                if (item.second) column.emplace_back(append(item.first),item.second);
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
    std::vector<int> values(n), solution(n);
    int rank=0;
    for (int i : core) {
        Row row=rows[i];
        for (auto& x : row) { x.first=indices[x.first]; x.second=mod(x.second,prime); }
        std::sort(row.begin(),row.end());
        int y=rhs[i];
        while (!row.empty() && !pivots[row.front().first].empty()) {
            int j=row.front().first, factor=row.front().second;
            eliminate(row,pivots[j],factor,prime);
            y=mod(y-std::int64_t(factor)*values[j],prime);
        }
        if (row.empty()) {
            if (y) { std::cerr << "inconsistent in this parent span\n"; return 2; }
        } else {
            int j=row.front().first, inv=inverse(row.front().second,prime);
            for (auto& x : row) x.second=mod(std::int64_t(x.second)*inv,prime);
            pivots[j]=std::move(row); values[j]=mod(std::int64_t(y)*inv,prime);
            ++rank;
        }
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
