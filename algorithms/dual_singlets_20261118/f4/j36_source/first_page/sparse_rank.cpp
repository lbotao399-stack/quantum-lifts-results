// Exact sparse pivoting over F_1000003. Input is a little-endian uint32
// row count, column count, then a length and (column, residue) pairs per row.
#include <cstdint>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

constexpr int prime = 1000003;
using Residues = std::unordered_map<int,int>;
using Neighbours = std::unordered_set<int>;
using Entry = std::pair<int,int>;

static uint32_t read_word() {
    unsigned char b[4];
    if (!std::cin.read(reinterpret_cast<char*>(b), 4)) throw std::runtime_error("short input");
    return uint32_t(b[0]) | uint32_t(b[1])<<8 | uint32_t(b[2])<<16 | uint32_t(b[3])<<24;
}

static int multiply(int x, int y) {
    return int(int64_t(x)*y%prime);
}

static int inverse(int x) {
    int y = 1;
    for (int power = prime-2; power; power >>= 1) {
        if (power&1) y = multiply(y,x);
        x = multiply(x,x);
    }
    return y;
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    const int height = read_word(), width = read_word();
    std::vector<Residues> rows(height);
    std::vector<Neighbours> cols(width);
    for (int i=0; i<height; ++i) {
        const auto size = read_word();
        for (uint32_t k=0; k<size; ++k) {
            const int j = read_word(), v = read_word();
            if (j<0 || j>=width || v<=0 || v>=prime || !rows[i].emplace(j,v).second)
                throw std::runtime_error("invalid coordinate");
            cols[j].insert(i);
        }
    }
    if (std::cin.peek()!=std::char_traits<char>::eof()) throw std::runtime_error("extra input");
    // A stale heap entry is discarded when the current degree has changed.
    std::priority_queue<Entry,std::vector<Entry>,std::greater<Entry>> queue;
    auto degree = [&](int id) -> int {
        return int(id<height ? rows[id].size() : cols[id-height].size());
    };
    auto enqueue = [&](int id) {
        int n=degree(id);
        if (n) queue.emplace(n,id);
    };
    for (int id=0; id<height+width; ++id) enqueue(id);
    int rank=0;
    while (!queue.empty()) {
        const auto [n,id] = queue.top();
        queue.pop();
        if (degree(id)!=n) continue;
        int r,c;
        if (id<height) {
            r=id;
            c=rows[r].begin()->first;
            for (auto [j,v]:rows[r])
                if (Entry(cols[j].size(),j)<Entry(cols[c].size(),c)) c=j;
        } else {
            c=id-height;
            r=*cols[c].begin();
            for (int i:cols[c])
                if (Entry(rows[i].size(),i)<Entry(rows[r].size(),r)) r=i;
        }
        const int inv = inverse(rows[r].at(c));
        Residues pivot = std::move(rows[r]);
        Neighbours incident = std::move(cols[c]);
        rows[r].clear();
        cols[c].clear();
        pivot.erase(c);
        incident.erase(r);
        for (auto [j,v]:pivot) cols[j].erase(r);
        for (int i:incident) {
            auto &row=rows[i];
            int factor=multiply(row.at(c),inv);
            row.erase(c);
            for (auto [j,v]:pivot) {
                auto it=row.find(j);
                int old = it==row.end() ? 0 : it->second;
                int updated = old-multiply(factor,v);
                if (updated<0) updated+=prime;
                if (!updated) {
                    row.erase(it);
                    cols[j].erase(i);
                } else if (it==row.end()) {
                    row.emplace(j,updated);
                    cols[j].insert(i);
                } else {
                    it->second=updated;
                }
            }
            enqueue(i);
        }
        for (auto [j,v]:pivot) enqueue(height+j);
        ++rank;
    }
    std::cout << rank << '\n';
}
