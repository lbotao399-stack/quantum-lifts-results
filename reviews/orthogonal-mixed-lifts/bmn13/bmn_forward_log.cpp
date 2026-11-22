// Expand each raw scalar trace of the half trace logarithm, then reverse it
// by adjacent swaps. Output the active primitive numerator over 4*length.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using Matrix = std::array<int64_t, 4>;
using Couplings = std::array<Matrix, 3>;
using Values = std::array<int64_t, 2>;

Matrix multiply(const Matrix& a, const Matrix& b) {
    return {a[0]*b[0]+a[1]*b[2], a[0]*b[1]+a[1]*b[3],
            a[2]*b[0]+a[3]*b[2], a[2]*b[1]+a[3]*b[3]};
}

std::set<std::string> active_traces(const char* path) {
    std::set<std::string> traces;
    std::ifstream input(path);
    assert(input);
    std::string line, word, trace;
    while (std::getline(input, line)) {
        std::istringstream row(line);
        row >> word;
        std::istringstream factors(word);
        while (std::getline(factors, trace, '/')) {
            if (trace.find_first_of("567") != std::string::npos) traces.insert(trace);
        }
    }
    return traces;
}

void expand(const std::string& raw, const Values& traces,
            const std::unordered_map<std::string, int>& requested, std::vector<Values>& result) {
    auto state = raw;
    for (int last = int(raw.size())-1; last > 0; --last) {
        for (int j = 0; j < last; ++j) {
            char a = state[j], b = state[j+1];
            if (a == b) continue;
            // 01, 02, 12 encode the three oriented matter second jets.
            char pair = '4' + (a-'0') + (b-'0');
            auto key = state.substr(0,j) + pair + state.substr(j+2);
            auto it = requested.find(key);
            if (it != requested.end()) {
                for (int seed = 0; seed < 2; ++seed)
                    result[it->second][seed] += (a<b ? traces[seed] : -traces[seed]);
            }
            std::swap(state[j], state[j+1]);
        }
    }
    assert(std::equal(raw.begin(), raw.end(), state.rbegin()));
}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "usage: forward_log parents.dat parameters.dat active.dat\n";
        return 2;
    }
    std::array<Couplings, 2> couplings;
    std::ifstream parameters(argv[2]);
    std::string line;
    for (int seed = 0; seed < 2; ++seed) {
        std::getline(parameters, line);
        std::istringstream row(line);
        int id;
        row >> id;
        assert(id == seed);
        for (auto& m : couplings[seed]) {
            row >> m[0] >> m[1] >> m[2];
            m[3] = -m[0];
        }
        assert(row);
    }
    auto active = active_traces(argv[1]);
    std::map<int, std::vector<std::string>> by_size;
    for (const auto& t : active) by_size[t.size()+1].push_back(t);
    std::ofstream out(argv[3]);
    assert(out);
    for (const auto& [size, words] : by_size) {
        assert(size % 2 && size <= 11);
        std::unordered_map<std::string, int> requested;
        for (int i = 0; i < int(words.size()); ++i) {
            for (auto s : {words[i], std::string(words[i].rbegin(), words[i].rend())}) {
                for (int k = 0; k < size-1; ++k) {
                    auto [it, fresh] = requested.emplace(s.substr(k)+s.substr(0,k), i);
                    assert(fresh || it->second == i);
                }
            }
        }
        int limit = 1;
        for (int i = 0; i < size; ++i) limit *= 3;
        std::vector<Values> result(words.size(), {0, 0});
        for (int id = 0; id < limit; ++id) {
            std::string raw(size, '0');
            int value = id;
            for (char& label : raw) {
                label += value % 3;
                value /= 3;
            }
            Values trace;
            for (int seed = 0; seed < 2; ++seed) {
                Matrix matrix{1,0,0,1};
                for (char c : raw) matrix = multiply(matrix, couplings[seed][c-'0']);
                trace[seed] = matrix[0]+matrix[3];
            }
            if (trace[0] || trace[1]) expand(raw, trace, requested, result);
        }
        for (int i = 0; i < int(words.size()); ++i)
            out << words[i] << ' ' << 4*size << ' ' << result[i][0] << ' ' << result[i][1] << '\n';
        std::cout << "scalar length " << size << ", raw words " << limit
                  << ", active traces " << words.size() << '\n';
    }
}
