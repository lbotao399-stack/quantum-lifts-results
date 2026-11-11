// spanp_dual.cpp: extends spanp.cpp with combos == 2 (dual witness mode):
// Runs fast combos == 0 sparse elimination on base columns, then for each independent ("F") test column,
// outputs "F_PHI i_star r1 v1 r2 v2 ..." giving a sparse functional phi over F_p satisfying:
//   phi^T * base_col == 0 for every base column, and phi^T * test_col == 1 (mod p),
// constructed in O(nnz(cols)) time by deterministic reverse pivot back-substitution!
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <queue>
#include <algorithm>
#include <functional>
#include <string>

typedef uint32_t u32;
typedef uint64_t u64;
typedef std::vector<std::pair<u32, u32>> SV;

static u64 P;
static u64 pw(u64 a, u64 e) { u64 r = 1; a %= P; while (e) { if (e & 1) r = r * a % P; a = a * a % P; e >>= 1; } return r; }

int main() {
    std::vector<u32> buf;
    {
        const size_t CH = 1 << 20;
        std::vector<u32> tmp(CH);
        size_t got;
        while ((got = fread(tmp.data(), 4, CH, stdin)) > 0) buf.insert(buf.end(), tmp.begin(), tmp.begin() + got);
    }
    size_t at = 0;
    auto rd = [&]() -> u32 { if (at >= buf.size()) { fprintf(stderr, "spanp_dual: truncated input\n"); exit(2); } return buf[at++]; };
    P = rd();
    u32 nrows = rd(), nbase = rd(), ntest = rd(), mode = rd();
    bool track_combos = (mode == 1);
    bool dual_phi = (mode == 2);
    u32 ncols = nbase + ntest;
    std::vector<u64> acc(nrows, 0);
    std::vector<char> mark(nrows, 0);
    std::vector<int> piv(nrows, -1);
    std::vector<SV> cols;
    std::vector<SV> cco;
    std::vector<u64> cacc(track_combos ? ncols : 0, 0);
    std::vector<char> cmark(track_combos ? ncols : 0, 0);
    std::vector<u32> ctouch;
    std::string out;
    char line[64];
    u32 base_rank = 0;

    for (u32 j = 0; j < ncols; ++j) {
        u32 nnz = rd();
        std::priority_queue<u32, std::vector<u32>, std::greater<u32>> h;
        for (u32 t = 0; t < nnz; ++t) {
            u32 r = rd(), v = rd();
            if (r >= nrows) { fprintf(stderr, "spanp_dual: row out of range\n"); return 2; }
            acc[r] = (acc[r] + v) % P;
            if (!mark[r]) { mark[r] = 1; h.push(r); }
        }
        ctouch.clear();
        if (track_combos) { cacc[j] = 1; cmark[j] = 1; ctouch.push_back(j); }
        bool fresh = false;
        u32 fresh_row = 0;
        u64 fresh_val = 0;

        while (!h.empty()) {
            u32 i = h.top(); h.pop(); mark[i] = 0;
            u64 a = acc[i];
            if (!a) continue;
            int r = piv[i];
            if (r < 0) {
                fresh_row = i;
                fresh_val = a;
                if (j >= nbase && dual_phi) {
                    // Do NOT mutate base cols if we only want dual phi for this test column:
                    acc[i] = 0;
                    while (!h.empty()) {
                        u32 k = h.top(); h.pop(); mark[k] = 0;
                        acc[k] = 0;
                    }
                    fresh = true;
                    break;
                }
                SV v; v.emplace_back(i, (u32)a); acc[i] = 0;
                while (!h.empty()) {
                    u32 k = h.top(); h.pop(); mark[k] = 0;
                    if (acc[k]) { v.emplace_back(k, (u32)acc[k]); acc[k] = 0; }
                }
                u64 ai = pw(a, P - 2);
                for (auto &e : v) e.second = (u32)(e.second * ai % P);
                SV c;
                if (track_combos) {
                    for (u32 k : ctouch) { if (cacc[k]) c.emplace_back(k, (u32)(cacc[k] * ai % P)); cacc[k] = 0; cmark[k] = 0; }
                }
                piv[i] = (int)cols.size();
                cols.push_back(std::move(v));
                if (track_combos) cco.push_back(std::move(c));
                fresh = true;
                break;
            }
            u64 na = P - a;
            acc[i] = 0;
            const SV &pc = cols[r];
            for (size_t t = 1; t < pc.size(); ++t) {
                u32 k = pc[t].first;
                acc[k] = (acc[k] + na * pc[t].second) % P;
                if (!mark[k]) { mark[k] = 1; h.push(k); }
            }
            if (track_combos) {
                for (auto &e : cco[r]) {
                    u32 k = e.first;
                    cacc[k] = (cacc[k] + na * e.second) % P;
                    if (!cmark[k]) { cmark[k] = 1; ctouch.push_back(k); }
                }
            }
        }
        if (j < nbase) {
            if (fresh) ++base_rank;
            if (!fresh && track_combos) for (u32 k : ctouch) { cacc[k] = 0; cmark[k] = 0; }
            continue;
        }
        if (fresh) {
            if (dual_phi) {
                // Reverse back-substitution over pivot rows i < fresh_row:
                std::vector<u64> phi(nrows, 0);
                phi[fresh_row] = pw(fresh_val, P - 2);
                for (int i = (int)fresh_row - 1; i >= 0; --i) {
                    int r = piv[i];
                    if (r < 0) continue;
                    const SV &pc = cols[r];
                    u64 s = 0;
                    for (size_t t = 1; t < pc.size(); ++t) {
                        u32 k = pc[t].first;
                        if (k > fresh_row) break; // phi[k] == 0 for all k > fresh_row!
                        if (phi[k]) s = (s + phi[k] * pc[t].second) % P;
                    }
                    if (s) phi[i] = (P - s) % P;
                }
                snprintf(line, sizeof line, "F %u", fresh_row);
                out += line;
                for (u32 i = 0; i <= fresh_row; ++i) {
                    if (phi[i]) {
                        snprintf(line, sizeof line, " %u %llu", i, (unsigned long long)phi[i]);
                        out += line;
                    }
                }
                out += "\n";
            } else {
                out += "F\n";
            }
            continue;
        }
        out += "D";
        if (track_combos) {
            std::sort(ctouch.begin(), ctouch.end());
            for (u32 k : ctouch) {
                if (cacc[k]) { snprintf(line, sizeof line, " %u %llu", k, (unsigned long long)cacc[k]); out += line; }
                cacc[k] = 0; cmark[k] = 0;
            }
        }
        out += "\n";
    }
    printf("rank %u\n", base_rank);
    fwrite(out.data(), 1, out.size(), stdout);
    printf("total_rank %zu\n", cols.size());
    return 0;
}
