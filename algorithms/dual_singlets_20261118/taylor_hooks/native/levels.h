// levels.h -- per-level HW bases (A-H generator selection), Q0 rows, cohomology data (kernel, Cartan restriction).
#pragma once
#include "trace.h"
#include "batch.h"
#include "linalg.h"

static int MARGIN = 3; static bool HWCHECK = false, KEEPRAW = false, SCALAR = false; static int VERB = 1;
// evaluate G at probes 0..K-1 of pr (batched in lanes, threads over chunks)
static void eval_all(const Sector &sec, const Skel &S, const Gen &G, const Probes &pr, int K, std::vector<u64> &v) {
  v.assign(K, 0);
  if (SCALAR) { par_for((size_t)K, [&](size_t k) { v[k] = eval_at(sec, S, G, pr, (int)k); }); return; }
  size_t nch = (K + LN - 1) / LN;
  par_for(nch, [&](size_t c) { int ks[LN]; for (int l = 0; l < LN; l++) { int k = (int)(c * LN + l); ks[l] = k < K ? k : -1; } u64 out[LN]; b_eval_at(sec, S, G, pr, ks, out); for (int l = 0; l < LN; l++) if (ks[l] >= 0) v[ks[l]] = out[l]; });
}
static void build_skel_basis(const Sector &sec, Skel &S) {
  S.K = (int)S.hw + MARGIN; add_probes(S, S.pr, S.K, false);
  IncEch E; E.N = S.K; size_t fails = 0, genfails = 0;
  while ((i64)E.R.size() < S.hw) {
    Gen G; if (!random_gen(sec, S, G)) { if (++genfails > 10000000) { fprintf(stderr, "no generator for %s\n", S.key.c_str()); exit(8); } continue; }
    S.tried++;
    std::vector<u64> v; eval_all(sec, S, G, S.pr, S.K, v);
    std::vector<u64> v0; if (KEEPRAW) v0 = v;
    if (E.insert(v)) { S.gens.push_back(G); fails = 0; if (KEEPRAW) S.raw.push_back(to32(v0)); }
    else if (++fails > (size_t)std::max<i64>(5000, 100 * S.hw)) { fprintf(stderr, "stall %s rank %zu / %lld after %zu tries\n", S.key.c_str(), E.R.size(), (long long)S.hw, S.tried); exit(9); }
  }
  S.piv = E.piv;
  if (HWCHECK) {   // E12, E23, J+ kill every selected generator (2 fresh probes, t-substitution)
    Probes T; add_probes(S, T, 2, false);
    for (auto &G : S.gens) for (int o = 0; o < 3; o++) for (int p = 0; p < 2; p++) {
      u64 v = eval_gen(sec, S, G, [&](int l, GMat &out) {
        letter_at(S, T, p, l, out); RTerm rt = RAISE[o][l]; if (rt.a < 0) return;
        GMat ex; letter_at(S, T, p, rt.a, ex); for (auto &e : ex) { GE x; x.m = e.m | F_T; x.a = e.a; mscale(x.a, rt.c); out.push_back(x); } }, S.top_te | F_T);
      if (v) { fprintf(stderr, "HW check failed: skeleton %s op %d\n", S.key.c_str(), o); exit(10); }
    }
  }
}
struct Lev { int Y = 0; std::vector<Skel> sk; size_t dim = 0; i64 exact = 0; };
static void build_lev(const Sector &sec, Lev &L, int Y) {
  L.Y = Y; L.sk.clear(); L.dim = 0; L.exact = 0;
  if (Y < 1 || Y > sec.Ymax) return;
  double t0 = now();
  for (auto &kv : hw_by_skel(sec, Y)) {
    Skel S = parse_skel(kv.first, Y, kv.second); L.exact += kv.second;
    build_skel_basis(sec, S);
    L.dim += S.gens.size(); L.sk.push_back(std::move(S));
  }
  if ((i64)L.dim != L.exact) { fprintf(stderr, "level dim %zu != exact %lld\n", L.dim, (long long)L.exact); exit(12); }
  if (VERB) fprintf(stderr, "#level [%d,%d,%d,%d,%d] Y=%d: %zu skeletons, dim %zu, %.1fs, evals %zu, rss %.0fMB\n", sec.WS[0], sec.WS[1], sec.WS[2], sec.WS[3], sec.WS[4], Y, L.sk.size(), L.dim, now() - t0, (size_t)NEVAL, maxrss_mb());
}
static std::vector<size_t> col_offsets(const Lev &A, size_t &nc) {
  std::vector<size_t> c0(A.sk.size()); nc = 0;
  for (size_t i = 0; i < A.sk.size(); i++) { c0[i] = nc; nc += A.sk[i].gens.size(); }
  return c0;
}
static bool merge_groups(const Skel &S, int gx, const Skel &T, std::vector<std::pair<int, int>> &gab) {
  std::map<std::pair<int, int>, int> m;
  for (auto &g : T.g) m[{g.type, g.N}] += g.k;
  for (size_t i = 0; i < S.g.size(); i++) m[{S.g[i].type, S.g[i].N}] -= S.g[i].k - ((int)i == gx);
  gab.clear();
  for (auto &kv : m) { if (kv.second < 0) return false; for (int q = 0; q < kv.second; q++) gab.push_back(kv.first); }
  std::sort(gab.begin(), gab.end());
  return gab.size() == 2;
}
// Q0-substituted letter (x + eps q(x)) of target skeleton T at probe k, for source skeleton S with channels chans
typedef std::vector<std::pair<int, std::vector<std::pair<int, int>>>> Chans;
static void q0_letter(const Skel &S, const Skel &T, const Chans &chans, int k, int l, std::map<int, GMat> &qcache, GMat &out) {
  letter_at(T, T.pr, k, l, out);
  int gi = S.lgrp[l];
  for (auto &ch : chans) {
    if (gi != ch.first) continue;
    auto it = qcache.find(l);
    if (it == qcache.end()) {
      GMat acc;
      for (auto &qt : QB[l]) {
        std::vector<std::pair<int, int>> pr2 = {{LT[qt.a].type, LT[qt.a].r + LT[qt.a].s}, {LT[qt.b].type, LT[qt.b].r + LT[qt.b].s}}; std::sort(pr2.begin(), pr2.end());
        if (pr2 != ch.second) continue;
        GMat Xa, Xb, AB, BA; letter_at(T, T.pr, k, qt.a, Xa); letter_at(T, T.pr, k, qt.b, Xb);
        gmul(Xa, Xb, AB); gmul(Xb, Xa, BA);
        bool oo = LT[qt.a].odd && LT[qt.b].odd;   // [A,B} = AB - (-1)^{|a||b|} BA
        for (auto &x : AB) { GE e; e.m = x.m | 1ULL; e.a = x.a; mscale(e.a, qt.c); acc.push_back(e); }
        for (auto &x : BA) { GE e; e.m = x.m | 1ULL; e.a = x.a; mscale(e.a, oo ? qt.c : negm(qt.c)); acc.push_back(e); }
      }
      gcompact(acc);
      it = qcache.emplace(l, acc).first;
    }
    for (auto &x : it->second) out.push_back(x);
  }
}
// rows of Q0: A(Y) -> B(Y+1) at the pivot probes of B; row length = dim A
static Mtx q0_rows(const Sector &sec, const Lev &A, const Lev &B) {
  size_t nc; std::vector<size_t> c0 = col_offsets(A, nc);
  std::vector<std::pair<int, int>> rowid; for (size_t ti = 0; ti < B.sk.size(); ti++) for (int k : B.sk[ti].piv) rowid.push_back({(int)ti, k});
  Mtx rows(rowid.size(), std::vector<u32>(nc, 0));
  // tasks: (target skeleton, chunk of <= LN consecutive rows of that skeleton)
  struct Task { int ti; size_t r0, n; };
  std::vector<Task> tasks;
  { size_t r = 0; for (size_t ti = 0; ti < B.sk.size(); ti++) { size_t n = B.sk[ti].piv.size(), step = SCALAR ? 1 : LN; for (size_t o = 0; o < n; o += step) tasks.push_back({(int)ti, r + o, std::min(step, n - o)}); r += n; } }
  par_for(tasks.size(), [&](size_t t) {
    const Skel &T = B.sk[tasks[t].ti];
    for (size_t si = 0; si < A.sk.size(); si++) {
      const Skel &S = A.sk[si];
      Chans chans; std::vector<std::pair<int, int>> gab;
      for (size_t gx = 0; gx < S.g.size(); gx++) if (merge_groups(S, (int)gx, T, gab)) chans.push_back({(int)gx, gab});
      if (chans.empty()) continue;
      if (SCALAR) {
        int k = rowid[tasks[t].r0].second; std::map<int, GMat> qc;
        Ctx ctx = [&](int l, GMat &out) { q0_letter(S, T, chans, k, l, qc, out); };
        for (size_t i = 0; i < S.gens.size(); i++) rows[tasks[t].r0][c0[si] + i] = (u32)eval_gen(sec, S, S.gens[i], ctx, T.top_te | 1ULL);
        continue;
      }
      std::vector<std::map<int, GMat>> qc(LN); std::map<int, BGM> lcache;
      BCtx ctx = [&](int l, BGM &out) {
        auto it = lcache.find(l); if (it != lcache.end()) { out = it->second; return; }
        std::vector<GMat> lanes(LN);
        for (size_t ln = 0; ln < tasks[t].n; ln++) q0_letter(S, T, chans, rowid[tasks[t].r0 + ln].second, l, qc[ln], lanes[ln]);
        b_from_lanes(lanes, out); lcache[l] = out;
      };
      u64 val[LN];
      for (size_t i = 0; i < S.gens.size(); i++) { b_eval_gen(sec, S, S.gens[i], ctx, T.top_te | 1ULL, val); for (size_t ln = 0; ln < tasks[t].n; ln++) rows[tasks[t].r0 + ln][c0[si] + i] = (u32)val[ln]; }
    }
  });
  return rows;
}
// cohomology data at one level
struct Cohom {
  Sector sec; int Y = 0;
  Lev Lm, L0, Lp;
  size_t dim = 0, rank_in = 0, rank_out = 0; long h = 0;
  Mtx K;                          // kernel of the out-map (cocycles), vectors over the basis of L0
  Mtx Min;                        // in-map rows (pivot probes of L0 x basis of Lm): columns span B
  std::vector<size_t> c0; size_t nc = 0;
  struct RB { int si; Mtx val; }; // diagonal-probe values of the basis generators of skeleton si
  std::vector<RB> res;
  size_t rho = 0; std::vector<std::pair<int, int>> witness; Mtx resK;
};
static void compute_res(Cohom &C) {
  C.res.clear();
  for (size_t si = 0; si < C.L0.sk.size(); si++) {
    Skel &S = C.L0.sk[si]; size_t nb = S.gens.size(); if (!nb) continue;
    int nd = (int)nb + MARGIN; S.dg.clear(); add_probes(S, S.dg, nd, true);
    Cohom::RB rb; rb.si = (int)si; rb.val.assign(nd, std::vector<u32>(nb, 0));
    if (SCALAR) par_for((size_t)nd * nb, [&](size_t t) { int k = (int)(t / nb); size_t i = t % nb; rb.val[k][i] = (u32)eval_at(C.sec, S, S.gens[i], S.dg, k); });
    else { size_t nch = (nd + LN - 1) / LN; par_for(nch * nb, [&](size_t t) { size_t c = t / nb, i = t % nb; int ks[LN]; for (int l = 0; l < LN; l++) { int k = (int)(c * LN + l); ks[l] = k < nd ? k : -1; } u64 out[LN]; b_eval_at(C.sec, S, S.gens[i], S.dg, ks, out); for (int l = 0; l < LN; l++) if (ks[l] >= 0) rb.val[ks[l]][i] = (u32)out[l]; }); }
    C.res.push_back(std::move(rb));
  }
}
// want_in: in-map rank; keep_min: keep the in-map rows (for exact-in-span tests); want_res: Cartan restriction
static void compute_cohom(Cohom &C, const Sector &sec, int Y, bool want_in, bool keep_min, bool want_res) {
  C.sec = sec; C.Y = Y; double t0 = now();
  build_lev(sec, C.Lm, want_in ? Y - 1 : 0); build_lev(sec, C.L0, Y); build_lev(sec, C.Lp, Y + 1);
  C.dim = C.L0.dim; C.c0 = col_offsets(C.L0, C.nc);
  Mtx Mout = q0_rows(sec, C.L0, C.Lp);
  C.K = kernel_of(Mout, C.nc); C.rank_out = C.nc - C.K.size();
  if (want_in && C.Lm.dim) { Mtx Min = q0_rows(sec, C.Lm, C.L0); size_t n2; col_offsets(C.Lm, n2); if (keep_min) C.Min = Min; C.rank_in = rank_of(Min, n2); }
  C.h = (long)C.dim - (long)C.rank_in - (long)C.rank_out;
  if (VERB) fprintf(stderr, "#Q0 [%d,%d,%d,%d,%d] Y=%d: dim %zu in %zu out %zu kernel %zu h %ld  %.1fs\n", sec.WS[0], sec.WS[1], sec.WS[2], sec.WS[3], sec.WS[4], Y, C.dim, C.rank_in, C.rank_out, C.K.size(), C.h, now() - t0);
  if (!want_res) return;
  compute_res(C);
  IncEch E; E.N = (int)C.K.size(); C.witness.clear(); C.resK.clear();
  for (auto &rb : C.res) {
    size_t o = C.c0[rb.si], nb = C.L0.sk[rb.si].gens.size();
    for (size_t k = 0; k < rb.val.size(); k++) {
      std::vector<u64> row(C.K.size(), 0);
      for (size_t z = 0; z < C.K.size(); z++) { u64 s = 0; int c = 0; for (size_t i = 0; i < nb; i++) { s += (u64)rb.val[k][i] * C.K[z][o + i]; if (++c == 3) { s = red(s); c = 0; } } row[z] = red(s); }
      std::vector<u64> rc = row;
      if (E.insert(rc)) { C.witness.push_back({rb.si, (int)k}); C.resK.push_back(to32(row)); }
    }
  }
  C.rho = C.witness.size();
  if (VERB) fprintf(stderr, "#res Y=%d: rho %zu (h %ld)  %.1fs\n", Y, C.rho, C.h, now() - t0);
}
