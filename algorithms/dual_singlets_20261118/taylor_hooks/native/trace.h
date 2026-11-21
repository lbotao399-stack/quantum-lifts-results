// trace.h -- letters, sectors, per-skeleton exact HW dimensions, Grassmann matrices, A-H trace-word generators,
// probes, Q0 letter table.  Ported from su3_engine/design/A_trace_eval/proto_AH2.cpp and su3_engine/q1/x7/q1ah.cpp
// (copies in ../ref/), with the colour algebra replaced by core.h (any N_F <= 7) and without the charge-conjugation
// split (C is inner for so(2n+1), sp(n); for SU(N) the engine reports the C-summed numbers).
#pragma once
#include "core.h"

// ---------------- letters (union box of all sectors of the run) ----------------
enum { GAM = 0, BET = 1, BF = 2, CG = 3 };
struct CL { int type, flav, r, s, odd; int w[5]; };
static std::vector<CL> LT; static std::map<std::array<int, 4>, int> LTID; static int BOX[5];
static void build_letters() {
  LT.clear(); LTID.clear();
  for (int t = 0; t < 4; t++) {
    int nf = t < 2 ? 3 : 1;
    for (int r = 0; r <= BOX[3] + 1; r++) for (int s = 0; s <= BOX[4] + 1; s++) {
      if (t == CG && r + s == 0) continue;
      for (int f = 0; f < nf; f++) {
        CL x; x.type = t; x.flav = f; x.r = r; x.s = s; x.odd = (t == BET || t == CG);
        for (int i = 0; i < 3; i++) x.w[i] = t == GAM ? (i == f) : t == BET ? (i != f) : t == BF ? 1 : 0;
        x.w[3] = r; x.w[4] = s;
        bool ok = true; for (int i = 0; i < 3; i++) if (x.w[i] > BOX[i] + 2) ok = false;
        if (!ok) continue;
        LTID[{t, f, r, s}] = (int)LT.size(); LT.push_back(x);
      }
    }
  }
}
static int ltid(int t, int f, int r, int s) { if (r < 0 || s < 0) return -1; if (t == CG && r + s == 0) return -1; auto it = LTID.find({t, f, r, s}); return it == LTID.end() ? -1 : it->second; }

// ---------------- sector data ----------------
struct Sector { int WS[5]; int COLS = 0, CLEN[16], ETAB[16][3], NBR = 0; u64 TOP_EZ = 0; int Theta = 0, n = 0, j = 0, Ymax = 0; };
static Sector make_sector(const int w[5]) {
  Sector s; for (int i = 0; i < 5; i++) s.WS[i] = w[i];
  s.COLS = w[0]; int eb = 24;
  if (w[0] > 16) { fprintf(stderr, "too many flavour columns\n"); exit(4); }
  for (int j = 0; j < s.COLS; j++) { s.CLEN[j] = (w[0] > j) + (w[1] > j) + (w[2] > j); for (int f = 0; f < s.CLEN[j]; f++) s.ETAB[j][f] = eb++; }
  s.NBR = w[4];
  if (eb > 48 || 48 + 2 * s.NBR > 64) { fprintf(stderr, "eta/zeta family overflow\n"); exit(4); }
  for (int b = 24; b < eb; b++) s.TOP_EZ |= 1ULL << b;
  for (int e = 0; e < s.NBR; e++) s.TOP_EZ |= 3ULL << (48 + 2 * e);
  s.Theta = w[0] + w[1] + w[2]; s.n = w[3] + w[4]; s.j = 2 * s.Theta + 3 * s.n; s.Ymax = s.Theta + s.n;
  return s;
}

// ---------------- per-skeleton HW counts (flavour SU(3) x dotted SL(2) alternation of singlet counts) ----------------
static std::string skel_of(const std::vector<std::pair<int, int>> &L) {
  std::map<std::pair<int, int>, int> c; for (auto &x : L) c[{LT[x.first].type, LT[x.first].r + LT[x.first].s}] += x.second;
  std::string s; for (auto &kv : c) { s += char('0' + kv.first.first); s += char('A' + kv.first.second); s += char('a' + kv.second); } return s;
}
static void skel_singlets(const int mu[5], int Y, int sg, std::map<std::string, i64> &acc) {
  for (int i = 0; i < 5; i++) if (mu[i] < 0) return;
  std::vector<std::pair<int, int>> cur; int w[5] = {0, 0, 0, 0, 0};
  std::function<void(size_t, int)> rec = [&](size_t i, int left) {
    if (!left) { for (int q = 0; q < 5; q++) if (w[q] != mu[q]) return; std::vector<std::pair<int, int>> pat; for (auto &x : cur) pat.push_back({x.second, LT[x.first].odd}); i64 d = pattern_sing(pat); if (d) acc[skel_of(cur)] += sg * d; return; }
    if (i == LT.size()) return;
    const CL &x = LT[i]; int mk = left; if (x.odd) mk = std::min(mk, DIMG);
    for (int k = mk; k >= 0; k--) {
      bool ok = true; for (int q = 0; q < 5; q++) if (w[q] + k * x.w[q] > mu[q]) ok = false;
      if (!ok) continue;
      for (int q = 0; q < 5; q++) w[q] += k * x.w[q];
      if (k) cur.push_back({(int)i, k});
      rec(i + 1, left - k);
      if (k) cur.pop_back();
      for (int q = 0; q < 5; q++) w[q] -= k * x.w[q];
    }
  };
  rec(0, Y);
}
static std::map<std::string, i64> hw_by_skel(const Sector &sec, int Y) {
  static const int T3[6][4] = {{0,0,0,1},{1,-1,0,-1},{0,1,-1,-1},{2,0,-2,-1},{1,1,-2,1},{2,-1,-1,1}};
  std::map<std::string, i64> out;
  for (int dd = 0; dd < 2; dd++) for (auto &t : T3) { int mu[5] = {sec.WS[0] + t[0], sec.WS[1] + t[1], sec.WS[2] + t[2], sec.WS[3] + dd, sec.WS[4] - dd}; skel_singlets(mu, Y, t[3] * (dd ? -1 : 1), out); }
  for (auto it = out.begin(); it != out.end();) { if (it->second < 0) { fprintf(stderr, "negative skeleton HW count %s\n", it->first.c_str()); exit(15); } if (!it->second) it = out.erase(it); else ++it; }
  return out;
}

// ---------------- Q0 letter table (engine convention of v6 / proto_AH2) ----------------
struct QTerm { int a, b; u64 c; };
static std::vector<std::vector<QTerm>> QB;
static u64 binom(int n, int k) { if (k < 0 || k > n) return 0; u64 r = 1; for (int i = 1; i <= k; i++) r = r * (n - k + i) / i; return r; }
static void build_q() {
  QB.assign(LT.size(), {}); u64 half = invm(2);
  for (size_t id = 0; id < LT.size(); id++) {
    const CL &x = LT[id]; int r = x.r, s = x.s;
    auto br = [&](int ta, int fa, int tb, int fb, u64 coef) {
      for (int u = 0; u <= r; u++) for (int v = 0; v <= s; v++) {
        int a = ltid(ta, fa, u, v), b = ltid(tb, fb, r - u, s - v); if (a < 0 || b < 0) continue;
        QB[id].push_back({a, b, mulm(coef, mulm(binom(r, u) % P, binom(s, v) % P))});
      }
    };
    if (x.type == CG) br(CG, 0, CG, 0, mulm(P - 1, half));                                 // -1/2 [c,c]
    else if (x.type == GAM) br(CG, 0, GAM, x.flav, P - 1);                                 // -[c,gamma^I]
    else if (x.type == BET) { br(CG, 0, BET, x.flav, P - 1); br(GAM, (x.flav + 1) % 3, GAM, (x.flav + 2) % 3, 1); }   // -[c,beta_I] + [g^J,g^K]
    else { br(CG, 0, BF, 0, P - 1); for (int J = 0; J < 3; J++) br(GAM, J, BET, J, P - 1); }   // -[c,b] - [gamma^J,beta_J]
  }
}
struct RTerm { int a; u64 c; };
static std::vector<RTerm> RAISE[3];   // E12, E23, J+ on letters
static void build_raise() {
  for (int o = 0; o < 3; o++) RAISE[o].assign(LT.size(), {-1, 0});
  for (size_t id = 0; id < LT.size(); id++) {
    const CL &x = LT[id];
    if (x.type == GAM && x.flav == 1) RAISE[0][id] = {ltid(GAM, 0, x.r, x.s), 1};
    if (x.type == BET && x.flav == 0) RAISE[0][id] = {ltid(BET, 1, x.r, x.s), P - 1};
    if (x.type == GAM && x.flav == 2) RAISE[1][id] = {ltid(GAM, 1, x.r, x.s), 1};
    if (x.type == BET && x.flav == 1) RAISE[1][id] = {ltid(BET, 2, x.r, x.s), P - 1};
    if (x.s > 0) RAISE[2][id] = {ltid(x.type, x.flav, x.r + 1, x.s - 1), (u64)x.s % P};
  }
}

// ---------------- family-graded Grassmann matrices ----------------
// bits 0..22 theta/eps (bit 0 = eps of Q0), bit 23 t (even nilpotent), bits 24..47 eta, bits 48..63 zeta;
// the three odd families commute with each other (sign computed per family)
static const u64 F_TE = (1ULL << 23) - 1, F_T = 1ULL << 23, F_ETA = ((1ULL << 48) - 1) & ~((1ULL << 24) - 1), F_Z = ~((1ULL << 48) - 1);
static inline int fsign(u64 S, u64 T) {   // sign of S*T -> sorted(S|T), per family
  int c = 0; const u64 fam[3] = {F_TE, F_ETA, F_Z};
  for (u64 F : fam) { u64 So = S & F, To = T & F; while (To) { int b = __builtin_ctzll(To); To &= To - 1; c += (b == 63) ? 0 : __builtin_popcountll(So >> (b + 1)); } }
  return c & 1;
}
struct GE { u64 m; Mat a; };
typedef std::vector<GE> GMat;
typedef std::vector<std::pair<u64, u64>> GSc;
static void gcompact(GMat &C) {
  if (C.size() < 2) { if (C.size() == 1 && miszero(C[0].a)) C.clear(); return; }
  std::vector<std::pair<u64, u32>> key(C.size()); for (size_t i = 0; i < C.size(); i++) key[i] = {C[i].m, (u32)i};
  std::sort(key.begin(), key.end());
  GMat out; out.reserve(C.size());
  for (size_t i = 0; i < key.size();) {
    size_t j = i; GE e = C[key[i].second]; j++;
    while (j < key.size() && key[j].first == key[i].first) { maddto(e.a, C[key[j].second].a); j++; }
    if (!miszero(e.a)) out.push_back(e);
    i = j;
  }
  C.swap(out);
}
static void gmul(const GMat &A, const GMat &B, GMat &C) {
  C.clear(); C.reserve(A.size() * B.size());
  for (auto &a : A) for (auto &b : B) {
    if (a.m & b.m) continue;
    GE e; e.m = a.m | b.m; mmul(a.a, b.a, e.a);
    if (fsign(a.m, b.m)) mneg(e.a);
    C.push_back(e);
  }
  gcompact(C);
}
static void scompact(GSc &C) {
  std::sort(C.begin(), C.end()); size_t w = 0;
  for (size_t i = 0; i < C.size(); i++) { if (w && C[w - 1].first == C[i].first) C[w - 1].second = addm(C[w - 1].second, C[i].second); else C[w++] = C[i]; }
  C.resize(w); size_t z = 0; for (size_t i = 0; i < C.size(); i++) if (C[i].second) C[z++] = C[i]; C.resize(z);
}
static void gtrmul(const GMat &A, const GMat &B, GSc &s) {   // tr(A B) as a Grassmann scalar
  s.clear();
  for (auto &a : A) for (auto &b : B) { if (a.m & b.m) continue; u64 v = mtrmul(a.a, b.a); if (!v) continue; if (fsign(a.m, b.m)) v = negm(v); s.push_back({a.m | b.m, v}); }
  scompact(s);
}
static void gtr(const GMat &A, GSc &s) { s.clear(); for (auto &a : A) { u64 v = mtrace(a.a); if (v) s.push_back({a.m, v}); } }
static void gsmul(const GSc &A, const GSc &B, GSc &C) {
  C.clear(); C.reserve(A.size() * B.size());
  for (auto &a : A) for (auto &b : B) { if (a.first & b.first) continue; u64 v = mulm(a.second, b.second); if (fsign(a.first, b.first)) v = negm(v); C.push_back({a.first | b.first, v}); }
  scompact(C);
}
static void sgmul(const GSc &A, const GMat &B, GMat &C) {   // scalar * matrix (order kept)
  C.clear();
  for (auto &a : A) for (auto &b : B) { if (a.first & b.m) continue; GE e; e.m = a.first | b.m; e.a = b.a; u64 c = fsign(a.first, b.m) ? negm(a.second) : a.second; mscale(e.a, c); C.push_back(e); }
  gcompact(C);
}
static void gmul_const(const GMat &A, const Mat &M, bool left, GMat &C) {   // M*A (left) or A*M, M even constant
  C.clear(); for (auto &a : A) { GE e; e.m = a.m; if (left) mmul(M, a.a, e.a); else mmul(a.a, M, e.a); C.push_back(e); }
  gcompact(C);
}
static void gtransp(const GMat &A, GMat &T) { T.clear(); for (auto &a : A) { GE e; e.m = a.m; mtranspose(a.a, e.a); T.push_back(e); } }
static void gadd(GMat &A, const GMat &B, bool neg = false) { for (auto &b : B) { GE e = b; if (neg) mneg(e.a); A.push_back(e); } gcompact(A); }
static void gsadd(GSc &A, const GSc &B, u64 c = 1) { for (auto &b : B) A.push_back({b.first, mulm(b.second, c)}); scompact(A); }
// product of the traces of the cycles, top coefficient
static u64 eval_traces(const std::vector<std::vector<int>> &tr, const std::vector<GMat> &pos, u64 top) {
  GSc tot = {{0, 1}}, s, ns; GMat cur, nx;
  for (auto &t : tr) {
    if (t.size() == 1) { gtr(pos[t[0]], s); }
    else {
      cur = pos[t[0]];
      for (size_t i = 1; i + 1 < t.size(); i++) { gmul(cur, pos[t[i]], nx); std::swap(cur, nx); if (cur.empty()) return 0; }
      gtrmul(cur, pos[t.back()], s);
    }
    if (s.empty()) return 0;
    gsmul(tot, s, ns); std::swap(tot, ns); if (tot.empty()) return 0;
  }
  for (auto &x : tot) if (x.first == top) return x.second;
  return 0;
}

// ---------------- skeletons, probes, generators ----------------
static const int THDEG[4] = {1, 2, 3, 0};
struct Grp { int type, N, k; };
struct Gen { std::vector<std::vector<int>> cols, edges; std::vector<int> mu; std::vector<std::vector<int>> tr; };
typedef std::vector<std::vector<std::vector<Mat>>> Probes;   // [probe][letter] -> matrices
struct IncEch {   // incremental reduced echelon of length-N vectors (u32 storage); rows kept fully reduced (RREF)
  int N = 0; std::vector<std::vector<u32>> R; std::vector<int> piv;
  bool insert(std::vector<u64> v) {
    size_t r = R.size();
    std::vector<u64> a(r); for (size_t k = 0; k < r; k++) a[k] = v[piv[k]];
    size_t nchunk = (N >= 512 && r >= 64) ? (size_t)std::min(64, (N + 255) / 256) : 1, cw = (N + nchunk - 1) / nchunk;
    par_for(nchunk, [&](size_t c) {
      size_t j0 = c * cw, j1 = std::min((size_t)N, j0 + cw);
      for (size_t k = 0; k < r; k++) { if (!a[k]) continue; u64 na = P - a[k]; const u32 *row = R[k].data(); for (size_t j = j0; j < j1; j++) if (row[j]) v[j] = red(v[j] + na * row[j]); }
    });
    int p = -1; for (int j = 0; j < N; j++) if (v[j]) { p = j; break; }
    if (p < 0) return false;
    u64 inv = invm(v[p]); std::vector<u32> w(N); for (int j = 0; j < N; j++) w[j] = (u32)mulm(v[j], inv);
    std::vector<u32> nz; for (int j = 0; j < N; j++) if (w[j]) nz.push_back((u32)j);
    par_for(r, [&](size_t k) { u64 c = R[k][p]; if (!c) return; u64 nc = P - c; u32 *row = R[k].data(); for (u32 j : nz) row[j] = (u32)red(row[j] + nc * w[j]); });
    R.push_back(std::move(w)); piv.push_back(p); return true;
  }
};
struct Skel {
  std::string key; int Y = 0; i64 hw = 0;
  std::vector<Grp> g; std::vector<int> pg, thbase; u64 top_te = 0;
  std::vector<int> lgrp;
  int K = 0; Probes pr, dg;
  std::vector<Gen> gens;            // selected generators = basis (all independent)
  std::vector<int> piv;             // pivot probes (injective coordinates on the span)
  std::vector<std::vector<u32>> raw, Einv;   // optional: raw values at all K probes; inverse of the pivot evaluation matrix
  size_t tried = 0;
};
static Skel parse_skel(const std::string &key, int Y, i64 hw) {
  Skel S; S.key = key; S.Y = Y; S.hw = hw;
  for (size_t i = 0; i < key.size(); i += 3) S.g.push_back({key[i] - '0', key[i + 1] - 'A', key[i + 2] - 'a'});
  for (size_t i = 0; i < S.g.size(); i++) for (int q = 0; q < S.g[i].k; q++) S.pg.push_back((int)i);
  int b = 1; S.thbase.assign(S.g.size(), -1);
  for (size_t i = 0; i < S.g.size(); i++) if (S.g[i].type == BET || S.g[i].type == CG) { S.thbase[i] = b; b += S.g[i].k; }
  if (b > 23) { fprintf(stderr, "theta family overflow\n"); exit(4); }
  for (int q = 1; q < b; q++) S.top_te |= 1ULL << q;
  S.lgrp.assign(LT.size(), -1);
  for (size_t l = 0; l < LT.size(); l++) for (size_t i = 0; i < S.g.size(); i++) if (S.g[i].type == LT[l].type && S.g[i].N == LT[l].r + LT[l].s) S.lgrp[l] = (int)i;
  return S;
}
static void add_probes(const Skel &S, Probes &pr, int n, bool diag) {
  for (int p = 0; p < n; p++) {
    std::vector<std::vector<Mat>> one(LT.size());
    for (size_t l = 0; l < LT.size(); l++) { int gi = S.lgrp[l]; if (gi < 0) continue; int m = LT[l].odd ? S.g[gi].k : 1; for (int r = 0; r < m; r++) one[l].push_back(diag ? rand_cartan() : rand_alg()); }
    pr.push_back(std::move(one));
  }
}
static void letter_at(const Skel &S, const Probes &pr, int p, int l, GMat &out) {
  out.clear(); if (l < 0) return; const auto &M = pr[p][l]; if (M.empty()) return;
  if (!LT[l].odd) { out.push_back({0, M[0]}); return; }
  int gi = S.lgrp[l];
  for (size_t r = 0; r < M.size(); r++) out.push_back({1ULL << (S.thbase[gi] + r), M[r]});
}
static bool random_gen(const Sector &sec, const Skel &S, Gen &G) {
  int Y = S.Y; G.cols.assign(Y, {}); G.edges.assign(Y, {}); G.mu.assign(Y, 0); G.tr.clear();
  std::vector<int> cl; for (int j = 0; j < sec.COLS; j++) for (int i = 0; i < sec.CLEN[j]; i++) cl.push_back(j);
  for (size_t i = cl.size(); i > 1; i--) std::swap(cl[i - 1], cl[rnd() % i]);
  size_t c = 0;
  for (int p = 0; p < Y; p++) for (int q = 0; q < THDEG[S.g[S.pg[p]].type]; q++) G.cols[p].push_back(cl[c++]);
  std::vector<int> cap(Y); for (int p = 0; p < Y; p++) cap[p] = S.g[S.pg[p]].N;
  for (int e = 0; e < sec.NBR; e++) {
    std::vector<int> av; for (int p = 0; p < Y; p++) if (cap[p]) av.push_back(p);
    if (av.size() < 2) return false;
    int a = av[rnd() % av.size()], b; do b = av[rnd() % av.size()]; while (b == a);
    G.edges[a].push_back(e); G.edges[b].push_back(e); cap[a]--; cap[b]--;
  }
  for (int p = 0; p < Y; p++) G.mu[p] = cap[p];
  std::vector<int> pm(Y);
  while (true) { std::iota(pm.begin(), pm.end(), 0); for (int i = Y; i > 1; i--) std::swap(pm[i - 1], pm[rnd() % i]); bool fp = false; for (int i = 0; i < Y; i++) if (pm[i] == i) fp = true; if (!fp) break; }
  std::vector<int> seen(Y, 0);
  for (int i = 0; i < Y; i++) if (!seen[i]) { std::vector<int> cyc; int x = i; while (!seen[x]) { seen[x] = 1; cyc.push_back(x); x = pm[x]; } G.tr.push_back(cyc); }
  return true;
}
struct FT { u64 m; int sg; int flav; };
static void flavour_terms(const Sector &sec, int type, const std::vector<int> &cols, std::vector<FT> &out) {
  out.clear();
  if (type == CG) { out.push_back({0, 0, 0}); return; }
  if (type == GAM) { int j = cols[0]; for (int f = 0; f < sec.CLEN[j]; f++) out.push_back({1ULL << sec.ETAB[j][f], 0, f}); return; }
  if (type == BET) {
    int j1 = cols[0], j2 = cols[1];
    for (int f1 = 0; f1 < sec.CLEN[j1]; f1++) for (int f2 = 0; f2 < sec.CLEN[j2]; f2++) {
      if (f1 == f2) continue;
      u64 m1 = 1ULL << sec.ETAB[j1][f1], m2 = 1ULL << sec.ETAB[j2][f2]; if (m1 & m2) continue;
      int I = 3 - f1 - f2, ws = (f1 == (I + 1) % 3 && f2 == (I + 2) % 3) ? 0 : 1;
      out.push_back({m1 | m2, fsign(m1, m2) ^ ws, I});
    }
    return;
  }
  int j1 = cols[0], j2 = cols[1], j3 = cols[2];
  int perm[6][3] = {{0,1,2},{1,2,0},{2,0,1},{0,2,1},{2,1,0},{1,0,2}};
  for (int q = 0; q < 6; q++) {
    int f1 = perm[q][0], f2 = perm[q][1], f3 = perm[q][2];
    if (f1 >= sec.CLEN[j1] || f2 >= sec.CLEN[j2] || f3 >= sec.CLEN[j3]) continue;
    u64 m1 = 1ULL << sec.ETAB[j1][f1], m2 = 1ULL << sec.ETAB[j2][f2], m3 = 1ULL << sec.ETAB[j3][f3];
    if ((m1 & m2) || ((m1 | m2) & m3)) continue;
    out.push_back({m1 | m2 | m3, fsign(m1, m2) ^ fsign(m1 | m2, m3) ^ (q >= 3), 0});
  }
}
struct DT { u64 m; int sg; int r, s; };
static void dotted_terms(const std::vector<int> &edges, int mu, std::vector<DT> &out) {
  out.clear(); int ne = (int)edges.size();
  for (int c = 0; c < (1 << ne); c++) {
    u64 m = 0; int sg = 0, nu = 0;
    for (int i = 0; i < ne; i++) { u64 b = 1ULL << (48 + 2 * edges[i] + ((c >> i) & 1)); sg ^= fsign(m, b); m |= b; nu += !((c >> i) & 1); }
    out.push_back({m, sg, mu + nu, ne - nu});
  }
}
struct Lab { u64 m; int sg; int l; };
static void position_labels(const Sector &sec, int type, const Gen &G, int p, std::vector<Lab> &out) {
  out.clear(); std::vector<FT> ft; std::vector<DT> dt; flavour_terms(sec, type, G.cols[p], ft); dotted_terms(G.edges[p], G.mu[p], dt);
  for (auto &f : ft) for (auto &d : dt) { int l = ltid(type, f.flav, d.r, d.s); if (l < 0) continue; out.push_back({f.m | d.m, f.sg ^ d.sg, l}); }
}
typedef std::function<void(int, GMat &)> Ctx;
static void position_mat(const Sector &sec, int type, const Gen &G, int p, const Ctx &ctx, GMat &out) {
  out.clear(); std::vector<Lab> lab; position_labels(sec, type, G, p, lab);
  GMat lm;
  for (auto &a : lab) {
    ctx(a.l, lm);
    for (auto &e : lm) { GE x; x.m = a.m | e.m; x.a = e.a; if (a.sg) mneg(x.a); out.push_back(x); }
  }
  gcompact(out);
}
static std::atomic<size_t> NEVAL(0);
static u64 eval_gen(const Sector &sec, const Skel &S, const Gen &G, const Ctx &ctx, u64 top_te) {
  NEVAL++;
  std::vector<GMat> pos(S.Y);
  for (int p = 0; p < S.Y; p++) { position_mat(sec, S.g[S.pg[p]].type, G, p, ctx, pos[p]); if (pos[p].empty()) return 0; }
  return eval_traces(G.tr, pos, top_te | sec.TOP_EZ);
}
static u64 eval_at(const Sector &sec, const Skel &S, const Gen &G, const Probes &pr, int k) {
  return eval_gen(sec, S, G, [&](int l, GMat &o) { letter_at(S, pr, k, l, o); }, S.top_te);
}
