// batch.h -- lane-batched Grassmann matrix evaluation: one generator at LN probes at once (same mask structure,
// lane-wise matrix arithmetic).  Bit-for-bit the same arithmetic as the scalar path (trace.h), checked by mode batchcheck.
#pragma once
#include "trace.h"

static const int LN = 8;
struct BGM { std::vector<u64> m; std::vector<u32> d; size_t size() const { return m.size(); } void clear() { m.clear(); d.clear(); } };
struct BSc { std::vector<u64> m; std::vector<u64> v; size_t size() const { return m.size(); } void clear() { m.clear(); v.clear(); } };
static inline u32 *bblk(BGM &A, size_t e) { return &A.d[e * (size_t)NN * LN]; }
static inline const u32 *bblk(const BGM &A, size_t e) { return &A.d[e * (size_t)NN * LN]; }
static inline bool bzero(const u32 *x, size_t n) { for (size_t i = 0; i < n; i++) if (x[i]) return false; return true; }
static void bgcompact(BGM &C) {
  size_t n = C.m.size(), B = (size_t)NN * LN; if (!n) return;
  std::vector<std::pair<u64, u32>> key(n); for (size_t i = 0; i < n; i++) key[i] = {C.m[i], (u32)i};
  std::sort(key.begin(), key.end());
  BGM out; out.m.reserve(n); out.d.reserve(n * B);
  for (size_t i = 0; i < n;) {
    size_t j = i + 1; size_t base = out.d.size();
    const u32 *src = bblk(C, key[i].second); out.d.insert(out.d.end(), src, src + B);
    while (j < n && key[j].first == key[i].first) { const u32 *s2 = bblk(C, key[j].second); u32 *dst = &out.d[base]; for (size_t q = 0; q < B; q++) dst[q] = (u32)addm(dst[q], s2[q]); j++; }
    if (bzero(&out.d[base], B)) out.d.resize(base); else out.m.push_back(key[i].first);
    i = j;
  }
  C.m.swap(out.m); C.d.swap(out.d);
}
template <int N> static void bmmulT(const u32 *A, const u32 *B, u32 *C, bool neg) {   // lane-wise C = +-A B
  for (int i = 0; i < N; i++) for (int k = 0; k < N; k++) {
    u64 acc[LN];
    for (int l = 0; l < LN; l++) acc[l] = (u64)A[(i * N) * LN + l] * B[(k) * LN + l];
    for (int j = 1; j < N; j++) {
      const u32 *a = &A[(i * N + j) * LN], *b = &B[(j * N + k) * LN];
      for (int l = 0; l < LN; l++) acc[l] += (u64)a[l] * b[l];
      if (j % 3 == 2) for (int l = 0; l < LN; l++) acc[l] = MERS ? (acc[l] & 2147483647ULL) + (acc[l] >> 31) : acc[l] % P;
    }
    u32 *c = &C[(i * N + k) * LN];
    for (int l = 0; l < LN; l++) { u64 v = red(acc[l]); c[l] = (u32)(neg ? negm(v) : v); }
  }
}
template <int N> static void btrmulT(const u32 *A, const u32 *B, u64 *out, bool neg) {   // lane-wise tr(A B)
  u64 s[LN]; for (int l = 0; l < LN; l++) s[l] = 0;
  for (int i = 0; i < N; i++) {
    u64 t[LN]; for (int l = 0; l < LN; l++) t[l] = 0;
    for (int j = 0; j < N; j++) { const u32 *a = &A[(i * N + j) * LN], *b = &B[(j * N + i) * LN]; for (int l = 0; l < LN; l++) t[l] += (u64)a[l] * b[l]; if (j % 3 == 2) for (int l = 0; l < LN; l++) t[l] = MERS ? (t[l] & 2147483647ULL) + (t[l] >> 31) : t[l] % P; }
    for (int l = 0; l < LN; l++) s[l] += red(t[l]);
  }
  for (int l = 0; l < LN; l++) { u64 v = red(s[l]); out[l] = neg ? negm(v) : v; }
}
typedef void (*BMMulF)(const u32 *, const u32 *, u32 *, bool); typedef void (*BTrF)(const u32 *, const u32 *, u64 *, bool);
static BMMulF BMMUL = bmmulT<7>; static BTrF BTRMUL = btrmulT<7>;
static void set_bmm() {
  switch (NF) { case 2: BMMUL = bmmulT<2>; BTRMUL = btrmulT<2>; break; case 3: BMMUL = bmmulT<3>; BTRMUL = btrmulT<3>; break; case 4: BMMUL = bmmulT<4>; BTRMUL = btrmulT<4>; break;
    case 5: BMMUL = bmmulT<5>; BTRMUL = btrmulT<5>; break; case 6: BMMUL = bmmulT<6>; BTRMUL = btrmulT<6>; break; default: BMMUL = bmmulT<7>; BTRMUL = btrmulT<7>; }
}
static void bgmul(const BGM &A, const BGM &B, BGM &C) {
  C.clear(); size_t Bs = (size_t)NN * LN;
  size_t cnt = 0; for (u64 a : A.m) for (u64 b : B.m) if (!(a & b)) cnt++;
  C.m.reserve(cnt); C.d.resize(cnt * Bs);
  size_t e = 0;
  for (size_t i = 0; i < A.m.size(); i++) for (size_t j = 0; j < B.m.size(); j++) {
    u64 a = A.m[i], b = B.m[j]; if (a & b) continue;
    BMMUL(bblk(A, i), bblk(B, j), &C.d[e * Bs], fsign(a, b)); C.m.push_back(a | b); e++;
  }
  bgcompact(C);
}
static void bscompact(BSc &C) {
  size_t n = C.m.size(); if (!n) return;
  std::vector<std::pair<u64, u32>> key(n); for (size_t i = 0; i < n; i++) key[i] = {C.m[i], (u32)i};
  std::sort(key.begin(), key.end());
  BSc out; out.m.reserve(n); out.v.reserve(n * LN);
  for (size_t i = 0; i < n;) {
    size_t j = i + 1, base = out.v.size();
    for (int l = 0; l < LN; l++) out.v.push_back(C.v[key[i].second * LN + l]);
    while (j < n && key[j].first == key[i].first) { for (int l = 0; l < LN; l++) out.v[base + l] = addm(out.v[base + l], C.v[key[j].second * LN + l]); j++; }
    bool z = true; for (int l = 0; l < LN; l++) if (out.v[base + l]) z = false;
    if (z) out.v.resize(base); else out.m.push_back(key[i].first);
    i = j;
  }
  C.m.swap(out.m); C.v.swap(out.v);
}
static void bgtrmul(const BGM &A, const BGM &B, BSc &s) {
  s.clear();
  for (size_t i = 0; i < A.m.size(); i++) for (size_t j = 0; j < B.m.size(); j++) {
    u64 a = A.m[i], b = B.m[j]; if (a & b) continue;
    size_t o = s.v.size(); s.v.resize(o + LN); BTRMUL(bblk(A, i), bblk(B, j), &s.v[o], fsign(a, b)); s.m.push_back(a | b);
  }
  bscompact(s);
}
static void bgtr(const BGM &A, BSc &s) {
  s.clear();
  for (size_t i = 0; i < A.m.size(); i++) { const u32 *x = bblk(A, i); size_t o = s.v.size(); s.v.resize(o + LN); for (int l = 0; l < LN; l++) { u64 t = 0; for (int q = 0; q < NF; q++) t += x[(q * NF + q) * LN + l]; s.v[o + l] = red(t); } s.m.push_back(A.m[i]); }
  bscompact(s);
}
static void bgsmul(const BSc &A, const BSc &B, BSc &C) {
  C.clear();
  for (size_t i = 0; i < A.m.size(); i++) for (size_t j = 0; j < B.m.size(); j++) {
    u64 a = A.m[i], b = B.m[j]; if (a & b) continue; bool ng = fsign(a, b);
    size_t o = C.v.size(); C.v.resize(o + LN);
    for (int l = 0; l < LN; l++) { u64 v = mulm(A.v[i * LN + l], B.v[j * LN + l]); C.v[o + l] = ng ? negm(v) : v; }
    C.m.push_back(a | b);
  }
  bscompact(C);
}
static void b_eval_traces(const std::vector<std::vector<int>> &tr, const std::vector<BGM> &pos, u64 top, u64 *out) {
  for (int l = 0; l < LN; l++) out[l] = 0;
  BSc tot; tot.m.push_back(0); tot.v.assign(LN, 1);
  BSc s, ns; BGM cur, nx;
  for (auto &t : tr) {
    if (t.size() == 1) bgtr(pos[t[0]], s);
    else {
      cur = pos[t[0]];
      for (size_t i = 1; i + 1 < t.size(); i++) { bgmul(cur, pos[t[i]], nx); std::swap(cur, nx); if (!cur.size()) return; }
      bgtrmul(cur, pos[t.back()], s);
    }
    if (!s.size()) return;
    bgsmul(tot, s, ns); std::swap(tot, ns); if (!tot.size()) return;
  }
  for (size_t i = 0; i < tot.m.size(); i++) if (tot.m[i] == top) for (int l = 0; l < LN; l++) out[l] = tot.v[i * LN + l];
}
// lanes -> BGM (union of masks, zero fill)
static void b_from_lanes(const std::vector<GMat> &lanes, BGM &out) {
  out.clear(); std::vector<u64> ms; for (auto &g : lanes) for (auto &e : g) ms.push_back(e.m);
  std::sort(ms.begin(), ms.end()); ms.erase(std::unique(ms.begin(), ms.end()), ms.end());
  out.m = ms; out.d.assign(ms.size() * (size_t)NN * LN, 0);
  for (size_t l = 0; l < lanes.size() && l < (size_t)LN; l++) for (auto &e : lanes[l]) {
    size_t i = std::lower_bound(ms.begin(), ms.end(), e.m) - ms.begin(); u32 *b = bblk(out, i);
    for (int x = 0; x < NN; x++) b[x * LN + l] = (u32)addm(b[x * LN + l], e.a[x]);
  }
}
typedef std::function<void(int, BGM &)> BCtx;
static void b_position_mat(const Sector &sec, int type, const Gen &G, int p, const BCtx &ctx, BGM &out) {
  out.clear(); std::vector<Lab> lab; position_labels(sec, type, G, p, lab);
  BGM lm; size_t Bs = (size_t)NN * LN;
  for (auto &a : lab) {
    ctx(a.l, lm);
    for (size_t i = 0; i < lm.m.size(); i++) {
      out.m.push_back(a.m | lm.m[i]); size_t o = out.d.size(); out.d.insert(out.d.end(), bblk(lm, i), bblk(lm, i) + Bs);
      if (a.sg) for (size_t q = 0; q < Bs; q++) out.d[o + q] = (u32)negm(out.d[o + q]);
    }
  }
  bgcompact(out);
}
static void b_eval_gen(const Sector &sec, const Skel &S, const Gen &G, const BCtx &ctx, u64 top_te, u64 *out) {
  NEVAL += LN;
  std::vector<BGM> pos(S.Y);
  for (int p = 0; p < S.Y; p++) { b_position_mat(sec, S.g[S.pg[p]].type, G, p, ctx, pos[p]); if (!pos[p].size()) { for (int l = 0; l < LN; l++) out[l] = 0; return; } }
  b_eval_traces(G.tr, pos, top_te | sec.TOP_EZ, out);
}
// letters of skeleton S at probes ks[0..LN) (ks[l] < 0: zero lane)
static void b_letter_at(const Skel &S, const Probes &pr, const int *ks, int l, BGM &out) {
  out.clear(); if (l < 0) return;
  int gi = S.lgrp[l]; if (gi < 0) return;
  int m = LT[l].odd ? S.g[gi].k : 1;
  out.m.resize(m); out.d.assign((size_t)m * NN * LN, 0);
  for (int r = 0; r < m; r++) {
    out.m[r] = LT[l].odd ? (1ULL << (S.thbase[gi] + r)) : 0;
    u32 *b = bblk(out, r);
    for (int ln = 0; ln < LN; ln++) { if (ks[ln] < 0) continue; const Mat &M = pr[ks[ln]][l][r]; for (int x = 0; x < NN; x++) b[x * LN + ln] = M[x]; }
  }
}
static void b_eval_at(const Sector &sec, const Skel &S, const Gen &G, const Probes &pr, const int *ks, u64 *out) {
  b_eval_gen(sec, S, G, [&](int l, BGM &o) { b_letter_at(S, pr, ks, l, o); }, S.top_te, out);
}
