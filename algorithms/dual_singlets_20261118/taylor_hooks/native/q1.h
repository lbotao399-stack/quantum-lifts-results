// q1.h -- one-loop operator Q1 on trace-word generators and the d1 driver.
//
// Pair table, jet kernel and word rule: su3_engine/q1/x7/q1ah.cpp (= classes_q1_su3 = qn4 FAST_1 with CENG[1] = 1,
// checked there: FAST_1 == q1_poly), colour-blind.  Colour factor K^{AB}_{DE} = kappa^{CC'} C^A_{CD} C^B_{C'E},
// kappa = tr_F, i.e. in matrix form   (Q1 g) = sum_pairs c * tau( .., [t_R, X_Lp] at p, .., [t^R, X_Lq] at q, .. ).
// Evaluators:
//   brute : explicit sum over R < dim g (validation);
//   fast  : sum_R (t_R)_ij (t^R)_kl = c1 d_il d_jk + c2 F_ik F_jl + c3 d_ij d_kl contracted on the trace words;
//   diag  : Cartan probes, only pairs inside one trace, weights CASD_ij = sum_R (t_R)_ij (t^R)_ji.
#pragma once
#include "levels.h"

static int eps3i(int i, int j, int k) { if (i == j || j == k || i == k) return 0; return ((i == 0 && j == 1 && k == 2) || (i == 1 && j == 2 && k == 0) || (i == 2 && j == 0 && k == 1)) ? 1 : -1; }
struct ChanT { int co, tU, fU, tV, fV; };
static inline bool chan_has(int tX, int tY) { return (tX == BF && (tY == CG || tY == BF || tY == GAM)) || (tX == BET && (tY == GAM || tY == BET || tY == BF)); }
static int MUTATE = 0;
static std::vector<ChanT> channel(int tX, int fX, int tY, int fY) {
  std::vector<ChanT> o;
  if (tX == BF && tY == CG) o.push_back({1, CG, 0, CG, 0});
  else if (tX == BF && tY == BF) { o.push_back({1, CG, 0, BF, 0}); o.push_back({-1, BF, 0, CG, 0}); for (int I = 0; I < 3; I++) { o.push_back({1, BET, I, GAM, I}); o.push_back({-1, GAM, I, BET, I}); } }
  else if (tX == BET && tY == GAM) { if (fX == fY) o.push_back({1, CG, 0, CG, 0}); }
  else if (tX == BET && tY == BET) { for (int K = 0; K < 3; K++) { int e = eps3i(fX, fY, K); if (e) { o.push_back({e, CG, 0, GAM, K}); o.push_back({-e, GAM, K, CG, 0}); } } }
  else if (tX == BF && tY == GAM) { o.push_back({1, CG, 0, GAM, fY}); o.push_back({-1, GAM, fY, CG, 0}); }
  else if (tX == BET && tY == BF) {
    int I = fX; o.push_back({1, BET, I, CG, 0}); o.push_back({1, CG, 0, BET, I});
    if (MUTATE != 2) for (int J = 0; J < 3; J++) for (int K = 0; K < 3; K++) { int e = eps3i(I, J, K); if (e) o.push_back({e, GAM, J, GAM, K}); }
  }
  return o;
}
static const int BKS[4] = {1, -1, 1, 1};
struct KerT { int p1, p2, q1, q2; u64 v; };
static std::map<std::array<int, 4>, std::vector<KerT>> KCACHE;
static u64 INVT[256];
static u64 cbm(int n, int k) { if (k < 0 || k > n) return 0; u64 r = 1; for (int i = 1; i <= k; i++) r = mulm(mulm(r, (u64)(n - k + i) % P), INVT[i]); return r; }
static const std::vector<KerT> &kernelB(int a1, int a2, int b1, int b2) {
  std::array<int, 4> key = {a1, a2, b1, b2};
  auto it = KCACHE.find(key); if (it != KCACHE.end()) return it->second;
  std::map<std::array<int, 4>, u64> acc;
  for (int j1 = 0; j1 <= a1; j1++) for (int j2 = 0; j2 <= a2; j2++) {
    u64 pref = mulm(cbm(a1, j1), cbm(a2, j2)); if ((j1 + j2) & 1) pref = (P - pref) % P;
    int c1 = b1 + j1, c2 = b2 + j2, e1 = a1 - j1, e2 = a2 - j2;
    for (int k1 = 0; k1 <= c1; k1++) for (int k2 = 0; k2 <= c2; k2++) {
      u64 g = mulm(mulm(cbm(c1, k1), cbm(c2, k2)), mulm(INVT[c1 + c2 + 2], INVT[k1 + k2 + 1]));
      for (int x1 = 0; x1 <= e1; x1++) for (int x2 = 0; x2 <= e2; x2++) {
        u64 base = mulm(mulm(pref, g), mulm(cbm(e1, x1), cbm(e2, x2)));
        int U1 = k1 + x1, U2 = k2 + x2, V1 = c1 - k1 + e1 - x1, V2 = c2 - k2 + e2 - x2;
        u64 &s1 = acc[{U1 + 1, U2, V1, V2 + 1}]; s1 = addm(s1, base);
        u64 &s2 = acc[{U1, U2 + 1, V1 + 1, V2}]; s2 = subm(s2, base);
      }
    }
  }
  std::vector<KerT> out;
  for (auto &kv : acc) if (kv.second) out.push_back({kv.first[0], kv.first[1], kv.first[2], kv.first[3], kv.second});
  return KCACHE.emplace(key, std::move(out)).first->second;
}
struct Q1T { int Lp, Lq; u64 c; };
static std::mutex Q1MU;
static std::map<std::pair<int, int>, std::vector<Q1T>> Q1CACHE;
static std::vector<Q1T> q1_forward(int X, int Y) {
  std::vector<Q1T> res; const CL &x = LT[X], &y = LT[Y];
  auto fl = channel(x.type, x.flav, y.type, y.flav); if (fl.empty()) return res;
  const std::vector<KerT> &ker = kernelB(x.s, x.r, y.s, y.r);
  std::map<std::pair<int, int>, u64> acc;
  int sXY = BKS[x.type] * BKS[y.type];
  for (auto &ch : fl) {
    u64 s = md(sXY * BKS[ch.tU] * BKS[ch.tV] * ch.co);
    if (MUTATE == 4 && ch.tU == CG && ch.tV == CG) s = mulm(s, 2);
    for (auto &k : ker) {
      int U = ltid(ch.tU, ch.fU, k.p2, k.p1), V = ltid(ch.tV, ch.fV, k.q2, k.q1);
      if (U < 0 || V < 0) { fprintf(stderr, "Q1 output letter outside the letter table\n"); exit(11); }
      u64 &a = acc[{U, V}]; a = addm(a, mulm(s, k.v));
    }
  }
  for (auto &kv : acc) if (kv.second) res.push_back({kv.first.first, kv.first.second, kv.second});
  return res;
}
static const std::vector<Q1T> &q1_pair(int X, int Y) {
  std::lock_guard<std::mutex> lk(Q1MU);
  auto key = std::make_pair(X, Y); auto it = Q1CACHE.find(key); if (it != Q1CACHE.end()) return it->second;
  std::vector<Q1T> res; const CL &x = LT[X], &y = LT[Y];
  if (chan_has(x.type, y.type)) res = q1_forward(X, Y);
  else if (chan_has(y.type, x.type)) {
    for (auto &t : q1_forward(Y, X)) { int sg = (x.odd & y.odd) ^ (LT[t.Lp].odd & LT[t.Lq].odd); res.push_back({t.Lq, t.Lp, sg ? negm(t.c) : t.c}); }
  }
  return Q1CACHE.emplace(key, std::move(res)).first->second;
}
static inline int gkey(int type, int N) { return type * 64 + N; }
static bool pair_target(const Skel &S, int gp, int gq, const Skel &T, std::pair<int, int> &need) {
  std::map<int, int> m;
  for (auto &g : T.g) m[gkey(g.type, g.N)] += g.k;
  for (size_t i = 0; i < S.g.size(); i++) m[gkey(S.g[i].type, S.g[i].N)] -= S.g[i].k - ((int)i == gp) - ((int)i == gq);
  std::vector<int> ab;
  for (auto &kv : m) { if (kv.second < 0) return false; for (int q = 0; q < kv.second; q++) ab.push_back(kv.first); }
  if (ab.size() != 2) return false;
  std::sort(ab.begin(), ab.end()); need = {ab[0], ab[1]}; return true;
}
static inline int lkey(int l) { return gkey(LT[l].type, LT[l].r + LT[l].s); }
static bool skel_connected(const Skel &S, const Skel &T) {
  std::pair<int, int> need;
  for (size_t a = 0; a < S.g.size(); a++) for (size_t b = a; b < S.g.size(); b++) {
    if (a == b && S.g[a].k < 2) continue;
    if (!chan_has(S.g[a].type, S.g[b].type) && !chan_has(S.g[b].type, S.g[a].type)) continue;
    if (pair_target(S, (int)a, (int)b, T, need)) return true;
  }
  return false;
}
static void init_q1() {
  INVT[0] = 0; for (int i = 1; i < 256; i++) INVT[i] = invm((u64)i);
  MUTATE = getenv("Q1_MUTATE") ? atoi(getenv("Q1_MUTATE")) : 0;
  KCACHE.clear(); Q1CACHE.clear();
  if (getenv("Q1_CMUT") && atoi(getenv("Q1_CMUT")) == 1) CC2 = negm(CC2);   // negative control: wrong invariant form term
  if (getenv("Q1_CMUT") && atoi(getenv("Q1_CMUT")) == 2) CC1 = mulm(CC1, 2);   // negative control: c1 doubled
}

// ---------------- pair terms of one generator ----------------
// For each ordered slot pair (p before q in evaluation order) with a channel into T: groups (lb, Lq) -> X-side GMat
struct PairGroup { GMat X; int Lq; u64 mq; };   // X: sum_{la,t} cf (la.m | X_Lp entries); Y-side: lb.m | X_Lq entries
struct PairTerms { int p, q; std::vector<PairGroup> grp; };
typedef std::function<void(int, GMat &)> LetterF;   // letter -> GMat at the probe
static void pair_terms(const Sector &sec, const Skel &S, const Gen &G, const Skel &T, const LetterF &lf, bool same_cycle_only, std::vector<PairTerms> &out) {
  out.clear();
  int Y = S.Y; std::vector<int> ord, cyc(Y); for (size_t c = 0; c < G.tr.size(); c++) for (int p : G.tr[c]) { ord.push_back(p); cyc[p] = (int)c; }
  std::vector<int> od(Y); for (int p = 0; p < Y; p++) od[p] = (S.g[S.pg[p]].type == BET || S.g[S.pg[p]].type == CG);
  std::vector<std::vector<Lab>> lab(Y); for (int p = 0; p < Y; p++) position_labels(sec, S.g[S.pg[p]].type, G, p, lab[p]);
  std::map<int, GMat> lcache; auto letter = [&](int l) -> const GMat & { auto it = lcache.find(l); if (it == lcache.end()) { GMat g; lf(l, g); it = lcache.emplace(l, g).first; } return it->second; };
  for (int a = 0; a < Y; a++) for (int b = a + 1; b < Y; b++) {
    int p = ord[a], q = ord[b];
    if (same_cycle_only && cyc[p] != cyc[q]) continue;
    int tp = S.g[S.pg[p]].type, tq = S.g[S.pg[q]].type;
    if (!chan_has(tp, tq) && !chan_has(tq, tp)) continue;
    std::pair<int, int> need; if (!pair_target(S, S.pg[p], S.pg[q], T, need)) continue;
    int Sp = 0; for (int i = 0; i < a; i++) Sp += od[ord[i]];
    int Sq = 0; for (int i = 0; i < b; i++) if (i != a) Sq += od[ord[i]];
    std::map<std::pair<size_t, int>, GMat> G2;
    for (auto &la : lab[p]) for (size_t ib = 0; ib < lab[q].size(); ib++) {
      const Lab &lb = lab[q][ib];
      if (la.m & lb.m) continue;
      const std::vector<Q1T> &tl = q1_pair(la.l, lb.l);
      for (auto &t : tl) {
        std::pair<int, int> got = {lkey(t.Lp), lkey(t.Lq)}; if (got.first > got.second) std::swap(got.first, got.second);
        if (got != need) continue;
        int E = ((od[p] + LT[t.Lp].odd) * Sp + (od[q] + LT[t.Lq].odd) * Sq + la.sg + lb.sg) & 1;
        u64 cf = E ? negm(t.c) : t.c;
        const GMat &XP = letter(t.Lp); if (XP.empty()) continue;
        if (letter(t.Lq).empty()) continue;
        GMat &acc = G2[{ib, t.Lq}];
        for (auto &x : XP) { GE e; e.m = la.m | x.m; e.a = x.a; mscale(e.a, cf); acc.push_back(e); }
      }
    }
    PairTerms pt; pt.p = p; pt.q = q;
    for (auto &kv : G2) { gcompact(kv.second); if (kv.second.empty()) continue; pt.grp.push_back({kv.second, kv.first.second, lab[q][kv.first.first].m}); }
    if (!pt.grp.empty()) out.push_back(std::move(pt));
  }
}
static void yside(const PairGroup &g, const LetterF &lf, GMat &Yq) { GMat XQ; lf(g.Lq, XQ); Yq.clear(); for (auto &x : XQ) { GE e; e.m = g.mq | x.m; e.a = x.a; Yq.push_back(e); } gcompact(Yq); }

// ---------------- brute-force evaluator (sum over the algebra basis) ----------------
static u64 q1_eval_brute(const Sector &sec, const Skel &S, const Gen &G, const Skel &T, const LetterF &lf, u64 top_te, bool same_cycle_only) {
  int Y = S.Y; std::vector<GMat> base(Y);
  Ctx ctx = [&](int l, GMat &o) { lf(l, o); };
  for (int p = 0; p < Y; p++) position_mat(sec, S.g[S.pg[p]].type, G, p, ctx, base[p]);
  std::vector<PairTerms> pts; pair_terms(sec, S, G, T, lf, same_cycle_only, pts);
  u64 tot = 0, top = top_te | sec.TOP_EZ;
  for (auto &pt : pts) for (auto &g : pt.grp) {
    GMat Yq; yside(g, lf, Yq);
    for (int R = 0; R < DIMG; R++) {
      GMat A, B, t1, t2;
      gmul_const(g.X, TB[R], true, t1); gmul_const(g.X, TB[R], false, t2); A = t1; gadd(A, t2, true);    // [t_R, X]
      gmul_const(Yq, TD[R], true, t1); gmul_const(Yq, TD[R], false, t2); B = t1; gadd(B, t2, true);      // [t^R, Y]
      std::vector<GMat> pos = base; pos[pt.p] = A; pos[pt.q] = B;
      tot = addm(tot, eval_traces(G.tr, pos, top));
    }
  }
  return tot;
}

// ---------------- fast evaluator (closed-form Casimir contraction) ----------------
static GMat gid() { GMat I(1); I[0].m = 0; mzero(I[0].a); for (int i = 0; i < NF; i++) I[0].a[i * NF + i] = 1; return I; }
static void gprod(const std::vector<GMat> &pos, const std::vector<int> &cy, int from, int to, GMat &out) {   // product pos[cy[from..to-1]]
  if (from >= to) { out = gid(); return; }
  out = pos[cy[from]]; GMat nx; for (int i = from + 1; i < to; i++) { gmul(out, pos[cy[i]], nx); std::swap(out, nx); }
}
static void gscale(GMat &A, u64 c) { for (auto &e : A) mscale(e.a, c); }
// Psi(U,V,W) = sum_R tr(U t_R V t^R W)
static void psi3(const GMat &U, const GMat &V, const GMat &W, GSc &out) {
  out.clear(); GSc s; GMat t1, t2;
  if (CC1) { GSc tv; gtr(V, tv); sgmul(tv, W, t1); gtrmul(U, t1, s); gsadd(out, s, CC1); }
  if (CC2) { GMat VT, A, B; gmul_const(U, FM, false, A); gtransp(V, VT); gmul_const(VT, FM, false, B); gmul(A, B, t1); gtrmul(t1, W, s); gsadd(out, s, CC2); }
  if (CC3) { gmul(U, V, t1); gtrmul(t1, W, s); gsadd(out, s, CC3); }
}
// Phi(Q1, S, Q2) = sum_{ijkl} C_{ij,kl} Q1_ij S Q2_kl
static void phi2(const GMat &Q1, const GSc &Smid, const GMat &Q2, GSc &out) {
  out.clear(); GSc s; GMat t1, t2;
  if (CC1) { sgmul(Smid, Q2, t1); gtrmul(Q1, t1, s); gsadd(out, s, CC1); }
  if (CC2) { GMat A, B, BT; gmul_const(Q2, FM, true, A); Mat FT; mtranspose(FM, FT); gmul_const(A, FT, false, B); gtransp(B, BT); sgmul(Smid, BT, t1); gtrmul(Q1, t1, s); gsadd(out, s, CC2); }
  if (CC3) { GSc a, b, ab; gtr(Q1, a); gtr(Q2, b); gsmul(a, Smid, ab); gsmul(ab, b, s); gsadd(out, s, CC3); }
}
static u64 q1_eval_fast(const Sector &sec, const Skel &S, const Gen &G, const Skel &T, const LetterF &lf, u64 top_te) {
  int Y = S.Y; std::vector<GMat> base(Y);
  Ctx ctx = [&](int l, GMat &o) { lf(l, o); };
  for (int p = 0; p < Y; p++) position_mat(sec, S.g[S.pg[p]].type, G, p, ctx, base[p]);
  std::vector<PairTerms> pts; pair_terms(sec, S, G, T, lf, false, pts);
  if (pts.empty()) return 0;
  size_t nc = G.tr.size(); std::vector<int> cyc(Y), cpos(Y);
  for (size_t c = 0; c < nc; c++) for (size_t i = 0; i < G.tr[c].size(); i++) { cyc[G.tr[c][i]] = (int)c; cpos[G.tr[c][i]] = (int)i; }
  std::vector<GSc> trc(nc);
  for (size_t c = 0; c < nc; c++) { GMat W; gprod(base, G.tr[c], 0, (int)G.tr[c].size() - 1, W); gtrmul(W, base[G.tr[c].back()], trc[c]); }
  auto trprod = [&](size_t from, size_t to) { GSc a = {{0, 1}}, nx; for (size_t c = from; c < to; c++) { gsmul(a, trc[c], nx); std::swap(a, nx); } return a; };
  u64 top = top_te | sec.TOP_EZ, tot = 0;
  for (auto &pt : pts) {
    int cp = cyc[pt.p], cq = cyc[pt.q];
    const std::vector<int> &C1 = G.tr[cp], &C2 = G.tr[cq];
    GSc pre = trprod(0, cp), post = trprod(cq + 1, nc), val, acc, t1, t2;
    if (cp == cq) {
      int x = cpos[pt.p], y = cpos[pt.q];
      GMat Pre, Mid, Suf; gprod(base, C1, 0, x, Pre); gprod(base, C1, x + 1, y, Mid); gprod(base, C1, y + 1, (int)C1.size(), Suf);
      for (auto &g : pt.grp) {
        GMat Yq; yside(g, lf, Yq);
        GMat XM, PX, MY, YS, XMY; gmul(g.X, Mid, XM); gmul(Pre, g.X, PX); gmul(Mid, Yq, MY); gmul(Yq, Suf, YS); gmul(XM, Yq, XMY);
        psi3(Pre, XM, YS, t1); gsadd(acc, t1, 1);
        psi3(Pre, XMY, Suf, t1); gsadd(acc, t1, P - 1);
        psi3(PX, Mid, YS, t1); gsadd(acc, t1, P - 1);
        psi3(PX, MY, Suf, t1); gsadd(acc, t1, 1);
      }
    } else {
      int x = cpos[pt.p], y = cpos[pt.q];
      GMat Pre1, Suf1, Pre2, Suf2; gprod(base, C1, 0, x, Pre1); gprod(base, C1, x + 1, (int)C1.size(), Suf1); gprod(base, C2, 0, y, Pre2); gprod(base, C2, y + 1, (int)C2.size(), Suf2);
      GSc mid = trprod(cp + 1, cq);
      for (auto &g : pt.grp) {
        GMat Yq; yside(g, lf, Yq);
        // tr(U1 t_R V1): (U1,V1) = (Pre1, X Suf1) with +, (Pre1 X, Suf1) with -;  Q = U^T . V^T
        GMat XS, PX, YS2, PY, a, b, Qa[2], Qb[2];
        gmul(g.X, Suf1, XS); gmul(Pre1, g.X, PX); gmul(Yq, Suf2, YS2); gmul(Pre2, Yq, PY);
        gtransp(Pre1, a); gtransp(XS, b); gmul(a, b, Qa[0]);
        gtransp(PX, a); gtransp(Suf1, b); gmul(a, b, Qa[1]);
        gtransp(Pre2, a); gtransp(YS2, b); gmul(a, b, Qb[0]);
        gtransp(PY, a); gtransp(Suf2, b); gmul(a, b, Qb[1]);
        for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) { phi2(Qa[i], mid, Qb[j], t1); gsadd(acc, t1, ((i + j) & 1) ? P - 1 : 1); }
      }
    }
    if (acc.empty()) continue;
    gsmul(pre, acc, t1); gsmul(t1, post, val);
    for (auto &x : val) if (x.first == top) tot = addm(tot, x.second);
  }
  return tot;
}

// ---------------- diagonal (Cartan) evaluator ----------------
static void gs_mono(const GSc &A, u64 m, u64 c, GSc &out) {
  out.clear(); for (auto &a : A) { if (a.first & m) continue; u64 v = mulm(a.second, c); if (fsign(a.first, m)) v = negm(v); out.push_back({a.first | m, v}); }
  std::sort(out.begin(), out.end());
}
static u64 q1_eval_diag(const Sector &sec, const Skel &S, const Gen &G, const Skel &T, const LetterF &lf, u64 top_te) {
  int Y = S.Y; std::vector<int> ord, cyc(Y), cpos(Y); for (size_t c = 0; c < G.tr.size(); c++) for (size_t i = 0; i < G.tr[c].size(); i++) { int p = G.tr[c][i]; ord.push_back(p); cyc[p] = (int)c; cpos[p] = (int)i; }
  std::vector<int> rk(Y); for (int i = 0; i < Y; i++) rk[ord[i]] = i;
  std::vector<int> od(Y); for (int p = 0; p < Y; p++) od[p] = (S.g[S.pg[p]].type == BET || S.g[S.pg[p]].type == CG);
  std::vector<std::vector<Lab>> lab(Y); for (int p = 0; p < Y; p++) position_labels(sec, S.g[S.pg[p]].type, G, p, lab[p]);
  std::vector<std::pair<int, int>> pairs, needs;
  for (size_t c = 0; c < G.tr.size(); c++) for (size_t a = 0; a < G.tr[c].size(); a++) for (size_t b = a + 1; b < G.tr[c].size(); b++) {
    int p = G.tr[c][a], q = G.tr[c][b]; std::pair<int, int> need;
    if (!chan_has(S.g[S.pg[p]].type, S.g[S.pg[q]].type) && !chan_has(S.g[S.pg[q]].type, S.g[S.pg[p]].type)) continue;
    if (!pair_target(S, S.pg[p], S.pg[q], T, need)) continue;
    pairs.push_back({p, q}); needs.push_back(need);
  }
  if (pairs.empty()) return 0;
  u64 top = top_te | sec.TOP_EZ;
  std::vector<std::vector<GSc>> Xd(Y, std::vector<GSc>(NF));
  GMat lm;
  for (int p = 0; p < Y; p++) {
    for (auto &a : lab[p]) { lf(a.l, lm); for (auto &e : lm) for (int i = 0; i < NF; i++) { u64 v = e.a[i * NF + i]; if (a.sg) v = negm(v); if (v) Xd[p][i].push_back({a.m | e.m, v}); } }
    for (int i = 0; i < NF; i++) scompact(Xd[p][i]);
  }
  size_t nc = G.tr.size(); std::vector<GSc> trc(nc);
  for (size_t c = 0; c < nc; c++) {
    GSc acc;
    for (int i = 0; i < NF; i++) { GSc cur = {{0, 1}}, nx; for (int p : G.tr[c]) { gsmul(cur, Xd[p][i], nx); std::swap(cur, nx); if (cur.empty()) break; } for (auto &x : cur) acc.push_back(x); }
    scompact(acc); trc[c] = acc;
  }
  std::vector<GSc> pre(nc + 1), suf(nc + 1); pre[0] = {{0, 1}}; suf[nc] = {{0, 1}};
  for (size_t c = 0; c < nc; c++) gsmul(pre[c], trc[c], pre[c + 1]);
  for (size_t c = nc; c-- > 0;) gsmul(trc[c], suf[c + 1], suf[c]);
  u64 tot = 0;
  for (size_t pi = 0; pi < pairs.size(); pi++) {
    int p = pairs[pi].first, q = pairs[pi].second; int c = cyc[p]; const std::vector<int> &cy = G.tr[c];
    int x = cpos[p], y = cpos[q];
    const GSc &Pre = pre[c], &Post = suf[c + 1];
    int Sp = 0; for (int i = 0; i < rk[p]; i++) Sp += od[ord[i]];
    int Sq = 0; for (int i = 0; i < rk[q]; i++) if (i != rk[p]) Sq += od[ord[i]];
    std::map<std::pair<u64, u64>, std::vector<u64>> W;
    GMat XP, XQ;
    for (auto &la : lab[p]) for (auto &lb : lab[q]) {
      if (la.m & lb.m) continue;
      const std::vector<Q1T> &tl = q1_pair(la.l, lb.l);
      for (auto &t : tl) {
        std::pair<int, int> got = {lkey(t.Lp), lkey(t.Lq)}; if (got.first > got.second) std::swap(got.first, got.second);
        if (got != needs[pi]) continue;
        int E = ((od[p] + LT[t.Lp].odd) * Sp + (od[q] + LT[t.Lq].odd) * Sq + la.sg + lb.sg) & 1;
        u64 cf = E ? negm(t.c) : t.c;
        lf(t.Lp, XP); lf(t.Lq, XQ);
        for (auto &xa : XP) for (auto &xb : XQ) {
          if (xa.m & xb.m) continue;
          u64 ma = la.m | xa.m, mb = lb.m | xb.m;
          auto &w = W[{ma, mb}]; if (w.empty()) w.assign(NN, 0);
          for (int i = 0; i < NF; i++) for (int j = 0; j < NF; j++) {
            u64 cw = CASD[i * NF + j]; if (i == j || !cw) continue;
            u64 A = subm(xa.a[j * NF + j], xa.a[i * NF + i]), B = subm(xb.a[i * NF + i], xb.a[j * NF + j]);
            w[i * NF + j] = addm(w[i * NF + j], mulm(cf, mulm(cw, mulm(A, B))));
          }
        }
      }
    }
    if (W.empty()) continue;
    for (int i = 0; i < NF; i++) for (int j = 0; j < NF; j++) {
      if (i == j || !CASD[i * NF + j]) continue;
      GSc Pi = {{0, 1}}, Mj = {{0, 1}}, Si = {{0, 1}}, nx;
      for (int t = 0; t < x; t++) { gsmul(Pi, Xd[cy[t]][i], nx); std::swap(Pi, nx); }
      for (int t = x + 1; t < y; t++) { gsmul(Mj, Xd[cy[t]][j], nx); std::swap(Mj, nx); }
      for (int t = y + 1; t < (int)cy.size(); t++) { gsmul(Si, Xd[cy[t]][i], nx); std::swap(Si, nx); }
      GSc L, R; gsmul(Pre, Pi, L); gsmul(Si, Post, R);
      if (L.empty() || Mj.empty() || R.empty()) continue;
      std::map<u64, std::unordered_map<u64, u64>> Z;
      for (auto &kv : W) {
        u64 w = kv.second[i * NF + j]; if (!w) continue;
        u64 ma = kv.first.first, mb = kv.first.second;
        auto &z = Z[ma];
        for (auto &r : R) { if (r.first & mb) continue; u64 v = mulm(w, r.second); if (fsign(mb, r.first)) v = negm(v); u64 &s = z[mb | r.first]; s = addm(s, v); }
      }
      GSc La, LaM;
      for (auto &zz : Z) {
        gs_mono(L, zz.first, 1, La); if (La.empty()) continue;
        gsmul(La, Mj, LaM);
        for (auto &m : LaM) {
          u64 rest = top & ~m.first; if ((m.first | rest) != top) continue;
          auto f = zz.second.find(rest); if (f == zz.second.end() || !f->second) continue;
          u64 v = mulm(m.second, f->second); if (fsign(m.first, rest)) v = negm(v);
          tot = addm(tot, v);
        }
      }
    }
  }
  return tot;
}

static int q1_selftest() { return 0; }
#include "d1.h"
