// d1.h -- drivers: d1 rank (Cartan shortcut / general exact-in-span), Q1 evaluator cross-checks, Q0Q1+Q1Q0 = 0 gate.
#pragma once

// coefficients of v (values at the pivot probes of S) in the basis S.gens: solve EVP c = v
static std::vector<u64> skel_coeffs(const Skel &S, const std::vector<u64> &vp) {
  size_t n = S.gens.size(); std::vector<u64> c(n, 0);
  for (size_t j = 0; j < n; j++) { u64 s = 0; for (size_t q = 0; q < n; q++) s = red(s + (u64)S.Einv[j][q] * vp[q]); c[j] = s; }
  return c;
}
static void make_einv(Skel &S) {
  size_t n = S.gens.size(); std::vector<std::vector<u64>> A(n, std::vector<u64>(n));
  for (size_t q = 0; q < n; q++) for (size_t j = 0; j < n; j++) A[q][j] = S.raw[j][S.piv[q]];
  std::vector<std::vector<u64>> I = inv_matrix(A);
  S.Einv.assign(n, std::vector<u32>(n)); for (size_t j = 0; j < n; j++) for (size_t q = 0; q < n; q++) S.Einv[j][q] = (u32)I[j][q];
}
static void ratrec(u64 a, long &nu, long &de) {   // a = nu/de mod P with |nu|, de <= sqrt(P/2)
  i64 r0 = (i64)P, r1 = (i64)a, s0 = 0, s1 = 1; i64 B = (i64)std::sqrt((double)P / 2);
  while (r1 > B) { i64 q = r0 / r1, t = r0 - q * r1; r0 = r1; r1 = t; t = s0 - q * s1; s0 = s1; s1 = t; }
  nu = (long)r1; de = (long)s1; if (de < 0) { nu = -nu; de = -de; }
}
static LetterF probe_letters(const Skel &T, const Probes &pr, int k) { return [&T, &pr, k](int l, GMat &o) { letter_at(T, pr, k, l, o); }; }

// Phi[row][col]: rows = (target skeleton, probe) list, cols = source basis (global); evaluator kind 0 fast, 1 diag, 2 brute
static Mtx q1_matrix(const Sector &ss, const Lev &A, const Sector &ts, const Lev &B, const std::vector<std::pair<int, int>> &rows, bool diag, int kind) {
  size_t nc; std::vector<size_t> c0 = col_offsets(A, nc);
  Mtx Phi(rows.size(), std::vector<u32>(nc, 0));
  struct Task { size_t r, si; };
  std::vector<Task> tasks;
  for (size_t r = 0; r < rows.size(); r++) for (size_t si = 0; si < A.sk.size(); si++) if (skel_connected(A.sk[si], B.sk[rows[r].first])) tasks.push_back({r, si});
  std::atomic<size_t> done(0); double t0 = now();
  par_for(tasks.size(), [&](size_t t) {
    size_t r = tasks[t].r, si = tasks[t].si; const Skel &T = B.sk[rows[r].first]; int k = rows[r].second;
    const Skel &S = A.sk[si];
    LetterF lf = probe_letters(T, diag ? T.dg : T.pr, k);
    for (size_t i = 0; i < S.gens.size(); i++) {
      u64 v = kind == 1 ? q1_eval_diag(ss, S, S.gens[i], T, lf, T.top_te) : kind == 2 ? q1_eval_brute(ss, S, S.gens[i], T, lf, T.top_te, false) : q1_eval_fast(ss, S, S.gens[i], T, lf, T.top_te);
      Phi[r][c0[si] + i] = (u32)v;
    }
    size_t d = ++done; if (VERB > 1 && d % 500 == 0) fprintf(stderr, "#q1 tasks %zu/%zu %.1fs\n", d, tasks.size(), now() - t0);
  });
  (void)ts;
  return Phi;
}
static std::vector<std::pair<int, int>> pivot_rows(const Lev &B) { std::vector<std::pair<int, int>> r; for (size_t ti = 0; ti < B.sk.size(); ti++) for (int k : B.sk[ti].piv) r.push_back({(int)ti, k}); return r; }
static Mtx concat_cols(const Mtx &A, size_t na, const Mtx &B, size_t nb) {
  size_t n = std::max(A.size(), B.size()); Mtx C(n, std::vector<u32>(na + nb, 0));
  for (size_t i = 0; i < A.size(); i++) for (size_t j = 0; j < na; j++) C[i][j] = A[i][j];
  for (size_t i = 0; i < B.size(); i++) for (size_t j = 0; j < nb; j++) C[i][na + j] = B[i][j];
  return C;
}

template <class F> static int q1_main(const std::string &mode, int argc, char **argv, int *w, const Sector &src, F setseed, double T0) {
  if (argc < 9) { fprintf(stderr, "need Y\n"); return 1; }
  int Y = atoi(argv[8]);
  int wt[5] = {w[0] - 1, w[1] - 1, w[2] - 1, w[3] + 1, w[4] + 1};
  if (wt[2] < 0) { fprintf(stderr, "no target sector (t3 = 0)\n"); return 1; }
  Sector tgt = make_sector(wt);
  int ngen = 0; u64 seed;
  if (mode == "q1check") { if (argc < 10) { fprintf(stderr, "need NGEN\n"); return 1; } ngen = atoi(argv[9]); seed = setseed(10); } else seed = setseed(9);
  init_all({src, tgt});
  print_head(mode.c_str(), w, src, seed);
  printf(",\"target\":[%d,%d,%d,%d,%d],\"Y\":%d", wt[0], wt[1], wt[2], wt[3], wt[4], Y);
  if (mode == "q1check") {
    Lev L; build_lev(src, L, Y); Lev LTg; build_lev(tgt, LTg, Y);
    std::vector<std::pair<int, int>> pick;
    for (int round = 0; (int)pick.size() < ngen && round < 1000; round++) { bool any = false; for (size_t si = 0; si < L.sk.size() && (int)pick.size() < ngen; si++) { bool conn = false; for (auto &T : LTg.sk) if (skel_connected(L.sk[si], T)) conn = true; if (conn && round < (int)L.sk[si].gens.size()) { pick.push_back({(int)si, round}); any = true; } } if (!any) break; }
    size_t tg = 0, bg = 0, nzg = 0, td = 0, bd = 0, nzd = 0; double tb = 0, tf = 0, tdg = 0;
    for (auto &pk : pick) {
      const Skel &S = L.sk[pk.first]; const Gen &G = S.gens[pk.second];
      for (auto &T : LTg.sk) {
        if (!skel_connected(S, T)) continue;
        Probes Pg, Pd; add_probes(T, Pg, 2, false); add_probes(T, Pd, 2, true);
        for (int k = 0; k < 2; k++) {
          LetterF lg = probe_letters(T, Pg, k), ld = probe_letters(T, Pd, k);
          double a = now(); u64 b = q1_eval_brute(src, S, G, T, lg, T.top_te, false); tb += now() - a;
          a = now(); u64 f = q1_eval_fast(src, S, G, T, lg, T.top_te); tf += now() - a;
          tg++; if (b != f) bg++; if (b) nzg++;
          u64 bd1 = q1_eval_brute(src, S, G, T, ld, T.top_te, false), bd2 = q1_eval_brute(src, S, G, T, ld, T.top_te, true), fd = q1_eval_fast(src, S, G, T, ld, T.top_te);
          a = now(); u64 dd = q1_eval_diag(src, S, G, T, ld, T.top_te); tdg += now() - a;
          td++; if (bd1 != bd2 || bd1 != fd || bd1 != dd) bd++; if (bd1) nzd++;
        }
      }
    }
    printf(",\"generators\":%zu,\"general_tests\":%zu,\"general_mismatch\":%zu,\"general_nonzero\":%zu,\"diag_tests\":%zu,\"diag_mismatch\":%zu,\"diag_nonzero\":%zu,\"t_brute\":%.3f,\"t_fast\":%.3f,\"t_diag\":%.3f,\"seconds\":%.1f}\n",
           pick.size(), tg, bg, nzg, td, bd, nzd, tb, tf, tdg, now() - T0);
    return (bg || bd) ? 3 : 0;
  }
  if (mode == "ratio") {
    // (Q1 g)(x) / g_t(x) for every source basis generator g and every 1-dimensional target skeleton, at 3 probes
    Lev L; build_lev(src, L, Y); Lev LTg; build_lev(tgt, LTg, Y);
    printf(",\"ratios\":[");
    bool first = true; size_t incons = 0;
    for (auto &S : L.sk) for (size_t i = 0; i < S.gens.size(); i++) for (auto &T : LTg.sk) {
      if (T.gens.size() != 1 || !skel_connected(S, T)) continue;
      Probes Pr; add_probes(T, Pr, 3, false);
      std::vector<u64> r;
      for (int k = 0; k < 3; k++) { LetterF lf = probe_letters(T, Pr, k); u64 a = q1_eval_fast(src, S, S.gens[i], T, lf, T.top_te), b = eval_at(tgt, T, T.gens[0], Pr, k); r.push_back(b ? mulm(a, invm(b)) : (u64)-1); }
      if (r[0] != r[1] || r[1] != r[2]) incons++;
      long nu, de; ratrec(r[0], nu, de);
      printf("%s{\"src\":\"%s\",\"tgt\":\"%s\",\"ratio\":\"%ld/%ld\"}", first ? "" : ",", S.key.c_str(), T.key.c_str(), nu, de); first = false;
    }
    printf("],\"inconsistent\":%zu}\n", incons);
    return 0;
  }
  if (mode == "anticomm") {
    // (Q0 Q1 + Q1 Q0) g' = 0 for every basis generator g' of V(s, Y-1), in the pivot coordinates of V(t, Y)
    KEEPRAW = true;
    Lev Sm, S0, Tm, T0l; build_lev(src, Sm, Y - 1); build_lev(src, S0, Y); build_lev(tgt, Tm, Y - 1); build_lev(tgt, T0l, Y);
    for (auto &S : S0.sk) make_einv(S); for (auto &S : Tm.sk) make_einv(S);
    std::vector<std::pair<int, int>> rT0 = pivot_rows(T0l), rTm = pivot_rows(Tm);
    Mtx Q0s = q0_rows(src, Sm, S0);              // rows: pivots of S0, cols: Sm basis
    Mtx Q0t = q0_rows(tgt, Tm, T0l);             // rows: pivots of T0, cols: Tm basis
    Mtx Q1a = q1_matrix(src, S0, tgt, T0l, rT0, false, 0);   // rows pivots T0, cols S0 basis
    Mtx Q1b = q1_matrix(src, Sm, tgt, Tm, rTm, false, 0);    // rows pivots Tm, cols Sm basis
    size_t ncm; col_offsets(Sm, ncm); size_t nc0; std::vector<size_t> c0s = col_offsets(S0, nc0); size_t nct; std::vector<size_t> c0t = col_offsets(Tm, nct);
    // coefficient matrices: C0s = coords of Q0 g' in S0 basis (nc0 x ncm), C1b = coords of Q1 g' in Tm basis (nct x ncm)
    Mtx C0s(nc0, std::vector<u32>(ncm, 0)), C1b(nct, std::vector<u32>(ncm, 0));
    { size_t r = 0; for (size_t si = 0; si < S0.sk.size(); si++) { size_t n = S0.sk[si].gens.size(); for (size_t col = 0; col < ncm; col++) { std::vector<u64> vp(n); for (size_t q = 0; q < n; q++) vp[q] = Q0s[r + q][col]; std::vector<u64> c = skel_coeffs(S0.sk[si], vp); for (size_t j = 0; j < n; j++) C0s[c0s[si] + j][col] = (u32)c[j]; } r += n; } }
    { size_t r = 0; for (size_t ti = 0; ti < Tm.sk.size(); ti++) { size_t n = Tm.sk[ti].gens.size(); for (size_t col = 0; col < ncm; col++) { std::vector<u64> vp(n); for (size_t q = 0; q < n; q++) vp[q] = Q1b[r + q][col]; std::vector<u64> c = skel_coeffs(Tm.sk[ti], vp); for (size_t j = 0; j < n; j++) C1b[c0t[ti] + j][col] = (u32)c[j]; } r += n; } }
    Mtx A1 = mul_abt(Q1a, transpose(C0s, ncm), nc0);   // Q1 Q0 g'
    Mtx A2 = mul_abt(Q0t, transpose(C1b, ncm), nct);   // Q0 Q1 g'
    size_t bad = 0, nz1 = 0, nz2 = 0, tot = 0;
    for (size_t r = 0; r < A1.size(); r++) for (size_t c = 0; c < ncm; c++) { tot++; if (A1[r][c]) nz1++; if (A2[r][c]) nz2++; if (addm(A1[r][c], A2[r][c])) bad++; }
    printf(",\"src_dims\":[%zu,%zu],\"tgt_dims\":[%zu,%zu],\"entries\":%zu,\"Q1Q0_nonzero\":%zu,\"Q0Q1_nonzero\":%zu,\"anticomm_violations\":%zu,\"seconds\":%.1f,\"maxrss_MB\":%.0f}\n",
           Sm.dim, S0.dim, Tm.dim, T0l.dim, tot, nz1, nz2, bad, now() - T0, maxrss_mb());
    return bad ? 3 : 0;
  }
  // ---- d1 ----
  std::string dm = getenv("D1_MODE") ? getenv("D1_MODE") : "auto";
  KEEPRAW = getenv("D1_CLOSURE") && atoi(getenv("D1_CLOSURE"));
  Cohom Ct; compute_cohom(Ct, tgt, Y, true, true, true);
  Cohom Cs; compute_cohom(Cs, src, Y, true, false, getenv("SRC_RES") != nullptr);
  printf(",\"src\":{\"dim\":%zu,\"rank_in\":%zu,\"rank_out\":%zu,\"h\":%ld%s},\"tgt\":{\"dim\":%zu,\"rank_in\":%zu,\"rank_out\":%zu,\"h\":%ld,\"rho\":%zu}",
         Cs.dim, Cs.rank_in, Cs.rank_out, Cs.h, getenv("SRC_RES") ? (",\"rho\":" + std::to_string(Cs.rho)).c_str() : "", Ct.dim, Ct.rank_in, Ct.rank_out, Ct.h, Ct.rho);
  bool do_cartan = dm == "cartan" || dm == "both" || (dm == "auto" && Ct.h == (long)Ct.rho);
  bool do_general = dm == "general" || dm == "both" || (dm == "auto" && Ct.h != (long)Ct.rho);
  long rc = -1, rg = -1;
  if (Cs.K.empty() || Ct.h == 0) { printf(",\"d1_rank\":0,\"reason\":\"%s\"", Cs.K.empty() ? "no source cocycles" : "target h = 0"); do_cartan = do_general = false; rc = rg = 0; }
  if (do_cartan) {
    double t1 = now();
    Mtx Phi = q1_matrix(src, Cs.L0, tgt, Ct.L0, Ct.witness, true, 1);
    Mtx D1 = mul_abt(Phi, Cs.K, Cs.nc);
    rc = (long)rank_of(D1, Cs.K.size());
    printf(",\"cartan\":{\"witness\":%zu,\"d1_rank\":%ld,\"exact\":%s,\"t\":%.1f}", Ct.witness.size(), rc, Ct.h == (long)Ct.rho ? "true" : "false", now() - t1);
  }
  if (do_general) {
    double t1 = now();
    std::vector<std::pair<int, int>> rows = pivot_rows(Ct.L0);
    Mtx Phi = q1_matrix(src, Cs.L0, tgt, Ct.L0, rows, false, 0);
    Mtx Q1Z = mul_abt(Phi, Cs.K, Cs.nc);                         // target coords x source cocycles
    size_t nm = Ct.Lm.dim; Mtx Cat = concat_cols(Ct.Min, nm, Q1Z, Cs.K.size());
    size_t rcat = rank_of(Cat, nm + Cs.K.size());
    rg = (long)rcat - (long)Ct.rank_in;
    // closure (Q0 Q1 z = 0) on all source cocycles, optional
    long closure_bad = -1;
    if (KEEPRAW) {
      for (auto &T : Ct.L0.sk) make_einv(T);
      Mtx Mout = q0_rows(tgt, Ct.L0, Ct.Lp);
      size_t ntc; std::vector<size_t> c0t = col_offsets(Ct.L0, ntc);
      Mtx Cz(ntc, std::vector<u32>(Cs.K.size(), 0));
      size_t r = 0; for (size_t ti = 0; ti < Ct.L0.sk.size(); ti++) { size_t n = Ct.L0.sk[ti].gens.size(); for (size_t z = 0; z < Cs.K.size(); z++) { std::vector<u64> vp(n); for (size_t q = 0; q < n; q++) vp[q] = Q1Z[r + q][z]; std::vector<u64> c = skel_coeffs(Ct.L0.sk[ti], vp); for (size_t j = 0; j < n; j++) Cz[c0t[ti] + j][z] = (u32)c[j]; } r += n; }
      Mtx QQ = mul_abt(Mout, transpose(Cz, Cs.K.size()), ntc);
      closure_bad = 0; for (auto &row : QQ) for (auto x : row) if (x) closure_bad++;
    }
    size_t nzq = 0; for (auto &row : Q1Z) for (auto x : row) if (x) nzq++;
    printf(",\"general\":{\"coords\":%zu,\"Q1Z_nonzero\":%zu,\"rank_B\":%zu,\"rank_B_plus_Q1Z\":%zu,\"d1_rank\":%ld,\"closure_violations\":%ld,\"t\":%.1f}", rows.size(), nzq, Ct.rank_in, rcat, rg, closure_bad, now() - t1);
  }
  long r = rg >= 0 ? rg : rc;
  printf(",\"d1_rank\":%ld,\"seconds\":%.1f,\"evals\":%zu,\"maxrss_MB\":%.0f}\n", r, now() - T0, (size_t)NEVAL, maxrss_mb());
  return 0;
}
