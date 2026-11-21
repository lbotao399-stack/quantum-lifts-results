// engine.cpp -- N=4 SYM HW Q0-cohomology and one-loop d1 for G = SU(N), SO(2n+1), Sp(n) (trace-word A-H engine).
//
// Modes (GAUGE in {SU2,SU3,SO3,SO5,SO7,Sp1,Sp2,Sp3}):
//   engine selftest GAUGE [seed] [prime]
//   engine dims  GAUGE t1 t2 t3 n+ n- [Ylo Yhi]                 exact per-Y HW cochain dimensions (+ Euler characteristic)
//   engine chain GAUGE t1 t2 t3 n+ n- Ylo Yhi [seed] [prime]    dim, Q0 ranks, h_Y, Cartan rank rho_Y for Ylo <= Y <= Yhi
//   engine d1    GAUGE t1 t2 t3 n+ n- Y [seed] [prime]          rank of d1 = [Q1]: H(s,Y) -> H(s-(1,1,1,-1,-1),Y)
//   engine q1check GAUGE t1 t2 t3 n+ n- Y NGEN [seed] [prime]   Q1 evaluator cross-checks (general fast vs brute force vs diag)
//   engine anticomm GAUGE t1 t2 t3 n+ n- Y [seed] [prime]       (Q0 Q1 + Q1 Q0) = 0 on V(s,Y-1) -> V(s-d,Y) (pivot coordinates)
// env: A_MARGIN (3), PAR_THREADS (1), A_CHECK=1 (HW check of every generator), D1_MODE=auto|cartan|general|both,
//      Q1_MUTATE=k (negative controls), A_VERB (1)
#include "levels.h"
static void init_q1();
static void init_all(const std::vector<Sector> &secs) {
  for (int i = 0; i < 5; i++) BOX[i] = 0;
  int ym = 0; for (auto &s : secs) { for (int i = 0; i < 5; i++) BOX[i] = std::max(BOX[i], s.WS[i]); ym = std::max(ym, s.Ymax); }
  build_algebra(); build_char_grid(ym + 1); build_letters(); build_q(); build_raise(); init_q1();
}
static void print_head(const char *mode, const int w[5], const Sector &s, u64 seed) {
  printf("{\"mode\":\"%s\",\"gauge\":\"%s\",\"sector\":[%d,%d,%d,%d,%d],\"flavour\":[%d,%d],\"spin2\":%d,\"j\":%d,\"prime\":%llu,\"seed\":%llu",
         mode, GNAME.c_str(), w[0], w[1], w[2], w[3], w[4], w[0] - w[1], w[1] - w[2], w[3] - w[4], s.j, (unsigned long long)P, (unsigned long long)seed);
}

#include "q1.h"

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: engine MODE GAUGE ...\n"); return 1; }
  std::string mode = argv[1];
  if (!parse_gauge(argv[2])) { fprintf(stderr, "bad gauge %s\n", argv[2]); return 1; }
  set_mm(); set_bmm();
  if (getenv("A_MARGIN")) MARGIN = atoi(getenv("A_MARGIN"));
  if (getenv("PAR_THREADS")) NTHREADS = atoi(getenv("PAR_THREADS"));
  if (getenv("A_CHECK")) HWCHECK = atoi(getenv("A_CHECK")) != 0;
  if (getenv("A_VERB")) VERB = atoi(getenv("A_VERB"));
  if (getenv("A_SCALAR")) SCALAR = atoi(getenv("A_SCALAR")) != 0;
  auto setseed = [&](int ai) -> u64 {
    u64 SEED = 1;
    if (argc > ai) SEED = strtoull(argv[ai], 0, 10);
    u64 z = SEED * 0xD1B54A32D192ED03ULL + 0x8CB92BA72F3D8DD7ULL; z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; RS = z ^ (z >> 31);
    if (argc > ai + 1) { P = strtoull(argv[ai + 1], 0, 10); MERS = (P == 2147483647ULL); if (P >= (1ULL << 31)) { fprintf(stderr, "prime must be < 2^31\n"); exit(1); } }
    return SEED;
  };
  double T0 = now();
  if (mode == "selftest") {
    u64 seed = setseed(3);
    int w[5] = {1, 1, 1, 1, 1}; Sector s = make_sector(w); init_all({s});
    printf("{\"mode\":\"selftest\",\"gauge\":\"%s\",\"NF\":%d,\"dimg\":%d,\"rank\":%d,\"C_F\":%llu,\"char_grid\":%d,\"char_prime\":%llu,\"prime\":%llu,\"seed\":%llu,\"q1_selftest\":%d}\n",
           GNAME.c_str(), NF, DIMG, RANK, (unsigned long long)CASF, CM, (unsigned long long)PC, (unsigned long long)P, (unsigned long long)seed, q1_selftest());
    return 0;
  }
  if (argc < 8) { fprintf(stderr, "need t1 t2 t3 n+ n-\n"); return 1; }
  int w[5]; for (int i = 0; i < 5; i++) w[i] = atoi(argv[3 + i]);
  if (w[0] < w[1] || w[1] < w[2] || w[2] < 0 || w[3] < w[4] || w[4] < 0) { fprintf(stderr, "sector must be dominant\n"); return 1; }
  Sector src = make_sector(w);
  if (mode == "dims") {
    int ylo = argc > 8 ? atoi(argv[8]) : 1, yhi = argc > 9 ? atoi(argv[9]) : src.Ymax;
    init_all({src});
    print_head("dims", w, src, 0); printf(",\"d_Y\":{");
    long chi = 0; bool first = true;
    for (int Y = ylo; Y <= std::min(yhi, src.Ymax); Y++) {
      i64 d = 0; size_t nsk = 0; i64 mx = 0; for (auto &kv : hw_by_skel(src, Y)) { d += kv.second; nsk++; mx = std::max(mx, kv.second); }
      if (getenv("DIMS_SKEL")) fprintf(stderr, "#Y=%d skeletons %zu max_skel_hw %lld\n", Y, nsk, (long long)mx);
      if (!d) continue;
      int F = ((Y + src.Theta) % 2 + 2) % 2; chi += F ? -d : d;
      printf("%s\"%d\":%lld", first ? "" : ",", Y, (long long)d); first = false;
    }
    printf("},\"chi\":%ld,\"seconds\":%.2f}\n", chi, now() - T0);
    return 0;
  }
  if (mode == "chain") {
    if (argc < 10) { fprintf(stderr, "chain needs Ylo Yhi\n"); return 1; }
    int ylo = std::max(1, atoi(argv[8])), yhi = std::min(src.Ymax, atoi(argv[9]));
    u64 seed = setseed(10);
    init_all({src});
    bool want_rho = !getenv("A_NORHO");
    int ya = std::max(1, ylo - 1), ytop = std::min(src.Ymax, yhi + 1);
    std::map<int, Lev> L; std::map<int, size_t> rout, rho;
    for (int Y = ya; Y <= ytop; Y++) build_lev(src, L[Y], Y);
    for (int Y = ya; Y <= yhi; Y++) {
      size_t nc; col_offsets(L[Y], nc);
      Mtx M; if (Y < src.Ymax) M = q0_rows(src, L[Y], L[Y + 1]);
      if (Y >= ylo && want_rho && nc) {
        Cohom C; C.sec = src; C.Y = Y; C.L0 = L[Y]; C.dim = L[Y].dim; C.c0 = col_offsets(C.L0, C.nc);
        if (Y < src.Ymax) C.K = kernel_of(M, nc); else for (size_t i = 0; i < nc; i++) { std::vector<u32> v(nc, 0); v[i] = 1; C.K.push_back(v); }
        rout[Y] = nc - C.K.size();
        compute_res(C);
        IncEch E; E.N = (int)C.K.size();
        for (auto &rb : C.res) { size_t o = C.c0[rb.si], nb = C.L0.sk[rb.si].gens.size();
          for (size_t k = 0; k < rb.val.size(); k++) { std::vector<u64> row(C.K.size(), 0); for (size_t z = 0; z < C.K.size(); z++) { u64 s = 0; int c = 0; for (size_t i = 0; i < nb; i++) { s += (u64)rb.val[k][i] * C.K[z][o + i]; if (++c == 3) { s = red(s); c = 0; } } row[z] = red(s); } E.insert(row); } }
        rho[Y] = E.R.size();
      } else rout[Y] = (Y < src.Ymax && nc) ? rank_of(M, nc) : 0;
      if (VERB) fprintf(stderr, "#Y=%d rank_out %zu  %.1fs rss %.0fMB\n", Y, rout[Y], now() - T0, maxrss_mb());
    }
    print_head("chain", w, src, seed); printf(",\"blocks\":[");
    bool first = true;
    for (int Y = ylo; Y <= yhi; Y++) {
      size_t d = L[Y].dim, ri = (Y >= 2 && rout.count(Y - 1)) ? rout[Y - 1] : 0, ro = rout[Y];
      long h = (long)d - (long)ri - (long)ro;
      int F = ((Y + src.Theta) % 2 + 2) % 2;
      printf("%s{\"Y\":%d,\"F\":%d,\"dim\":%zu,\"rank_in\":%zu,\"rank_out\":%zu,\"h\":%ld", first ? "" : ",", Y, F, d, ri, ro, h);
      if (rho.count(Y)) printf(",\"rho\":%zu", rho[Y]);
      printf("}"); first = false;
    }
    printf("],\"seconds\":%.1f,\"evals\":%zu,\"maxrss_MB\":%.0f}\n", now() - T0, (size_t)NEVAL, maxrss_mb());
    return 0;
  }
  if (mode == "d1" || mode == "q1check" || mode == "anticomm" || mode == "ratio") return q1_main(mode, argc, argv, w, src, setseed, T0);
  fprintf(stderr, "unknown mode\n"); return 1;
}
