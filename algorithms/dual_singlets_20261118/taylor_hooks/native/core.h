// core.h -- modular arithmetic, gauge algebra (SU(N), SO(2n+1), Sp(n)), exact colour-singlet counts.
// Part of the so7/sp3 trace-word engine (engine.cpp).  All matrices act on the fundamental (N_F x N_F):
//   SU(N): N_F = N, traceless;  SO(2n+1): N_F = 2n+1, A^T eta + eta A = 0, eta_{i,ibar} = 1 (ibar = N_F-1-i);
//   Sp(n): N_F = 2n, A^T Om + Om A = 0, Om_{i,ibar} = s_i (s_i = +1 for i < n, -1 for i >= n).
// In these bases the Cartan subalgebra is diagonal.  kappa_{RS} = tr_F(t_R t_S); t^R = kappa-dual basis.
#pragma once
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <map>
#include <mutex>
#include <numeric>
#include <set>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>
#include <sys/resource.h>

typedef uint64_t u64; typedef int64_t i64; typedef uint32_t u32;

// ---------------- rank prime ----------------
static u64 P = 2147483647ULL; static bool MERS = true;
static inline u64 red(u64 x) { if (MERS) { x = (x & 2147483647ULL) + (x >> 31); x = (x & 2147483647ULL) + (x >> 31); return x >= 2147483647ULL ? x - 2147483647ULL : x; } return x % P; }
static inline u64 md(i64 x) { i64 r = x % (i64)P; return r < 0 ? r + (i64)P : (u64)r; }
static inline u64 mulm(u64 a, u64 b) { return red(a * b); }
static inline u64 addm(u64 a, u64 b) { u64 s = a + b; return s >= P ? s - P : s; }
static inline u64 subm(u64 a, u64 b) { return a >= b ? a - b : a + P - b; }
static inline u64 negm(u64 a) { return a ? P - a : 0; }
static u64 powm(u64 a, u64 e) { u64 r = 1; a %= P; while (e) { if (e & 1) r = mulm(r, a); a = mulm(a, a); e >>= 1; } return r; }
static inline u64 invm(u64 a) { return powm(a, P - 2); }
static u64 RS = 0x9E3779B97F4A7C15ULL;
static u64 rnd() { u64 z = (RS += 0x9E3779B97F4A7C15ULL); z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; return z ^ (z >> 31); }
static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
static double maxrss_mb() { struct rusage ru; getrusage(RUSAGE_SELF, &ru);
#ifdef __APPLE__
  return ru.ru_maxrss / 1048576.0;
#else
  return ru.ru_maxrss / 1024.0;
#endif
}
static int NTHREADS = 1;
#include <condition_variable>
// persistent thread pool: par_for(n, f) runs f(i), i < n, dynamically scheduled; nested calls run serially
struct Pool {
  std::vector<std::thread> th; std::mutex mu; std::condition_variable cv, done_cv;
  std::function<void(size_t)> job; size_t n = 0; std::atomic<size_t> next{0}; int active = 0; u64 gen = 0; bool stop = false;
  void start(int nt) {
    for (int t = 0; t < nt; t++) th.emplace_back([this]() {
      u64 seen = 0;
      while (true) {
        std::unique_lock<std::mutex> lk(mu); cv.wait(lk, [&]() { return stop || gen != seen; }); if (stop) return; seen = gen; active++; lk.unlock();
        tl_in() = true; size_t i; while ((i = next++) < n) job(i); tl_in() = false;
        lk.lock(); if (--active == 0) done_cv.notify_all();
      }
    });
  }
  static bool &tl_in() { static thread_local bool v = false; return v; }
  void run(size_t nn, const std::function<void(size_t)> &f) {
    { std::unique_lock<std::mutex> lk(mu); done_cv.wait(lk, [&]() { return active == 0; }); job = f; n = nn; next = 0; gen++; }
    cv.notify_all();
    tl_in() = true; size_t i; while ((i = next++) < n) job(i); tl_in() = false;
    std::unique_lock<std::mutex> lk(mu); done_cv.wait(lk, [&]() { return active == 0 && next >= n; });
  }
  ~Pool() { { std::lock_guard<std::mutex> lk(mu); stop = true; } cv.notify_all(); for (auto &x : th) x.join(); }
};
static Pool *POOL = nullptr;
static std::mutex POOLMU;
template <class F> static void par_for(size_t n, F f) {
  if (NTHREADS <= 1 || n < 2 || Pool::tl_in()) { for (size_t i = 0; i < n; i++) f(i); return; }
  if (!POOL) { POOL = new Pool(); POOL->start(NTHREADS - 1); }
  std::lock_guard<std::mutex> lk(POOLMU);
  POOL->run(n, std::function<void(size_t)>(f));
}

// ---------------- gauge group ----------------
enum { G_SU = 0, G_SO = 1, G_SP = 2 };
static int GT = G_SO, GR = 3, NF = 7, NN = 49, RANK = 3, DIMG = 21;
static std::string GNAME = "SO7";
static const int MAXNF = 7;
typedef std::array<u32, MAXNF * MAXNF> Mat;   // row-major, stride NF, first NN entries used
static int SG[MAXNF];                          // Sp signs s_i

static bool parse_gauge(const std::string &g) {
  if (g.size() < 3) return false;
  std::string t = g.substr(0, 2); int k = atoi(g.c_str() + 2);
  if (t == "SU" && k >= 2 && k <= 7) { GT = G_SU; GR = k; NF = k; RANK = k - 1; DIMG = k * k - 1; }
  else if (t == "SO" && k % 2 == 1 && k >= 3 && k <= 7) { GT = G_SO; GR = (k - 1) / 2; NF = k; RANK = GR; DIMG = GR * (2 * GR + 1); }
  else if ((t == "SP" || t == "Sp") && k >= 1 && k <= 3) { GT = G_SP; GR = k; NF = 2 * k; RANK = k; DIMG = k * (2 * k + 1); }
  else return false;
  NN = NF * NF; GNAME = g; GNAME[1] = (GT == G_SP) ? 'p' : GNAME[1]; if (GT == G_SP) GNAME = "Sp" + std::to_string(k);
  for (int i = 0; i < NF; i++) SG[i] = (i < NF / 2) ? 1 : -1;
  return true;
}
static inline int bar(int i) { return NF - 1 - i; }
static void set_mm();

// ---- matrix helpers ----
static inline void mzero(Mat &M) { for (int i = 0; i < NN; i++) M[i] = 0; }
// C = A B, specialised on the fundamental dimension (products < 2^62; fold every 3 terms)
template <int N> static void mmulT(const Mat &A, const Mat &B, Mat &C) {
  for (int i = 0; i < N; i++) {
    u64 acc[N];
    for (int k = 0; k < N; k++) acc[k] = (u64)A[i * N] * B[k];
    for (int j = 1; j < N; j++) {
      u64 a = A[i * N + j];
      for (int k = 0; k < N; k++) acc[k] += a * B[j * N + k];
      if (j % 3 == 2) for (int k = 0; k < N; k++) acc[k] = MERS ? (acc[k] & 2147483647ULL) + (acc[k] >> 31) : acc[k] % P;
    }
    for (int k = 0; k < N; k++) C[i * N + k] = (u32)red(acc[k]);
  }
}
template <int N> static u64 mtrmulT(const Mat &A, const Mat &B) {
  u64 s = 0;
  for (int i = 0; i < N; i++) { u64 t = 0; for (int j = 0; j < N; j++) { t += (u64)A[i * N + j] * B[j * N + i]; if (j % 3 == 2) t = MERS ? (t & 2147483647ULL) + (t >> 31) : t % P; } s += red(t); }
  return red(s);
}
typedef void (*MMulF)(const Mat &, const Mat &, Mat &); typedef u64 (*MTrF)(const Mat &, const Mat &);
static MMulF MMUL = mmulT<7>; static MTrF MTRMUL = mtrmulT<7>;
static void set_mm() {
  switch (NF) { case 2: MMUL = mmulT<2>; MTRMUL = mtrmulT<2>; break; case 3: MMUL = mmulT<3>; MTRMUL = mtrmulT<3>; break; case 4: MMUL = mmulT<4>; MTRMUL = mtrmulT<4>; break;
    case 5: MMUL = mmulT<5>; MTRMUL = mtrmulT<5>; break; case 6: MMUL = mmulT<6>; MTRMUL = mtrmulT<6>; break; default: MMUL = mmulT<7>; MTRMUL = mtrmulT<7>; }
}
static inline void mmul(const Mat &A, const Mat &B, Mat &C) { MMUL(A, B, C); }
static inline u64 mtrace(const Mat &M) { u64 s = 0; for (int i = 0; i < NF; i++) s += M[i * NF + i]; return red(s); }
static inline u64 mtrmul(const Mat &A, const Mat &B) { return MTRMUL(A, B); }
static inline void mneg(Mat &M) { for (int i = 0; i < NN; i++) M[i] = (u32)negm(M[i]); }
static inline void mscale(Mat &M, u64 c) { for (int i = 0; i < NN; i++) M[i] = (u32)mulm(M[i], c); }
static inline void maddto(Mat &C, const Mat &A) { for (int i = 0; i < NN; i++) C[i] = (u32)addm(C[i], A[i]); }
static inline void mtranspose(const Mat &A, Mat &T) { for (int i = 0; i < NF; i++) for (int j = 0; j < NF; j++) T[i * NF + j] = A[j * NF + i]; }
static inline bool miszero(const Mat &M) { for (int i = 0; i < NN; i++) if (M[i]) return false; return true; }

// projection of an arbitrary matrix M onto the algebra (up to the factor 1/2 for SO, Sp)
static Mat to_alg(const Mat &M) {
  Mat A; mzero(A);
  if (GT == G_SU) { A = M; u64 s = 0; for (int i = 0; i < NF - 1; i++) s += M[i * NF + i]; A[NN - 1] = (u32)negm(red(s)); return A; }
  for (int a = 0; a < NF; a++) for (int b = 0; b < NF; b++) {
    u64 x = M[a * NF + b], y = M[bar(b) * NF + bar(a)];
    if (GT == G_SO) A[a * NF + b] = (u32)subm(x, y);
    else A[a * NF + b] = (u32)((SG[a] * SG[b] > 0) ? subm(x, y) : addm(x, y));
  }
  return A;
}
static bool in_alg(const Mat &A) {
  if (GT == G_SU) return mtrace(A) == 0;
  for (int a = 0; a < NF; a++) for (int b = 0; b < NF; b++) {
    u64 x = A[a * NF + b], y = A[bar(b) * NF + bar(a)];   // need y = -s_a s_b x (Sp) or y = -x (SO)
    if (GT == G_SO) { if (addm(x, y)) return false; }
    else { if (SG[a] * SG[b] > 0 ? addm(x, y) : subm(x, y)) return false; }
  }
  return true;
}
static Mat rand_alg() { Mat M; mzero(M); for (int i = 0; i < NN; i++) M[i] = (u32)(rnd() % P); return to_alg(M); }
static Mat rand_cartan() {
  Mat M; mzero(M);
  if (GT == G_SU) { u64 s = 0; for (int i = 0; i < NF - 1; i++) { u64 v = rnd() % P; M[i * NF + i] = (u32)v; s += v; } M[NN - 1] = (u32)negm(red(s)); return M; }
  for (int i = 0; i < GR; i++) { u64 v = rnd() % P; M[i * NF + i] = (u32)v; M[bar(i) * NF + bar(i)] = (u32)negm(v); }
  return M;
}

// ---- algebra basis, dual basis, Casimir tensor ----
static std::vector<Mat> TB, TD;     // t_R, t^R (R < DIMG)
static u64 CC1, CC2, CC3; static Mat FM; static bool HAVEF = false;   // sum_R (t_R)_ab (t^R)_cd = c1 d_ad d_bc + c2 F_ac F_bd + c3 d_ab d_cd
static std::vector<u64> CASD;       // CASD[i*NF+j] = sum_R (t_R)_ij (t^R)_ji
static u64 CASF;                    // C_F: sum_R t_R t^R = C_F 1
static std::vector<std::vector<u64>> inv_matrix(std::vector<std::vector<u64>> A) {
  size_t n = A.size();
  for (size_t i = 0; i < n; i++) { A[i].resize(2 * n, 0); A[i][n + i] = 1; }
  for (size_t c = 0; c < n; c++) {
    size_t p = c; while (p < n && !A[p][c]) p++;
    if (p == n) { fprintf(stderr, "singular Gram matrix\n"); exit(20); }
    std::swap(A[p], A[c]); u64 iv = invm(A[c][c]); for (auto &x : A[c]) x = mulm(x, iv);
    for (size_t r = 0; r < n; r++) if (r != c && A[r][c]) { u64 a = P - A[r][c]; for (size_t k = 0; k < 2 * n; k++) A[r][k] = red(A[r][k] + a * A[c][k]); }
  }
  std::vector<std::vector<u64>> I(n, std::vector<u64>(n)); for (size_t i = 0; i < n; i++) for (size_t j = 0; j < n; j++) I[i][j] = A[i][n + j];
  return I;
}
static void build_algebra() {
  TB.clear(); TD.clear();
  // spanning set: projections of matrix units; greedy independent selection
  std::vector<std::vector<u64>> ech; std::vector<int> piv;
  for (int i = 0; i < NF; i++) for (int j = 0; j < NF; j++) {
    Mat E; mzero(E); E[i * NF + j] = 1;
    if (GT == G_SU) { if (i == j) { if (i == NF - 1) continue; E[(i + 1) * NF + i + 1] = (u32)(P - 1); } }
    Mat A = (GT == G_SU) ? E : to_alg(E);
    std::vector<u64> v(NN); for (int q = 0; q < NN; q++) v[q] = A[q];
    for (size_t r = 0; r < ech.size(); r++) { u64 a = v[piv[r]]; if (!a) continue; u64 na = P - a; for (int q = 0; q < NN; q++) v[q] = red(v[q] + na * ech[r][q]); }
    int p = -1; for (int q = 0; q < NN; q++) if (v[q]) { p = q; break; }
    if (p < 0) continue;
    u64 iv = invm(v[p]); for (auto &x : v) x = mulm(x, iv); ech.push_back(v); piv.push_back(p); TB.push_back(A);
  }
  if ((int)TB.size() != DIMG) { fprintf(stderr, "algebra dimension %zu != %d\n", TB.size(), DIMG); exit(20); }
  for (auto &A : TB) if (!in_alg(A)) { fprintf(stderr, "basis element outside algebra\n"); exit(20); }
  std::vector<std::vector<u64>> G(DIMG, std::vector<u64>(DIMG));
  for (int a = 0; a < DIMG; a++) for (int b = 0; b < DIMG; b++) G[a][b] = mtrmul(TB[a], TB[b]);
  std::vector<std::vector<u64>> Gi = inv_matrix(G);
  for (int a = 0; a < DIMG; a++) { Mat D; mzero(D); for (int b = 0; b < DIMG; b++) if (Gi[a][b]) { Mat X = TB[b]; mscale(X, Gi[a][b]); maddto(D, X); } TD.push_back(D); }
  // Casimir tensor and its closed form
  std::vector<u64> C((size_t)NN * NN, 0);
  for (int R = 0; R < DIMG; R++) for (int x = 0; x < NN; x++) if (TB[R][x]) for (int y = 0; y < NN; y++) if (TD[R][y]) C[(size_t)x * NN + y] = addm(C[(size_t)x * NN + y], mulm(TB[R][x], TD[R][y]));
  u64 h = invm(2);
  mzero(FM); HAVEF = (GT != G_SU);
  if (GT == G_SU) { CC1 = 1; CC2 = 0; CC3 = negm(invm(NF)); }
  else { CC1 = h; CC2 = negm(h); CC3 = 0; for (int i = 0; i < NF; i++) FM[i * NF + bar(i)] = (u32)(GT == G_SO ? 1 : (SG[i] > 0 ? 1 : P - 1)); }
  for (int a = 0; a < NF; a++) for (int b = 0; b < NF; b++) for (int c = 0; c < NF; c++) for (int d = 0; d < NF; d++) {
    u64 v = 0;
    if (a == d && b == c) v = addm(v, CC1);
    if (HAVEF && FM[a * NF + c] && FM[b * NF + d]) v = addm(v, mulm(CC2, mulm(FM[a * NF + c], FM[b * NF + d])));
    if (a == b && c == d) v = addm(v, CC3);
    if (v != C[(size_t)(a * NF + b) * NN + c * NF + d]) { fprintf(stderr, "Casimir closed form mismatch at %d%d%d%d\n", a, b, c, d); exit(20); }
  }
  CASD.assign(NN, 0); for (int i = 0; i < NF; i++) for (int j = 0; j < NF; j++) CASD[i * NF + j] = C[(size_t)(i * NF + j) * NN + j * NF + i];
  Mat S; mzero(S); for (int R = 0; R < DIMG; R++) { Mat X; mmul(TB[R], TD[R], X); maddto(S, X); }
  CASF = S[0];
  for (int i = 0; i < NF; i++) for (int j = 0; j < NF; j++) if (S[i * NF + j] != (i == j ? CASF : 0)) { fprintf(stderr, "Casimir not scalar\n"); exit(20); }
  u64 expect = GT == G_SU ? mulm(md((i64)NF * NF - 1), invm(NF)) : GT == G_SO ? mulm(md(NF - 1), h) : mulm(md(NF + 1), h);
  if (CASF != expect) { fprintf(stderr, "C_F mismatch\n"); exit(20); }
  // closure: [rand, rand] in the algebra; Cartan probes commute and lie in the algebra
  for (int t = 0; t < 4; t++) { Mat A = rand_alg(), B = rand_alg(), X, Y; mmul(A, B, X); mmul(B, A, Y); for (int q = 0; q < NN; q++) X[q] = (u32)subm(X[q], Y[q]); if (!in_alg(X) || !in_alg(A)) { fprintf(stderr, "algebra not closed\n"); exit(20); } Mat H = rand_cartan(); if (!in_alg(H)) { fprintf(stderr, "Cartan probe outside algebra\n"); exit(20); } }
}

// ---------------- exact colour-singlet counts (Weyl integration on a finite torus grid mod PC) ----------------
// #Inv( (x)_i S^{k_i}/Lambda^{k_i}(adj) ) = CT_z[ prod_i h_{k_i}/e_{k_i}(z^{wt adj}) * prod_{a>0}(1 - z^{-a}) ]
// evaluated exactly as (1/M^r) sum over z in mu_M^r (M > max |exponent|), arithmetic mod a prime PC = 1 mod M.
static u64 PC = 0; static int CM = 0, NPT = 0, KMAXC = 0;
static std::vector<std::array<int, 3>> ADJW, POSR;
static std::vector<u64> CH, CE, CDEL;   // CH[k*NPT+pt], CE[k*NPT+pt], CDEL[pt]
static u64 mulc(u64 a, u64 b) { return (u64)((unsigned __int128)a * b % PC); }
static u64 powc(u64 a, u64 e) { u64 r = 1; a %= PC; while (e) { if (e & 1) r = mulc(r, a); a = mulc(a, a); e >>= 1; } return r; }
static bool isprime(u64 n) { if (n < 2) return false; for (u64 d = 2; d * d <= n; d++) if (n % d == 0) return false; return true; }
static void build_weights() {
  ADJW.clear(); POSR.clear();
  auto add = [&](std::vector<std::array<int, 3>> &v, int a, int b, int c) { v.push_back({a, b, c}); };
  if (GT == G_SU) {
    int N = NF;
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if (i != j) {
      std::array<int, 3> w = {0, 0, 0}; if (i < N - 1) w[i] += 1; if (j < N - 1) w[j] -= 1; ADJW.push_back(w);
      if (i < j) POSR.push_back(w);
    }
    for (int i = 0; i < N - 1; i++) add(ADJW, 0, 0, 0);
  } else {
    int n = GR, lr = (GT == G_SO) ? 1 : 2;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) for (int si = -1; si <= 1; si += 2) for (int sj = -1; sj <= 1; sj += 2) {
      std::array<int, 3> w = {0, 0, 0}; w[i] = si; w[j] = sj; ADJW.push_back(w); if (si > 0) POSR.push_back(w);
    }
    for (int i = 0; i < n; i++) for (int s = -1; s <= 1; s += 2) { std::array<int, 3> w = {0, 0, 0}; w[i] = s * lr; ADJW.push_back(w); if (s > 0) POSR.push_back(w); }
    for (int i = 0; i < n; i++) add(ADJW, 0, 0, 0);
  }
  if ((int)ADJW.size() != DIMG) { fprintf(stderr, "adjoint weights %zu != %d\n", ADJW.size(), DIMG); exit(20); }
}
static void build_char_grid(int Ymax) {
  build_weights();
  int amax = 0; for (auto &w : ADJW) for (int i = 0; i < RANK; i++) amax = std::max(amax, std::abs(w[i]));
  int dsum = 0; for (int i = 0; i < RANK; i++) { int s = 0; for (auto &a : POSR) s += std::abs(a[i]); dsum = std::max(dsum, s); }
  int D = amax * Ymax + dsum; CM = D + 1; KMAXC = Ymax + 1;
  NPT = 1; for (int i = 0; i < RANK; i++) NPT *= CM;
  // prime PC = 1 mod CM, PC < 2^31
  u64 t = ((1ULL << 31) - 1) / CM; while (true) { u64 q = t * CM + 1; if (isprime(q)) { PC = q; break; } t--; }
  // primitive CM-th root of unity
  std::vector<u64> fac; { u64 m = PC - 1; for (u64 d = 2; d * d <= m; d++) if (m % d == 0) { fac.push_back(d); while (m % d == 0) m /= d; } if (m > 1) fac.push_back(m); }
  u64 g = 2; while (true) { bool ok = true; for (u64 f : fac) if (powc(g, (PC - 1) / f) == 1) { ok = false; break; } if (ok) break; g++; }
  u64 om = powc(g, (PC - 1) / CM);
  std::vector<u64> opow(CM); opow[0] = 1; for (int i = 1; i < CM; i++) opow[i] = mulc(opow[i - 1], om);
  CH.assign((size_t)(KMAXC + 1) * NPT, 0); CE.assign((size_t)(KMAXC + 1) * NPT, 0); CDEL.assign(NPT, 0);
  std::vector<int> a(RANK);
  for (int pt = 0; pt < NPT; pt++) {
    int x = pt; for (int i = 0; i < RANK; i++) { a[i] = x % CM; x /= CM; }
    auto zp = [&](const std::array<int, 3> &w) { long e = 0; for (int i = 0; i < RANK; i++) e += (long)w[i] * a[i]; e %= CM; if (e < 0) e += CM; return opow[e]; };
    std::vector<u64> h(KMAXC + 1, 0), e(KMAXC + 1, 0); h[0] = e[0] = 1;
    for (auto &w : ADJW) {
      u64 y = zp(w);
      for (int k = 1; k <= KMAXC; k++) h[k] = (h[k] + mulc(y, h[k - 1])) % PC;
      for (int k = KMAXC; k >= 1; k--) e[k] = (e[k] + mulc(y, e[k - 1])) % PC;
    }
    for (int k = 0; k <= KMAXC; k++) { CH[(size_t)k * NPT + pt] = h[k]; CE[(size_t)k * NPT + pt] = e[k]; }
    u64 d = 1; for (auto &r : POSR) { std::array<int, 3> m = {-r[0], -r[1], -r[2]}; d = mulc(d, (1 + PC - zp(m)) % PC); }
    CDEL[pt] = d;
  }
}
static std::map<std::vector<std::pair<int, int>>, i64> PATC;
static std::mutex PATMU;
// pattern = list of (multiplicity k, odd); odd -> Lambda^k, even -> S^k
static i64 pattern_sing(std::vector<std::pair<int, int>> pat) {
  for (auto &x : pat) if (x.first == 1) x.second = 0;
  std::sort(pat.begin(), pat.end());
  { std::lock_guard<std::mutex> lk(PATMU); auto it = PATC.find(pat); if (it != PATC.end()) return it->second; }
  for (auto &x : pat) if (x.first > KMAXC) { fprintf(stderr, "pattern k %d > KMAXC %d\n", x.first, KMAXC); exit(21); }
  u64 s = 0;
  for (int pt = 0; pt < NPT; pt++) {
    u64 v = CDEL[pt];
    for (auto &x : pat) { v = mulc(v, (x.second ? CE : CH)[(size_t)x.first * NPT + pt]); if (!v) break; }
    s += v; if (s >= PC) s -= PC;
  }
  s = mulc(s, powc(powc(CM, RANK), PC - 2));
  if (s > (1ULL << 40)) { fprintf(stderr, "singlet count out of range (aliasing?)\n"); exit(21); }
  std::lock_guard<std::mutex> lk(PATMU); PATC[pat] = (i64)s; return (i64)s;
}
