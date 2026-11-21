// linalg.h -- dense linear algebra mod P (u32 storage), multithreaded row elimination.
#pragma once
#include "core.h"
typedef std::vector<std::vector<u32>> Mtx;

// reduced row echelon form in place; returns pivot columns; M trimmed to rank rows
static std::vector<int> rref(Mtx &M, size_t n) {
  std::vector<int> piv; size_t r = 0;
  for (size_t c = 0; c < n && r < M.size(); c++) {
    size_t p = r; while (p < M.size() && !M[p][c]) p++;
    if (p == M.size()) continue;
    std::swap(M[p], M[r]);
    u64 inv = invm(M[r][c]); for (auto &x : M[r]) x = (u32)mulm(x, inv);
    const std::vector<u32> &pr = M[r];
    std::vector<u32> nz; for (size_t j = c; j < n; j++) if (pr[j]) nz.push_back((u32)j);
    par_for(M.size(), [&](size_t i) {
      if (i == r || !M[i][c]) return; u64 a = P - M[i][c]; auto &row = M[i];
      for (u32 j : nz) row[j] = (u32)red(row[j] + a * pr[j]);
    });
    piv.push_back((int)c); r++;
  }
  M.resize(r);
  return piv;
}
// forward elimination only (row echelon, not reduced): faster rank
static size_t rank_of(Mtx M, size_t n) {
  size_t r = 0;
  for (size_t c = 0; c < n && r < M.size(); c++) {
    size_t p = r; while (p < M.size() && !M[p][c]) p++;
    if (p == M.size()) continue;
    std::swap(M[p], M[r]);
    u64 inv = invm(M[r][c]); for (size_t j = c; j < n; j++) M[r][j] = (u32)mulm(M[r][j], inv);
    const std::vector<u32> &pr = M[r];
    std::vector<u32> nz; for (size_t j = c; j < n; j++) if (pr[j]) nz.push_back((u32)j);
    size_t base = r + 1, cnt = M.size() - base;
    par_for(cnt, [&](size_t ii) {
      auto &row = M[base + ii]; if (!row[c]) return; u64 a = P - row[c];
      for (u32 j : nz) row[j] = (u32)red(row[j] + a * pr[j]);
    });
    std::vector<u32>().swap(M[r]);   // row no longer needed
    r++;
  }
  return r;
}
// basis of {x : M x = 0} (x of length n)
static Mtx kernel_of(Mtx M, size_t n) {
  std::vector<int> piv = rref(M, n);
  std::vector<char> isp(n, 0); for (int c : piv) isp[c] = 1;
  Mtx K;
  for (size_t f = 0; f < n; f++) if (!isp[f]) {
    std::vector<u32> v(n, 0); v[f] = 1;
    for (size_t r = 0; r < piv.size(); r++) v[piv[r]] = (u32)negm(M[r][f]);
    K.push_back(std::move(v));
  }
  return K;
}
// C = A * B^T : C[i][j] = sum_k A[i][k] Bt[j][k]
static Mtx mul_abt(const Mtx &A, const Mtx &Bt, size_t n) {
  Mtx C(A.size(), std::vector<u32>(Bt.size(), 0));
  par_for(A.size(), [&](size_t i) {
    const u32 *a = A[i].data();
    for (size_t j = 0; j < Bt.size(); j++) { const u32 *b = Bt[j].data(); u64 s = 0; int c = 0; for (size_t k = 0; k < n; k++) { s += (u64)a[k] * b[k]; if (++c == 3) { s = red(s); c = 0; } } C[i][j] = (u32)red(s); }
  });
  return C;
}
static Mtx transpose(const Mtx &A, size_t nc) {
  Mtx T(nc, std::vector<u32>(A.size(), 0));
  for (size_t i = 0; i < A.size(); i++) for (size_t j = 0; j < nc; j++) T[j][i] = A[i][j];
  return T;
}
static std::vector<u32> to32(const std::vector<u64> &v) { std::vector<u32> w(v.size()); for (size_t i = 0; i < v.size(); i++) w[i] = (u32)v[i]; return w; }
