// Q0 on G2-invariant BMN cochains, with actual adjoint matrices in the 7.
#include "core.h"
#include <stdexcept>

static std::vector<Mat> g2_basis;

static void set_g2_basis() {
  // The compact octonion three-form. Each pair of axes occurs once.
  const int triples[][4] = {{0,1,2,1},{0,3,4,1},{0,5,6,1},{1,3,5,1},
                            {1,4,6,-1},{2,3,6,-1},{2,4,5,-1}};
  int phi[7][7][7] = {};
  for (const auto &t : triples) {
    for (int q = 0; q < 3; ++q) {
      int i = t[q], j = t[(q+1)%3], k = t[(q+2)%3];
      phi[i][j][k] = t[3];
      phi[j][i][k] = -t[3];
    }
  }
  for (int k = 0; k < 7; ++k) {
    std::vector<std::pair<int,int>> pairs;
    for (int i = 0; i < 7; ++i)
      for (int j = i+1; j < 7; ++j)
        if (phi[i][j][k]) pairs.emplace_back(i,j);
    if (pairs.size() != 3) throw std::logic_error("invalid octonion triples");
    for (int q = 0; q < 2; ++q) {
      Mat matrix{};
      for (int h : {q, 2}) {
        auto [i,j] = pairs[h];
        int c = (h == q ? 1 : -1)*phi[i][j][k];
        matrix[7*i+j] = md(c);
        matrix[7*j+i] = md(-c);
      }
      // Check the stabilizer equations themselves, not only the seven contractions.
      for (int i = 0; i < 7; ++i)
        for (int j = i+1; j < 7; ++j)
          for (int l = j+1; l < 7; ++l) {
            i64 derivative = 0;
            for (int h = 0; h < 7; ++h)
              derivative += i64(matrix[7*h+i])*phi[h][j][l]
                          + i64(matrix[7*h+j])*phi[i][h][l]
                          + i64(matrix[7*h+l])*phi[i][j][h];
            if (md(derivative)) throw std::logic_error("matrix does not preserve three-form");
          }
      g2_basis.push_back(matrix);
    }
  }
}

static Mat g2_probe() {
  Mat matrix{};
  for (const Mat &basis : g2_basis) {
    u64 coefficient = rnd()%P;
    for (int k = 0; k < 49; ++k) {
      if (basis[k] == 1) matrix[k] = addm(matrix[k], coefficient);
      else if (basis[k]) matrix[k] = subm(matrix[k], coefficient);
    }
  }
  return matrix;
}

// Only the classical trace evaluator is used. All its non-Cartan probes use
// the stabilizer above. Cartan evaluation is unnecessary for the adjacent ranks.
#define rand_alg g2_probe
#include "trace.h"
#undef rand_alg
#include "bmn_sampler.h"
#define random_gen bmn_gen
#include "levels.h"
#undef random_gen

// At the last length only the image of Q0 is used. Evaluations on the target
// are linear functionals, whether or not they separate the whole target.
static void terminal_evaluations(const Sector &sec, Lev &level, int length, size_t image_cap) {
  level = Lev{};
  level.Y = length;
  if (length > sec.Ymax) return;
  for (const auto &entry : hw_by_skel(sec, length)) {
    Skel block = parse_skel(entry.first, length, entry.second);
    level.exact += block.hw;
    level.dim += block.hw;
    block.K = static_cast<int>(std::min<size_t>(block.hw, image_cap) + MARGIN);
    add_probes(block, block.pr, block.K, false);
    for (int k = 0; k < block.K; ++k) block.piv.push_back(k);
    level.sk.push_back(std::move(block));
  }
}

int main(int argc, char **argv) {
  bool probe_last = argc == 10 && std::string(argv[9]) == "--probe-last";
  if (argc != 9 && !probe_last) {
    fputs("usage: g2_chain w1 w2 w3 first_length last_length seed threads prime [--probe-last]\n", stderr);
    return 1;
  }
  int w[5] = {std::stoi(argv[1]), std::stoi(argv[2]), std::stoi(argv[3]), 0, 0};
  int lo = std::stoi(argv[4]), hi = std::stoi(argv[5]);
  u64 seed = std::stoull(argv[6]);
  NTHREADS = std::stoi(argv[7]);
  P = std::stoull(argv[8]);
  if ((P != 2147483647 && P != 1000003 && P != 1000033) || NTHREADS < 1 ||
      w[0] < w[1] || w[1] < w[2] || w[2] < 0 || lo < 1 || hi < lo) return 2;
  MERS = P == 2147483647;
  RS = seed;
  NF = 7; NN = 49; RANK = 2; DIMG = 14; GNAME = "G2";
  set_mm(); set_bmm(); set_g2_basis();
  // In the torus x,y,-x-y the short roots are the six nonzero weights
  // of the 7, and the six long roots are their A2 differences.
  POSR = {{1,0,0},{0,1,0},{1,1,0},{1,-1,0},{2,1,0},{1,2,0}};
  ADJW = {{0,0,0},{0,0,0}};
  for (auto r : POSR) {
    ADJW.push_back(r); ADJW.push_back({-r[0],-r[1],0});
  }
  Sector sec = make_sector(w);
  if (hi > sec.Ymax) return 2;
  std::copy(w,w+5,BOX);
  build_char_grid(sec.Ymax+1); build_letters(); build_q(); build_raise();
  HWCHECK = true;
  Lev prev, curr;
  size_t incoming = 0;
  std::vector<std::array<size_t,4>> blocks;
  for (int y = std::max(1,lo-1); y <= hi+1; ++y) {
    if (probe_last && y == hi+1) terminal_evaluations(sec,curr,y,prev.dim-incoming);
    else build_lev(sec,curr,y);
    size_t outgoing = prev.dim && curr.dim ? rank_of(q0_rows(sec,prev,curr),prev.dim) : 0;
    if (y > lo) {
      if (outgoing+incoming > prev.dim) throw std::logic_error("incompatible adjacent Q0 ranks");
      blocks.push_back({size_t(y-1),prev.dim,incoming,outgoing});
      fprintf(stderr,"Y=%d dim=%zu in=%zu out=%zu\n", y-1,prev.dim,incoming,outgoing);
    }
    incoming = outgoing;
    prev = std::move(curr);
  }
  if (probe_last) fputs("{\"terminal_probes_only\":true,", stdout);
  else fputc('{', stdout);
  printf("\"weight\":[%d,%d,%d],\"prime\":%llu,\"seed\":%llu,\"blocks\":[",
         w[0],w[1],w[2],(unsigned long long)P,(unsigned long long)seed);
  for (size_t k=0;k<blocks.size();++k) {
    const auto &b=blocks[k];
    printf("%s{\"length\":%zu,\"dimension\":%zu,\"incoming\":%zu,\"outgoing\":%zu}",
           k?",":"",b[0],b[1],b[2],b[3]);
  }
  puts("]}");
}
