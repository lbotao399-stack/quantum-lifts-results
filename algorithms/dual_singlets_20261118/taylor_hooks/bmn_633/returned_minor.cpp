// Test the returned product and its linear combination with omega in one USp(6) frame.
#define N4_EXPLICIT_IMAGE_HELPERS_ONLY
#include "../su2_levi_image.cpp"

int main(int argc, char **argv) {
  if (argc != 3 || !parse_gauge("Sp3")) return 1;
  P = 1000033; MERS = false; RS = 510731;
  if (const char *threads = std::getenv("PAR_THREADS")) NTHREADS = std::max(1, std::stoi(threads));
  set_mm(); set_bmm();
  int weights[5] = {5,0,0,3,3};
  const Sector sector = make_sector(weights);
  init_all({sector});
  const auto returned = read_image(argv[1],sector,9);
  const auto difference = read_image(argv[2],sector,9);
  Lev previous, target;
  build_lev(sector, previous, 8);
  build_lev(sector, target, 9);
  if (previous.dim!=2559 || target.dim!=2689) return 2;
  Mtx boundaries = q0_rows(sector,previous,target);
  std::vector<std::pair<size_t,int>> points;
  for (size_t j=0; j<target.sk.size(); ++j)
    for (int probe:target.sk[j].piv) points.emplace_back(j,probe);
  if (points.size()!=target.dim) return 3;
  std::vector<u32> product(points.size()), combination(points.size());
  par_for(points.size(),[&](size_t i) {
    const auto [s, probe] = points[i];
    product[i] = static_cast<u32>(evaluate_image(returned,target.sk[s],probe));
    combination[i] = static_cast<u32>(evaluate_image(difference,target.sk[s],probe));
  });
  print_head("symplectic_returned_product", weights,sector,510731);
  std::cout << ",\"length\":9,\"cochains\":[" << previous.dim << ',' << target.dim << "],\"incoming\":";
  print_minor(boundaries,previous.dim,true);
  for (size_t i=0; i<points.size(); ++i) boundaries[i].push_back(product[i]);
  std::cout << ",\"returned\":"; print_minor(boundaries,previous.dim+1,true);
  for (size_t i=0; i<points.size(); ++i) boundaries[i].back() = combination[i];
  std::cout << ",\"sum112\":"; print_minor(std::move(boundaries),previous.dim+1,true);
  std::cout << "}\n";
}
