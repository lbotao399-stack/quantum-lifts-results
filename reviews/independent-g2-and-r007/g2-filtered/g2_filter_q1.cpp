// Apply the trace-dual G2 pair map to selected cycles or correcting primitives.
#define FILTERED_QUANTUM
#include "g2_filtered_sector.cpp"
#include "g2_filter_kernel.h"
#include <iostream>

static Mtx load_vectors(const fs::path &path,size_t width) {
  std::ifstream f(path,std::ios::binary);
  f.exceptions(std::ios::failbit|std::ios::badbit);
  u32 header[4]; f.read(reinterpret_cast<char *>(header),sizeof(header));
  if(header[0]!=0x47325130 || header[1]!=P || header[3]!=width) throw std::logic_error("vectors");
  Mtx v(header[2],std::vector<u32>(width));
  for(auto &r:v) f.read(reinterpret_cast<char *>(r.data()),width*sizeof(u32));
  return v;
}

int main(int argc,char **argv) {
  if(argc!=9) return 1;
  fs::path dir=argv[1]; int grade=std::stoi(argv[2]),y=std::stoi(argv[3]);
  std::string vector=argv[4]; int base=std::stoi(argv[5]); P=std::stoul(argv[6]);
  MERS=P==2147483647; NTHREADS=std::stoi(argv[7]); bool dual=std::stoi(argv[8]);
  auto family=charge_family(); initialize_exceptional(family); set_colour_maps(); init_q1();
  const Sector &s=family.at(grade+1), &t=family.at(grade+2);
  Lev source,target; fixed_level(s,y,base,source); fixed_level(t,y,base,target);
  auto name=[&](std::string tag,int grade_id) {
    return dir/(tag+"_"+std::to_string(grade_id)+"_"+std::to_string(y)+".u32");
  };
  auto inputs=load_vectors(name(vector,grade),source.dim);
  Mtx tests;
  if(dual) tests=load_vectors(name("dual",grade+1),target.dim);
  std::vector<std::tuple<size_t,size_t,size_t>> probes;
  size_t offset=0;
  for(size_t k=0;k<target.sk.size();k++) {
    for(int j:target.sk[k].piv) {
      if(!dual || std::any_of(tests.begin(),tests.end(),[&](const auto &r){return r[offset];}))
        probes.push_back({offset,k,j});
      offset++;
    }
  }
  size_t source_width=0;
  auto origins=col_offsets(source,source_width);
  Mtx result(target.dim,std::vector<u32>(inputs.size()));
  std::atomic<size_t> completed{0},compared{0},compared_nonzero{0};
  // Check both colour routings on nonzero generator values before retaining a pairing.
  for(size_t q=0;q<probes.size();q++) {
    if(compared_nonzero>=6) break;
    auto [r,k,j]=probes[q]; const Skel &v=target.sk[k];
    auto lf=[&](int letter,GMat &z){letter_at(v,v.pr,j,letter,z);};
    for(const Skel &u:source.sk) {
      if(!skel_connected(u,v)) continue;
      for(size_t i=0;i<std::min<size_t>(u.gens.size(),20);i++) {
        u64 fast=q1_eval_fast(s,u,u.gens[i],v,lf,v.top_te);
        if(!fast && compared>=24) continue;
        if(fast!=q1_eval_brute(s,u,u.gens[i],v,lf,v.top_te,false))
          throw std::logic_error("direct colour contraction differs");
        compared++; if(fast) compared_nonzero++;
        if(compared_nonzero>=6) break;
      }
    }
  }
  std::cerr<<"direct generator tests "<<compared<<" nonzero "<<compared_nonzero<<"\n";
  par_for(probes.size(),[&](size_t qi) {
    auto [r,k,j]=probes[qi]; const Skel &v=target.sk[k];
    auto lf=[&](int letter,GMat &z){letter_at(v,v.pr,j,letter,z);};
    for(size_t g=0;g<source.sk.size();g++) {
      const auto &u=source.sk[g];
      if(!skel_connected(u,v)) continue;
      for(size_t i=0;i<u.gens.size();i++) {
        size_t index=origins[g]+i;
        if(std::none_of(inputs.begin(),inputs.end(),[&](const auto &c){return c[index];})) continue;
        auto value=q1_eval_fast(s,u,u.gens[i],v,lf,v.top_te);
        if(!value) continue;
        for(size_t h=0;h<inputs.size();h++)
          result[r][h]=addm(result[r][h],mulm(value,inputs[h][index]));
      }
    }
    size_t done=++completed;
    if(done%25==0) fprintf(stderr,"quantum probes %zu/%zu\n",done,probes.size());
  });
  Mtx output=dual ? mul_abt(tests,transpose(result,inputs.size()),target.dim):result;
  store_matrix(name(std::string(dual?"paired_":"acted_")+vector,grade),output,inputs.size());
  std::cout<<"finished "<<grade<<" "<<y<<" rows "<<output.size()<<" columns "<<inputs.size()<<"\n";
}
