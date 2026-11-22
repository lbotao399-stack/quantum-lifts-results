// Exact rank over the two prime fields. Every floating-point operation in a
// BLAS panel is on integers bounded by 128*(131-1)^2+131 < 2^24.
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>
extern "C" void cblas_sgemm(int, int, int, int, int, int, float, const float*, int,
                            const float*, int, float, float*, int);
constexpr int width = 128;
template<int p> static int residue(int x) {
    int y = x % p;
    return y < 0 ? y+p : y;
}
template<int p> static int rank_field(float* mtx, int rows, int cols) {
    std::vector<int> inv(p);
    inv[1]=1;
    for (int i=2; i<p; ++i) inv[i]=p-int((p/i)*inv[p%i]%p);
    int rank=0;
    auto row=[&](int i){return mtx+static_cast<std::size_t>(i)*cols;};
    for (int start=0; start<cols && rank<rows; start+=width) {
        int end=std::min(cols,start+width), size=end-start;
        int echelon[width][width]={}, piv[width]={}, v[width];
        std::vector<int> sites, indices;
        for(int i=rank;i<rows && int(sites.size())<size;++i) {
            for(int j=0;j<size;++j) v[j]=int(row(i)[start+j]);
            for(int j=0;j<size;++j) {
                if(!v[j]) continue;
                if(!piv[j]) {
                    piv[j]=1; sites.push_back(i); indices.push_back(start+j);
                    int a=inv[v[j]];
                    for(int k=j;k<size;++k) echelon[j][k]=residue<p>(v[k]*a);
                    break;
                }
                int a=v[j];
                for(int k=j+1;k<size;++k) v[k]=residue<p>(v[k]-a*echelon[j][k]);
            }
        }
        int n=static_cast<int>(sites.size());
        if(!n) continue;
        for(int i=0;i<n;++i) if(sites[i]!=rank+i)
            std::swap_ranges(row(rank+i)+start,row(rank+i)+cols,row(sites[i])+start);
        if(rank+n==rows || end==cols) {rank+=n; continue;}
        int b[width][2*width]={};
        for(int i=0;i<n;++i){
            for(int j=0;j<n;++j) b[i][j]=int(row(rank+i)[indices[j]]);
            b[i][n+i]=1;
        }
        for(int j=0;j<n;++j){
            int i=j;
            while(i<n && !b[i][j]) ++i;
            assert(i<n);
            if(i!=j) std::swap_ranges(b[i],b[i]+2*n,b[j]);
            int a=inv[b[j][j]];
            for(int k=j;k<2*n;++k) b[j][k]=residue<p>(a*b[j][k]);
            for(int t=0;t<n;++t) if(t!=j && b[t][j]){
                int coef=b[t][j];
                for(int k=j;k<2*n;++k) b[t][k]=residue<p>(b[t][k]-coef*b[j][k]);
            }
        }
        int remaining=rows-rank-n;
        std::vector<float> inverse(n*n), left(remaining*n), factors(remaining*n);
        for(int i=0;i<n;++i) for(int j=0;j<n;++j) inverse[i*n+j]=float(b[i][n+j]);
        for(int i=0;i<remaining;++i) for(int j=0;j<n;++j) left[i*n+j]=row(rank+n+i)[indices[j]];
        cblas_sgemm(101,111,111,remaining,n,n,1,left.data(),n,inverse.data(),n,0,factors.data(),n);
        for(float &v:factors) v=float(residue<p>(int(v)));
        cblas_sgemm(101,111,111,remaining,cols-end,n,-1,factors.data(),n,row(rank)+end,cols,
                    1,row(rank+n)+end,cols);
        #pragma omp parallel for if(remaining*(cols-end)>1000000)
        for(int i=rank+n;i<rows;++i) for(int j=end;j<cols;++j)
            row(i)[j]=float(residue<p>(int(row(i)[j])));
        rank+=n;
    }
    return rank;
}
extern "C" int panel_rank(float* a,int rows,int cols,int prime){
    if(prime==127) return rank_field<127>(a,rows,cols);
    if(prime==131) return rank_field<131>(a,rows,cols);
    return -1;
}
