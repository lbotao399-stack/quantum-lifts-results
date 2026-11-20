// Matrix prefixes and exterior subset ordering come from the archived trace evaluator.
#include "trace_frame_bmn_rank_threadsafe_v1.cpp"
#ifndef F4_COLOUR_COUNT
#define F4_COLOUR_COUNT 52
#endif
constexpr int colour_count = F4_COLOUR_COUNT;

extern "C" int f4_open(int k, int points, int degree, int prime, int letters_count,
    const int* parity, const std::uint64_t* offsets, const I* letters, int count,
    const std::uint64_t* starts, const unsigned char* words, int entries,
    const int* generators, std::uint32_t* output) {
    try {
        if (prime != modulus || degree > k) return 1;
        Context ctx(k,points,degree,parity,offsets,letters);
        const std::size_t coords = std::size_t(points)*ctx.masks[degree].size();
        std::vector<Polynomial> stack;
        std::string prev;
        for (int h=0; h<count; ++h) {
            std::string w(reinterpret_cast<const char*>(words+starts[h]),starts[h+1]-starts[h]);
            std::size_t c=0;
            while (c<prev.size() && c<w.size() && prev[c]==w[c]) ++c;
            stack.resize(c);
            if (stack.empty()) stack.push_back(ctx.first(static_cast<unsigned char>(w[0])));
            for (std::size_t i=stack.size(); i<w.size(); ++i)
                stack.push_back(ctx.matrix_letter(stack.back(),static_cast<unsigned char>(w[i])));
            const auto& v=stack.back();
            if (v.degree!=degree) return 2;
            for (std::size_t n=0; n<coords; ++n) {
                I acc[colour_count]={};
                const I* mat=v.values.data()+n*matrix_entries;
                for (int e=0; e<entries; ++e) {
                    const int* g=generators+3*e;
                    acc[g[0]]+=g[2]*mat[g[1]];
                }
                for (int i=0; i<colour_count; ++i)
                    output[(std::size_t(h)*coords+n)*colour_count+i]=residue(acc[i]);
            }
            prev=std::move(w);
        }
        return 0;
    } catch (const std::exception& e) {last_error=e.what();return 3;}
}

extern "C" void f4_pair(int count,int points,int nleft,int nright,int nout,
    const int* ids,const int* left_masks,const int* right_masks,const int* signs,const int* starts,
    const std::uint32_t* left,const std::uint32_t* right,I* output) {
    for (int w=0;w<count;++w) for (int p=0;p<points;++p) {
        const auto* a0=left+(std::size_t(ids[2*w])*points+p)*nleft*colour_count;
        const auto* b0=right+(std::size_t(ids[2*w+1])*points+p)*nright*colour_count;
        for (int n=0;n<nout;++n) {
            I total=0;
            for (int j=starts[n];j<starts[n+1];++j) {
                const auto* a=a0+std::size_t(left_masks[j])*colour_count;
                const auto* b=b0+std::size_t(right_masks[j])*colour_count;
                I dot=0;
                for (int i=0;i<colour_count;++i) dot+=I(a[i])*b[i];
                total+=signs[j]*dot;
            }
            output[(std::size_t(w)*points+p)*nout+n]=residue(total);
        }
    }
}
