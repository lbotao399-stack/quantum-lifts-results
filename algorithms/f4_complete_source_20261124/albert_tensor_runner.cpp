// Literal finite-field contractions of the public four-leg Albert tensors.
// Every input and intermediate entry is in [0,p).  A one-leg contraction
// sums at most 26 products; a two-leg contraction at most 676 products.
// For p=1009 or 1013 both sums fit in uint32_t.  Scalar contractions use
// uint64_t and have at most 26^4 products.  No floating arithmetic is used.
#include <algorithm>
#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr std::size_t n = 26, square = n*n, fourth = square*square;
using Tensor = std::vector<std::uint32_t>;
using Matrix = std::array<std::uint32_t, square>;
constexpr int pairs[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
constexpr int orders[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
constexpr int order_signs[6] = {1,-1,-1,1,1,-1};

int sign3(int a,int b,int c) { return ((a>b)+(a>c)+(b>c))%2 ? -1 : 1; }
int sign4(int a,int b,int c,int d) {
    return ((a>b)+(a>c)+(a>d)+(b>c)+(b>d)+(c>d))%2 ? -1 : 1;
}

struct Context {
    std::uint32_t prime;
    Tensor lower, upper;
    struct Sparse {
        std::vector<std::size_t> starts;
        std::vector<std::uint32_t> columns, values;
    } exchange, exchange_transposed;
    std::array<std::vector<std::uint32_t>,6> maps;
    std::vector<std::int64_t> product;

    Context(std::uint32_t p, const std::int64_t* lo, const std::int64_t* hi,
            const std::int64_t* ex, const std::int64_t* terms, std::size_t count)
        : prime(p), lower(fourth), upper(fourth), product(terms,terms+count*8) {
        if (p!=1009 && p!=1013) throw std::runtime_error("unsupported prime");
        if (std::uint64_t(square)*(p-1)*(p-1)>UINT32_MAX)
            throw std::runtime_error("two-leg accumulator bound");
        if (std::uint64_t(fourth)*(p-1)*(p-1)>UINT64_MAX)
            throw std::runtime_error("scalar accumulator bound");
        auto residue=[p](std::int64_t value) {
            value%=p; if(value<0)value+=p; return std::uint32_t(value);
        };
        for(std::size_t i=0;i<fourth;++i) {lower[i]=residue(lo[i]);upper[i]=residue(hi[i]);}
        for(int transposed=0;transposed<2;++transposed) {
            auto& result=transposed ? exchange_transposed : exchange;
            for(std::size_t i=0;i<square;++i) {
                result.starts.push_back(result.values.size());
                for(std::size_t j=0;j<square;++j) {
                    auto value=residue(ex[transposed ? j*square+i : i*square+j]);
                    if(value) {result.columns.push_back(j);result.values.push_back(value);}
                }
            }
            result.starts.push_back(result.values.size());
        }
        for(int pair=0;pair<6;++pair) {
            int axes[4]={pairs[pair][0],pairs[pair][1],0,0},k=2;
            for(int j=0;j<4;++j)if(j!=axes[0] && j!=axes[1])axes[k++]=j;
            maps[pair].resize(fourth);
            for(std::size_t i=0;i<fourth;++i) {
                std::size_t v=i; int digits[4];
                for(int j=3;j>=0;--j){digits[axes[j]]=v%n;v/=n;}
                maps[pair][i]=((digits[0]*n+digits[1])*n+digits[2])*n+digits[3];
            }
        }
    }

    template<std::uint32_t P>
    void add(Tensor& out,const Tensor& value,int sign) const {
        if(sign>0) {
            for(std::size_t i=0;i<fourth;++i){auto v=out[i]+value[i];out[i]=v>=P?v-P:v;}
        } else {
            for(std::size_t i=0;i<fourth;++i){auto v=out[i]+P-value[i];out[i]=v>=P?v-P:v;}
        }
    }

    template<std::uint32_t P>
    Tensor apply(const Matrix& matrix,const Tensor& value,int axis) const {
        const std::size_t strides[4]={n*n*n,n*n,n,1};
        const auto suffix=strides[axis], blocks=fourth/(n*suffix);
        Tensor out(fourth,0);
        if(axis==3) {
            Matrix transposed;
            for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)transposed[j*n+i]=matrix[i*n+j];
            for(std::size_t block=0;block<blocks;++block) {
                auto* target=out.data()+block*n;
                const auto* source=value.data()+block*n;
                for(std::size_t j=0;j<n;++j) {
                    const auto multiplier=source[j];
                    const auto* row=transposed.data()+j*n;
                    #pragma GCC ivdep
                    for(std::size_t i=0;i<n;++i)target[i]+=multiplier*row[i];
                }
            }
        } else {
            for(std::size_t block=0;block<blocks;++block)for(std::size_t i=0;i<n;++i) {
                auto* target=out.data()+(block*n+i)*suffix;
                for(std::size_t j=0;j<n;++j) {
                    const auto multiplier=matrix[i*n+j];
                    if(!multiplier)continue;
                    const auto* source=value.data()+(block*n+j)*suffix;
                    #pragma GCC ivdep
                    for(std::size_t k=0;k<suffix;++k)target[k]+=multiplier*source[k];
                }
            }
        }
        for(auto& entry:out)entry%=P;
        return out;
    }

    template<std::uint32_t P>
    Tensor transfer(const Tensor& value,int pair,bool transposed) const {
        const auto& mapping=maps[pair];
        const auto& sparse=transposed ? exchange_transposed : exchange;
        Tensor reordered(fourth), acted(fourth,0), out(fourth);
        for(std::size_t i=0;i<fourth;++i)reordered[i]=value[mapping[i]];
        for(std::size_t i=0;i<square;++i) {
            auto* target=acted.data()+i*square;
            for(std::size_t k=sparse.starts[i];k<sparse.starts[i+1];++k) {
                const auto multiplier=sparse.values[k];
                const auto* source=reordered.data()+sparse.columns[k]*square;
                #pragma GCC ivdep
                for(std::size_t j=0;j<square;++j)target[j]+=multiplier*source[j];
            }
        }
        for(std::size_t i=0;i<fourth;++i)out[mapping[i]]=acted[i]%P;
        return out;
    }

    template<std::uint32_t P>
    std::uint32_t dot(const Tensor& left,const Tensor& right) const {
        std::uint64_t sum=0;
        for(std::size_t i=0;i<fourth;++i)sum+=std::uint64_t(left[i])*right[i];
        return sum%P;
    }

    template<std::uint32_t P>
    Tensor half(const std::array<Matrix,4>& scalars,const std::array<Matrix,6>& marks,
                const Tensor& corner,const std::array<int,3>& end,
                const std::array<int,3>& middle) const {
        Tensor total(fourth,0);
        for(int pair=0;pair<6;++pair) {
            const auto i=pairs[pair][0],j=pairs[pair][1];
            int other[2],number=0;for(int k=0;k<4;++k)if(k!=i&&k!=j)other[number++]=k;
            const auto k=other[0],l=other[1];
            const auto mark=apply<P>(marks[pair],corner,end[0]);
            const auto forward=apply<P>(scalars[l],apply<P>(scalars[k],mark,end[1]),end[2]);
            const auto reverse=apply<P>(scalars[k],apply<P>(scalars[l],mark,end[1]),end[2]);
            const auto sign=sign4(i,j,k,l); add<P>(total,forward,sign);add<P>(total,reverse,-sign);
        }
        Tensor out(fourth,0);
        for(int order=0;order<6;++order) {
            auto value=apply<P>(scalars[orders[order][0]],total,middle[0]);
            value=apply<P>(scalars[orders[order][1]],value,middle[1]);
            value=apply<P>(scalars[orders[order][2]],value,middle[2]);
            add<P>(out,value,order_signs[order]);
        }
        return out;
    }

    template<std::uint32_t P>
    std::array<Tensor,4> alternating(const Tensor& corner,const std::array<Matrix,4>& matrices) const {
        std::array<Tensor,4> out;
        for(auto& tensor:out)tensor.assign(fourth,0);
        for(int i=0;i<4;++i) {
            const auto first=apply<P>(matrices[i],corner,0);
            for(int j=0;j<4;++j)if(j!=i) {
                const auto second=apply<P>(matrices[j],first,1);
                for(int k=0;k<4;++k)if(k!=i&&k!=j) {
                    const auto omitted=6-i-j-k;
                    const auto third=apply<P>(matrices[k],second,2);
                    add<P>(out[omitted],third,sign3(i,j,k)*(omitted%2?-1:1));
                }
            }
        }
        return out;
    }

    template<std::uint32_t P>
    void run(const std::uint32_t* input,std::uint32_t* output) const {
        std::array<Matrix,4> x,xt;
        std::array<Matrix,6> y,yt,b,bt;
        for(std::size_t i=0;i<4;++i)for(std::size_t j=0;j<n;++j)for(std::size_t k=0;k<n;++k) {
            x[i][j*n+k]=input[i*square+j*n+k];xt[i][k*n+j]=x[i][j*n+k];
        }
        for(std::size_t pair=0;pair<6;++pair)for(std::size_t j=0;j<n;++j)for(std::size_t k=0;k<n;++k) {
            y[pair][j*n+k]=input[(4+pair)*square+j*n+k];yt[pair][k*n+j]=y[pair][j*n+k];
            std::int64_t value=0;const auto a=pairs[pair][0],c=pairs[pair][1];
            for(std::size_t l=0;l<n;++l)value+=std::int64_t(x[a][j*n+l])*x[c][l*n+k]-std::int64_t(x[c][j*n+l])*x[a][l*n+k];
            value%=P;if(value<0)value+=P;b[pair][j*n+k]=value;bt[pair][k*n+j]=value;
        }
        auto left_bracket=half<P>(xt,bt,lower,{1,0,3},{0,1,3});
        auto left_mark=half<P>(xt,yt,lower,{1,0,3},{0,1,3});
        auto right_bracket=half<P>(x,b,upper,{0,2,3},{1,2,3});
        auto right_mark=half<P>(x,y,upper,{0,2,3},{1,2,3});
        std::array<Tensor,7> a,c,d,e;
        a[0]=std::move(left_bracket);c[0]=std::move(left_mark);
        d[0]=std::move(right_bracket);e[0]=std::move(right_mark);
        for(int pair=0;pair<6;++pair) {
            a[pair+1]=transfer<P>(a[0],pair,true);c[pair+1]=transfer<P>(c[0],pair,true);
            d[pair+1]=transfer<P>(d[0],pair,false);e[pair+1]=transfer<P>(e[0],pair,false);
        }
        for(int i=0;i<7;++i)for(int j=0;j<7;++j)output[i*7+j]=(dot<P>(a[i],e[j])+P-dot<P>(c[i],d[j]))%P;
        for(auto* bank:{&a,&c,&d,&e})for(auto& tensor:*bank)Tensor().swap(tensor);
        const auto right=alternating<P>(upper,x),left=alternating<P>(lower,xt);
        std::uint32_t triples[4][4][4];
        for(int k=0;k<4;++k) {
            const auto middle=alternating<P>(right[k],x);
            for(int j=0;j<4;++j)for(int i=0;i<4;++i)triples[i][j][k]=dot<P>(left[i],middle[j]);
        }
        std::uint32_t gram[4][4],marked[6][4];
        for(int i=0;i<4;++i)for(int j=0;j<4;++j) {
            std::uint64_t sum=0;for(std::size_t k=0;k<square;++k)sum+=std::uint64_t(x[i][k])*xt[j][k];gram[i][j]=sum%P;
        }
        for(int i=0;i<6;++i)for(int j=0;j<4;++j) {
            std::uint64_t sum=0;for(std::size_t k=0;k<square;++k)sum+=std::uint64_t(y[i][k])*xt[j][k];marked[i][j]=sum%P;
        }
        std::int64_t sum=0;
        for(std::size_t i=0;i<product.size();i+=8) {
            const auto* term=product.data()+i;
            const int ids[3]={int(term[1]),int(term[2]),int(term[3])};
            std::uint32_t scalar=0;
            for(const auto& order:orders)scalar+=triples[ids[order[0]]][ids[order[1]]][ids[order[2]]];
            const auto value=std::uint64_t(scalar%P)*gram[term[4]][term[5]]%P*marked[term[6]][term[7]]%P;
            sum=(sum+term[0]*std::int64_t(value))%P;
        }
        if(sum<0)sum+=P;output[49]=sum;
    }
};
thread_local std::string error;
}

extern "C" {
void* f4_albert_create(std::uint32_t prime,const std::int64_t* lower,const std::int64_t* upper,
                      const std::int64_t* exchange,const std::int64_t* terms,std::size_t count) {
    try{return new Context(prime,lower,upper,exchange,terms,count);}
    catch(const std::exception& exception){error=exception.what();return nullptr;}
}
int f4_albert_run(void* pointer,const std::uint32_t* input,std::uint32_t* output) {
    try {const auto& context=*static_cast<Context*>(pointer);
        if(context.prime==1009)context.run<1009>(input,output);else context.run<1013>(input,output);return 0;}
    catch(const std::exception& exception){error=exception.what();return 1;}
}
const char* f4_albert_error(){return error.c_str();}
void f4_albert_destroy(void* pointer){delete static_cast<Context*>(pointer);}
}
