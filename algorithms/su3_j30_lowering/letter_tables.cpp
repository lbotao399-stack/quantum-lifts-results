// Export the actual letter tables used by lower.cpp, without applying the parent projector.
#define main lowering_driver_main
#include "lower.cpp"
#undef main
int main(int argc,char** argv){try{
    require(argc==2,"letter_tables alphabet.txt"); std::ifstream f(argv[1]);std::vector<int> table;int x;
    while(f>>x)table.push_back(x);require(f.eof() && table.size()==119*3,"incomplete alphabet");setup(table);
    for(int i=1;i<=119;++i){std::cout<<i;for(const auto& d:deriv)std::cout<<' '<<d[i].first<<' '<<d[i].second;
        for(const auto& e:transformations)std::cout<<' '<<e[i].first<<' '<<e[i].second;
        std::cout<<' '<<swapped[i]<<'\n';}return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
