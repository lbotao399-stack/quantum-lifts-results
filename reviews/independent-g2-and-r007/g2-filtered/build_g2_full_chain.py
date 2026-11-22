"""Enable both dotted Taylor degrees in the public G2 classical chain executable."""
import argparse
from pathlib import Path
import subprocess
import sys


def build(repo, out):
    alg=repo.resolve()/'algorithms/dual_singlets_20261118'
    sys.path.insert(0,str(alg/'bmn_primary'))
    from build_g2_chain import build
    build(alg/'taylor_hooks/native',alg/'b2_level36/bmn_sampler.h',out)
    source=out/'g2_chain.cpp'
    data=source.read_text()
    changes={
        'int main(int argc, char **argv) {':'''int main(int argc, char **argv) {
  if (argc < 11) return 1;
  int dotted_plus = std::stoi(argv[argc-2]), dotted_minus = std::stoi(argv[argc-1]);
  if (dotted_minus < 0 || dotted_plus < dotted_minus) return 2;
  argc -= 2;''',
        'std::stoi(argv[3]), 0, 0};':'std::stoi(argv[3]), dotted_plus, dotted_minus};',
        'else build_lev(sec,curr,y);':'''else {
      fprintf(stderr,"Building length %d at dotted degrees %d,%d\\n",y,dotted_plus,dotted_minus);
      build_lev(sec,curr,y);
    }''',
        r'printf("\"weight\":[%d,%d,%d],':r'printf("\"weight\":[%d,%d,%d,%d,%d],',
        'w[0],w[1],w[2],(unsigned long long)P':'w[0],w[1],w[2],w[3],w[4],(unsigned long long)P',
        '  puts("]}");':'  puts("]}");\n  return 0;'
    }
    for old,new in changes.items():
        if data.count(old)!=1:
            raise ValueError('unexpected upstream source at '+old)
        data=data.replace(old,new)
    newfile=out/'g2_dotted_chain.cpp'
    newfile.write_text(data)
    binary=out/'g2_dotted_chain'
    subprocess.run(['c++','-std=c++17','-O3','-pthread',str(newfile),'-o',str(binary)],check=True)
    print(binary)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    build(a.repo,a.out)
