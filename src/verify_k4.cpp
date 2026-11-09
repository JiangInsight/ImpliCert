#include <algorithm>
#include <array>
#include <bit>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <vector>
#include <omp.h>
#include <xmmintrin.h>
using Id=uint32_t;
using Count=uint64_t;
using Edge=std::array<Id,2>;
template<class T> T read(std::ifstream& f) { T x; f.read(reinterpret_cast<char*>(&x),sizeof x); return x; }
template<class T> std::vector<T> load(std::ifstream& f,size_t n){
  std::vector<T> v(n);f.read(reinterpret_cast<char*>(v.data()),sizeof(T)*n);return v;
}

double downsum(double a,double b){return std::max({a,b,std::nextafter(a+b,-INFINITY)});}
double upsum(double a,double b){return b==0?a:std::nextafter(a+b,INFINITY);}
void verify(const char* path){
  static_assert(std::endian::native==std::endian::little && std::numeric_limits<double>::is_iec559 && sizeof(Edge)==8);
  if(std::fegetround()!=FE_TONEAREST||(_mm_getcsr() & 0x8040)) throw std::runtime_error("requires nearest rounding and gradual underflow");
  std::ifstream f(path,std::ios::binary);f.exceptions(std::ios::failbit|std::ios::badbit);
  auto magic=load<char>(f,8);if(std::string(magic.begin(),magic.end())!="ICERT001")throw std::runtime_error("magic");
  Id n=read<Id>(f),k=read<Id>(f);Count m=read<Count>(f);Id q=read<Id>(f);Count claim=read<Count>(f);
  double raw=read<double>(f),eps=read<double>(f);
  if(k!=4||!q||q>n||!std::isfinite(raw)||raw<0||!std::isfinite(eps)||eps<=0)throw std::runtime_error("header");
  auto edges=load<Edge>(f,m);auto w=load<double>(f,n);auto witness=load<Id>(f,q);
  f.exceptions(std::ios::badbit);if(f.peek()!=EOF)throw std::runtime_error("trailing data");
  for(double x:w)if(!std::isfinite(x)||x<=0)throw std::runtime_error("weight");
  std::vector<unsigned char> selected(n);
  for(Id v:witness){if(v>=n||selected[v])throw std::runtime_error("witness");selected[v]=1;}
  std::vector<Id> degree(n);Edge prior{0,0};
  for(auto e:edges){if(e[0]>=e[1]||e[1]>=n||e<=prior)throw std::runtime_error("edge order");prior=e;++degree[e[0]];++degree[e[1]];}
  std::vector<std::vector<Id>> out(n);
  for(auto [a,b]:edges){if(std::pair{degree[a],a}>std::pair{degree[b],b})std::swap(a,b);out[a].push_back(b);}
  std::vector<Edge>().swap(edges);

  for(auto& row:out)std::sort(row.begin(),row.end());

  Count bound=0;
  for(auto& row:out){Count d=row.size();if(d<3)continue;
    if(d>2000000)throw std::runtime_error("count range");
    Count term=d*(d-1)*(d-2)/6;
    if(term>UINT64_MAX-bound)throw std::runtime_error("count range");bound+=term;}
  int T=omp_get_max_threads();
  std::vector<std::vector<double>> loads(T,std::vector<double>(n));
  Count total=0,internal=0;int invalid=0;
#pragma omp parallel reduction(+:total,internal) reduction(|:invalid)
  {
    if(std::fegetround()!=FE_TONEAREST||(_mm_getcsr() & 0x8040))invalid=1;
    auto& ub=loads[omp_get_thread_num()];std::vector<Id> common,last;
#pragma omp for schedule(dynamic,16)
    for(long long root=0;root<n;++root){
      Id a=static_cast<Id>(root);
      for(Id b:out[a]){
        common.clear();
        std::set_intersection(out[a].begin(),out[a].end(),out[b].begin(),out[b].end(),std::back_inserter(common));
        double ab=downsum(w[a],w[b]);
        for(Id c:common){
          last.clear();std::set_intersection(common.begin(),common.end(),out[c].begin(),out[c].end(),std::back_inserter(last));
          if(last.empty())continue;
          double abc=downsum(ab,w[c]);
          for(Id d:last){const double den=downsum(abc,w[d]);
            for(Id v:{a,b,c,d})ub[v]=upsum(ub[v],std::min(1.,std::nextafter(w[v]/den,INFINITY)));
            ++total;internal+=selected[a]&selected[b]&selected[c]&selected[d];}
        }
      }
    }
  }
  if(invalid||internal!=claim) throw std::runtime_error("floating environment or witness recount");
  double safe=0;for(Id v=0;v<n;++v){double x=0;for(auto& row:loads)x=upsum(x,row[v]);safe=std::max(safe,x);}
  if(!std::isfinite(safe))throw std::runtime_error("nonfinite upper endpoint");
  std::printf("checked n=%u k=%u edges=%llu cliques=%llu count=%llu size=%u raw=%.17g safe=%.17g\n",n,k,m,total,internal,q,raw,safe);
}
int main(int argc,char** argv){try{if(argc!=2)return 2;verify(argv[1]);return 0;}catch(const std::exception& e){std::fprintf(stderr,"invalid: %s\n",e.what());return 1;}}
