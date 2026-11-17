#define KCDS_SOLVE_MAX_IT 1
#define main guarded_driver_main
#include "../src/implicert.cpp"
#undef main
#include <iostream>

static void check(bool ok, const char* why) {
  if (!ok) throw std::runtime_error(why);
}
static Graph graph_of(u32 n, const std::vector<std::pair<u32,u32>>& edges) {
  std::vector<std::vector<u32>> adj(n);
  for (auto [u,v] : edges) { adj[u].push_back(v); adj[v].push_back(u); }
  for (auto& a : adj) std::sort(a.begin(), a.end());
  return csr_from_lists(adj);
}
static Graph complete(u32 n) {
  std::vector<std::pair<u32,u32>> edges;
  for (u32 i=0;i<n;++i) for(u32 j=i+1;j<n;++j) edges.push_back({i,j});
  return graph_of(n,edges);
}
static std::vector<std::array<u32,3>> triangles(const Graph& g) {
  std::vector<std::array<u32,3>> out;
  auto edge = [&](u32 a,u32 b) { return std::binary_search(g.nb(a),g.nb(a)+g.deg(a),b); };
  for(u32 a=0;a<g.n;++a) for(u32 b=a+1;b<g.n;++b) for(u32 c=b+1;c<g.n;++c)
    if(edge(a,b)&&edge(a,c)&&edge(b,c)) out.push_back({a,b,c});
  return out;
}
static Solver compact_solver(const Graph& g) {
  Solver s; s.k=3; s.G.init(g);
  s.hstamp.assign(g.n,0); s.hpos.assign(g.n,0);
  s.ts.resize(1); s.ts[0].init(g.n,3,true);
  return s;
}
struct Answer {
  std::vector<double> lam,w,r;
  std::vector<u32> witness;
  double lb=0,prep=0,iter=0; u64 count=0;
  PrResult result{PrStop::NumericalFailure,0};
  explicit Answer(u32 n):lam(n,0),w(n,1),r(n,0){}
};
static void solve(const Graph& g, Answer& a, bool compact, int budget, bool rescue=false,
                  size_t prefix=100000, double time=0, double target=1.005) {
  if(compact) {
    Solver s=compact_solver(g);
    a.result=s.global_pr_c(target,a.lam,a.w,a.r,a.lb,a.witness,a.count,budget,&a.iter,0,time,rescue,prefix);
  } else {
    a.result=global_pr(g,3,target,a.lam,0,a.w,a.r,a.lb,a.witness,budget,&a.prep,&a.iter,0,time,&a.count,rescue,prefix);
  }
}
static void verify(const Graph& g,const Answer& a) {
  auto cs=triangles(g);
  std::vector<double> loads(g.n,0); u64 count=0;
  for(auto c:cs) {
    double den=0; bool inside=true;
    for(u32 v:c) { den+=a.w[v]; inside &= std::find(a.witness.begin(),a.witness.end(),v)!=a.witness.end(); }
    for(u32 v:c) loads[v]+=a.w[v]/den;
    count+=inside;
  }
  for(u32 v=0;v<g.n;++v) {
    check(std::abs(a.w[v]-weightOf(a.lam[v],0))<1e-12*std::max(1.,a.w[v]),"weights/potentials mismatch");
    check(std::abs(loads[v]-a.r[v])<1e-11,"returned loads do not match returned potentials");
  }
  check(count==a.count,"incorrect witness count");
  double optimum=0;
  for(unsigned mask=1;mask<(1u<<g.n);++mask) {
    u64 nc=0; unsigned nv=0;
    for(u32 v=0;v<g.n;++v) nv+=(mask>>v)&1;
    for(auto c:cs) nc+=((mask>>c[0])&1)&&((mask>>c[1])&1)&&((mask>>c[2])&1);
    optimum=std::max(optimum,(double)nc/nv);
  }
  double ub=*std::max_element(a.r.begin(),a.r.end());
  check(a.lb<=optimum+1e-11 && ub+1e-11>=optimum,"bounds violate exhaustive optimum");
  if(a.result.certified()) {
    check(!a.witness.empty(),"certified result has empty witness");
    check(ub<=1.005*a.lb+1e-11,"false certified status");
  }
}
int main(int argc,char** argv) {
  if(argc>1 && std::string(argv[1])=="--driver") return guarded_driver_main(argc-1,argv+1);
  try {
    omp_set_num_threads(1);
    Graph bow=graph_of(5,{{0,1},{0,2},{1,2},{0,3},{0,4},{3,4}});
    for(bool compact:{false,true}) {
      Answer capped(bow.n); solve(bow,capped,compact,1);
      check(capped.result.stop==PrStop::IterationLimit,"expected iteration limit"); verify(bow,capped);
      Answer timed(bow.n); solve(bow,timed,compact,20,false,100000,1e-12);
      check(timed.result.stop==PrStop::TimeLimit,"expected time limit"); verify(bow,timed);
      Graph full=complete(8);
      Answer truncated(8); solve(full,truncated,compact,1,false,4);
      check(truncated.result.stop==PrStop::IterationLimit,"prefix cap must not imply success");verify(full,truncated);
      Answer rescued(8); solve(full,rescued,compact,1,true,4);
      check(rescued.result.certified() && rescued.witness.size()==8,"full-prefix rescue failed");verify(full,rescued);
      Answer missed(8); solve(full,missed,compact,3,false,1);
      check(missed.result.certified(),"zero lower bound must recover without dividing by zero");verify(full,missed);
      Graph empty=graph_of(5,{}); Answer zero(5);solve(empty,zero,compact,1);
      check(zero.result.certified()&&zero.witness.size()==1&&certificate_gap(0,0)==1,"zero-clique case");verify(empty,zero);
      Answer nan(5);nan.lam[0]=std::numeric_limits<double>::quiet_NaN();solve(bow,nan,compact,1);
      check(nan.result.stop==PrStop::NumericalFailure,"NaN was hidden by max");
    }
    check(!certificate_closed(INFINITY,1,1.01)&&!certificate_closed(0,NAN,1.01)&&
          !certificate_closed(2,1,1.01),"invalid bounds accepted");
    bool refused=false;try{require_certificate(2,1,1.01,3,"test");}catch(const std::exception&){refused=true;}
    check(refused,"output guard did not fail");
    bool badtree=false;try{MaxTree t;t.build({1,NAN});}catch(const std::exception&){badtree=true;}
    check(badtree,"tree masked a non-finite load");
     
    std::vector<std::pair<u32,u32>> possible;
    for(u32 a=0;a<5;++a)for(u32 b=a+1;b<5;++b)possible.push_back({a,b});
    for(unsigned mask=0;mask<1024;++mask){
      std::vector<std::pair<u32,u32>> edges;
      for(unsigned i=0;i<possible.size();++i)if(mask&(1u<<i))edges.push_back(possible[i]);
      Graph g=graph_of(5,edges);
      for(bool compact:{false,true}) {Answer a(5);solve(g,a,compact,150,true);verify(g,a);
        check(a.result.certified(),"small graph unexpectedly failed to close");}
    }

    const u32 q=100001;std::vector<std::pair<u32,u32>> es;
    for(u32 i=0;i<q;++i){es.push_back({i,i+q});es.push_back({i,i+2*q});es.push_back({i+q,i+2*q});}
    Graph large=graph_of(3*q,es);Answer a(large.n);solve(large,a,true,1,true);
    check(a.result.certified()&&a.count>0&&std::abs(a.lb-1./3)<1e-12,"real prefix threshold regression");
    std::cout<<"PASS: cap/time state consistency, prefix rescue, zero cliques, NaN rejection, output refusal, "
                "2048 exhaustive graph/mode checks, 300003-vertex prefix regression\n";
    return 0;
  }catch(const std::exception& e){std::cerr<<"TEST FAILED: "<<e.what()<<"\n";return 1;}
}
