#pragma once
#include <condition_variable>
#include <mutex>
#include <thread>

class WorkingSetSampler {
  std::mutex lock;
  std::condition_variable wake;
  bool stopping=false;
  std::thread worker;
  size_t largest=0;
  unsigned samples=0;
  void sample() {
    PROCESS_MEMORY_COUNTERS counters{};
    if (GetProcessMemoryInfo(GetCurrentProcess(),&counters,sizeof(counters))) {
      largest=std::max(largest,static_cast<size_t>(counters.WorkingSetSize));
      ++samples;
    }
  }
public:
  WorkingSetSampler() {
    sample();
    worker=std::thread([this] {
      std::unique_lock<std::mutex> guard(lock);
      while (!wake.wait_for(guard,std::chrono::milliseconds(25),[this]{return stopping;})) sample();
    });
  }
  void stop() {
    if (!worker.joinable()) return;
    {
      std::lock_guard<std::mutex> guard(lock);
      stopping=true;
    }
    wake.notify_one();
    worker.join();
    sample();
  }
  ~WorkingSetSampler(){stop();}
  size_t peak() const {return largest;}
  unsigned count() const {return samples;}
};

struct FriendStress {
  using Edge = std::pair<u32,u32>;
  static constexpr int deletions = 5, insertions = 5;
  const double initial_ub;
  std::mt19937_64 rng;
  std::unordered_set<u32> old_witness;
  std::vector<u32> planted;
  size_t clique_size = 0;

  FriendStress(double ub, unsigned seed) : initial_ub(ub), rng(seed) {}
  static bool adjacent(const Solver& s, u32 u, u32 v) {
    const auto& nb=s.G.adj[u];
    return std::binary_search(nb.begin(),nb.end(),v);
  }
  void remember(const Solver& s) { old_witness.insert(s.S.begin(),s.S.end()); }
  static void log_edges(int b,char op,const std::vector<Edge>& edges) {
    std::printf("trace batch=%d op=%c edges=",b,op);
    for (const auto& [u,v] : edges) std::printf("%u:%u,",u,v);
    std::printf("\n");
  }
  void select_plant(const Solver& s) {
     
    size_t q=4;
    while (double(guard_audit::choose_bound(q,4))/q <= (1+s.eps)*initial_ub)
      if (++q>500) throw std::runtime_error("plant exceeds configured workload size");
    const size_t n=s.G.adj.size();
    size_t attempts=0, min_deg=n, max_deg=0;
    while (planted.size()<q) {
      if (++attempts>1000000) throw std::runtime_error("cannot select independent nonwitness plant");
      const u32 v=static_cast<u32>(rng()%n);
      if (old_witness.count(v) || s.G.adj[v].empty()) continue;
      if (std::any_of(planted.begin(),planted.end(),[&](u32 u){return u==v || adjacent(s,u,v);})) continue;
      planted.push_back(v);
      const size_t deg=s.G.adj[v].size();
      min_deg=std::min(min_deg,deg);
      max_deg=std::max(max_deg,deg);
    }
    std::printf("plant_design n=%zu q=%zu prior_witness_union=%zu initial_UB=%.17g preplant_UB=%.17g preplant_LB=%.17g min_degree=%zu max_degree=%zu selected=",n,q,old_witness.size(),initial_ub,s.tree.max(),s.rhoS(),min_deg,max_deg);
    for (const auto v:planted) std::printf("%u,",v);
    std::printf("\n");
  }
  void batch(const Solver& s,int b,std::vector<Edge>& del,std::vector<Edge>& ins) {
    if (b<=deletions) {
      remember(s);
      for (u32 u:s.S) for (u32 v:s.S) if (u<v && adjacent(s,u,v)) del.emplace_back(u,v);
      std::sort(del.begin(),del.end());
      std::shuffle(del.begin(),del.end(),rng);
      if (del.size()>32) del.resize(32);
      if (del.empty()) throw std::runtime_error("prior witness has no remaining edge");
      log_edges(b,'-',del);
    } else {
      if (planted.empty()) { remember(s); select_plant(s); }
      const int stage=b-deletions;
      const size_t next=stage==insertions ? planted.size()
        : static_cast<size_t>(planted.size()*std::sqrt(double(stage)/insertions));
      for (size_t j=clique_size;j<next;++j) for (size_t i=0;i<j;++i) {
        u32 a=planted[i], c=planted[j];
        if (adjacent(s,a,c)) throw std::runtime_error("planned insertion already present");
        ins.emplace_back(std::min(a,c),std::max(a,c));
      }
      clique_size=next;
      log_edges(b,'+',ins);
    }
    std::fflush(stdout);
  }
  void verify_plant(const Solver& s) const {
    for (size_t j=0;j<planted.size();++j)
      for (size_t i=0;i<j;++i)
        if (adjacent(s,planted[i],planted[j])!=(j<clique_size))
          throw std::runtime_error("post-update induced plant does not match K_q");
  }
  void finish(const Solver& s,int b,int attempted,double service, u64 prior_count,size_t prior_size,const WorkingSetSampler& sampler) const {
    PROCESS_MEMORY_COUNTERS pm{};
    if (!GetProcessMemoryInfo(GetCurrentProcess(),&pm,sizeof(pm))) throw std::runtime_error("process memory query failed");
    size_t overlap=0;
    for (size_t i=0;i<clique_size;++i) overlap+=s.inS[planted[i]]!=0;
    std::printf("stress_post batch=%d workload=%s local_rounds=%d prior_size=%zu prior_count=%llu returned_count=%.0f plant_size=%zu plant_count=%llu returned_plant_overlap=%zu service_s=%.9f sampled_service_peak_ws_bytes=%zu ws_samples=%u current_ws_bytes=%zu cumulative_peak_ws_bytes=%zu\n",
        b,b<=deletions?"delete":"plant",attempted,prior_size,prior_count,s.cS,clique_size,guard_audit::choose_bound(clique_size,4),overlap,service,sampler.peak(),sampler.count(),static_cast<size_t>(pm.WorkingSetSize),static_cast<size_t>(pm.PeakWorkingSetSize));
  }
};
