#define main implicert_cli
#include "../src/implicert.cpp"
#undef main

static void require(bool pass, const char* message) {
  if (!pass) throw std::runtime_error(message);
}

static void simultaneous_claims(int threads) {
  constexpr int slots = 64, epochs = 512;
  std::vector<char> flags(slots);
  std::vector<std::vector<int>> owners(threads);
  for (int epoch = 0; epoch < epochs; ++epoch) {
#pragma omp parallel
    {
      const int tid = omp_get_thread_num();
      for (int i = 0; i < slots; ++i) {
        const int slot = (i + epoch) % slots;
        if (claim_absorber(flags[slot])) owners[tid].push_back(slot);
      }
    }
    std::vector<int> frequency(slots);
    for (auto& local : owners) {
      for (int slot : local) ++frequency[slot];
      local.clear();
    }
    for (int i = 0; i < slots; ++i) {
      require(frequency[i] == 1 && flags[i] == 1, "expected one owner per touched vertex");
      flags[i] = 0;
    }
  }
}

static void hessian_and_gather(int threads) {

  constexpr u32 roots = 64, n = roots + 4;
  std::vector<std::vector<u32>> adj(n);
  for (u32 a = 0; a < n; ++a) {
    for (u32 b = std::max(a + 1, roots); b < n; ++b) {
      adj[a].push_back(b);
      adj[b].push_back(a);
    }
  }
  Solver s;
  s.G.init(csr_from_lists(adj));
  s.k = 4;
  s.eps = .01;
  s.q_abs = 4;
  s.inS.assign(n, 0);
  s.r.assign(n, 1.0);
  s.lam.assign(n, 0);
  s.w.assign(n, 1.0);
  s.xidx.assign(n, -1);
  s.vwork.assign(n, 0);
  s.absflag.assign(n, 0);
  s.absorb.assign(n, 0);
  s.abslist.resize(threads);
  s.ts.resize(threads);
  for (auto& t : s.ts) t.init(n, 4, false);
  std::vector<u32> block(roots), outside{roots, roots+1, roots+2, roots+3};
  std::iota(block.begin(), block.end(), 0);
  for (u32 v : block) s.inS[v] = 1;
  for (u32 v : outside) s.r[v] = 3.0 * roots / 4.0 + .25;
  const std::vector<double> targets(roots, 1.0);
  for (int epoch = 0; epoch < 64; ++epoch) {
    s.newton_lams(block, targets);
    std::vector<u32> claimed;
    for (const auto& local : s.abslist) claimed.insert(claimed.end(), local.begin(), local.end());
    std::sort(claimed.begin(), claimed.end());
    require(claimed == outside, "Hessian failed to record each outside vertex once");
    for (u32 v : outside) require(s.absorb[v] == 3.0*roots/4, "wrong crossing score");
    auto selected = s.gather_absorbers(2 * roots);
    std::sort(selected.begin(), selected.end());
    require(selected == outside, "wrong gathered candidates");
    for (u32 v : outside) require(s.absflag[v]==0 && s.absorb[v]==0, "gather failed to reset");
  }
}

static void clique_size_boundaries(u32 n, int k, u64 count, u64 through_pair) {
  std::vector<std::vector<u32>> adj(n);
  for (u32 v=0; v<n; ++v)
    for (u32 u=0; u<n; ++u) if (u!=v) adj[v].push_back(u);
  Solver s;
  s.G.init(csr_from_lists(adj));
  s.k=k;
  s.sat_frac=0;
  s.vwork.assign(n, 0);
  s.ts.resize(2);
  for (auto& t:s.ts) t.init(n,k);
  s.xidx.assign(n,-1);
  s.w.assign(n,1);
  s.lam.assign(n,0);
  s.inS.assign(n,1);
  std::vector<u32> block(n);
  std::iota(block.begin(),block.end(),0);
  require(s.count_inside(block)==count,"wrong exact complete-graph clique count");
  s.r.assign(n,0);
  s.dense_loads(s.w,s.r);
  for (double r:s.r) require(std::abs(r-double(count)/n)<1e-12,"wrong complete-graph loads");
  if (k<=18) {

    std::vector<double> targets=s.r;
    const double shift=.125*s.r[0];
    targets[0]-=shift;
    targets[1]+=shift;
    const double eigenvalue=double(n)*through_pair/(k*k);
    const double diagonal=eigenvalue*(n-1)/n;
    const double expected=shift/(eigenvalue+1e-6*std::max(1.0,diagonal));
    const auto next=s.newton_lams(block,targets);
    for (u32 i=0; i<n; ++i) {
      const double truth=i==0 ? expected : i==1 ? -expected : 0;
      require(std::abs(next[i]-truth)<1e-10,"wide Newton block differs from analytic direction");
    }
  }
  std::printf("passed complete_graph=%u k=%d witness_cliques=%llu\n",n,k,(unsigned long long)count);
}

static void split_singletons() {
  constexpr u32 leaves=10001, n=leaves+1;
  std::vector<std::vector<u32>> adj(n);
  for (u32 v=1; v<n; ++v) { adj[0].push_back(v); adj[v].push_back(0); }
  Solver s;
  s.G.init(csr_from_lists(adj));
  s.k=2;
  s.r.assign(n,.5); s.r[0]=leaves*.5;
  s.w.assign(n,1); s.lam.assign(n,0);
  s.xidx.assign(n,-1); s.inS.assign(n,0); s.inS[0]=1;
  s.absorb.assign(n,0); s.absflag.assign(n,0); s.abslist.resize(2);
  s.ts.resize(2);
  for (auto& t:s.ts) t.init(n,2);
  s.vwork.assign(n,0); s.vwork[0]=leaves;
  s.newton_lams({0},{s.r[0]});
  require(s.n_split>0 && s.vwork[0]==leaves,"split root missed edges");
  for (u32 v=1; v<n; ++v) require(s.absorb[v]==.5 && s.absflag[v]==1,"wrong split edge score");
  std::printf("passed split_edge_root_neighbors=%u k=2\n",leaves);
}

int main() {
  try {
    omp_set_dynamic(0);
    for (int threads : {2, 8, 24}) {
      omp_set_num_threads(threads);
      int active = 0;
#pragma omp parallel
      {
#pragma omp single
        active = omp_get_num_threads();
      }
      require(active == threads, "requested team size unavailable");
      simultaneous_claims(threads);
      hessian_and_gather(threads);
      std::printf("passed threads=%d claim_epochs=512 slots_per_epoch=64 Hessian_epochs=64 roots=64 outsiders=4\n", threads);
    }
    omp_set_num_threads(2);
    split_singletons();
    clique_size_boundaries(16,16,1,1);
    clique_size_boundaries(17,17,1,1);
    clique_size_boundaries(18,17,18,16);
    clique_size_boundaries(18,18,1,1);
    clique_size_boundaries(257,257,1,1);
    return 0;
  } catch (const std::exception& error) {
    std::fprintf(stderr, "failed: %s\n", error.what());
    return 1;
  }
}
