#include <atomic>
#include <cstdint>
#include <cstring>
#include "graph.hpp"

static Graph load_txt(const char* path) {
  FILE* f = std::fopen(path, "rb");
  if (!f) { std::fprintf(stderr, "cannot open %s\n", path); std::exit(1); }
  std::fseek(f, 0, SEEK_END);
  const long long sz = _ftelli64(f);
  std::fseek(f, 0, SEEK_SET);
  std::vector<char> buf((size_t)sz + 1);
  if (std::fread(buf.data(), 1, (size_t)sz, f) != (size_t)sz) { std::fprintf(stderr, "read error\n"); std::exit(1); }
  std::fclose(f);
  buf[(size_t)sz] = 0;
  const char* p = buf.data();
  auto next = [&]() -> u64 {
    while (*p && (*p < '0' || *p > '9')) ++p;
    u64 x = 0;
    while (*p >= '0' && *p <= '9') x = x * 10 + (u64)(*p++ - '0');
    return x;
  };
  const u32 n = (u32)next();
  const u64 m = next();
  std::vector<u32> deg(n, 0);
  std::vector<std::pair<u32, u32>> E;
  E.reserve(m);
  for (u64 i = 0; i < m; ++i) {
    const u32 a = (u32)next(), b = (u32)next();
    if (a == b || a >= n || b >= n) continue;
    E.push_back({a, b});
    deg[a]++;
    deg[b]++;
  }
  Graph g;
  g.n = n;
  g.off.assign(n + 1, 0);
  for (u32 v = 0; v < n; ++v) g.off[v + 1] = g.off[v] + deg[v];
  g.adj.resize(g.off[n]);
  std::vector<u64> pos(g.off.begin(), g.off.end() - 1);
  for (auto [a, b] : E) { g.adj[pos[a]++] = b; g.adj[pos[b]++] = a; }
  for (u32 v = 0; v < n; ++v) {
    auto* s = g.adj.data() + g.off[v];
    std::sort(s, s + deg[v]);
  }
  g.m = g.off[n] / 2;
  return g;
}

static std::atomic<bool> g_timeout{false};
static double g_limit = 60.0;
static Timer* g_timer = nullptr;

struct BB {
  int m = 0, W = 0;
  std::vector<uint64_t> adj;   
  std::atomic<int>* best = nullptr;
  u64 nodes = 0;
  static bool empty(const uint64_t* a, int W) {
    for (int i = 0; i < W; ++i) if (a[i]) return false;
    return true;
  }
   
  void expand(std::vector<uint64_t>& P, int size) {
    ++nodes;
    if ((nodes & 1023) == 0 && g_timer->sec() > g_limit) g_timeout = true;
    if (g_timeout) return;
    std::vector<int> order, color;
    std::vector<uint64_t> U = P, Q(W);
    int k = 0;
    while (!empty(U.data(), W)) {
      ++k;
      Q = U;
      for (int w = 0; w < W; ++w) {
        while (Q[w]) {
          const int b = (int)_tzcnt_u64(Q[w]);
          const int v = w * 64 + b;
          Q[w] &= Q[w] - 1;
          U[w] &= ~(1ULL << b);
          const uint64_t* av = &adj[(size_t)v * W];
          for (int x = w; x < W; ++x) Q[x] &= ~av[x];
          order.push_back(v);
          color.push_back(k);
        }
      }
    }
    std::vector<uint64_t> NP(W);
    for (int i = (int)order.size() - 1; i >= 0; --i) {
      if (size + color[i] <= best->load(std::memory_order_relaxed)) return;
      const int v = order[i];
      const uint64_t* av = &adj[(size_t)v * W];
      for (int x = 0; x < W; ++x) NP[x] = P[x] & av[x];
      if (empty(NP.data(), W)) {
        int cur = best->load();
        while (size + 1 > cur && !best->compare_exchange_weak(cur, size + 1)) {}
      } else {
        expand(NP, size + 1);
      }
      P[v >> 6] &= ~(1ULL << (v & 63));
    }
  }
};

int main(int argc, char** argv) {
  if (argc < 2) { std::fprintf(stderr, "usage: maxclique graph.{txt|bin}\n"); return 1; }
  Timer t, tall;
  g_timer = &tall;
  if (argc > 2) g_limit = std::atof(argv[2]);
  const std::string path = argv[1];
  Graph g = path.size() > 4 && path.substr(path.size() - 4) == ".bin" ? load_bin(path) : load_txt(argv[1]);
  const double t_load = t.sec();
  t.reset();
  std::vector<u32> order;
  const std::vector<u32> core = core_decomposition(g, &order);
  std::vector<u32> rank(g.n);
  for (u32 i = 0; i < g.n; ++i) rank[order[i]] = i;
  u32 degen = 0;
  for (u32 v = 0; v < g.n; ++v) degen = std::max(degen, core[v]);
  std::atomic<int> best{g.m > 0 ? 2 : (g.n > 0 ? 1 : 0)};
  std::atomic<u64> nodes{0};
#pragma omp parallel
  {
    BB bb;
    bb.best = &best;
    std::vector<u32> cand, loc(0);
    std::vector<int> where;
#pragma omp for schedule(dynamic, 1)
    for (long long i = (long long)g.n - 1; i >= 0; --i) {
      if (g_timeout) continue;
      const u32 v = order[(size_t)i];
      const int b0 = best.load();
      if ((int)core[v] + 1 <= b0) continue;
      cand.clear();
      for (u64 e = g.off[v]; e < g.off[v + 1]; ++e) {
        const u32 u = g.adj[e];
        if (rank[u] > rank[v] && (int)core[u] >= b0) cand.push_back(u);
      }
      if ((int)cand.size() + 1 <= b0) continue;
      const int m = (int)cand.size(), W = (m + 63) / 64;
      bb.m = m;
      bb.W = W;
      bb.adj.assign((size_t)m * W, 0);
       
      std::vector<u32> sc = cand;
      std::sort(sc.begin(), sc.end());
      for (int a = 0; a < m; ++a) {
        const u32 x = sc[a];
        size_t p = 0;
        for (u64 e = g.off[x]; e < g.off[x + 1] && p < sc.size(); ++e) {
          const u32 y = g.adj[e];
          while (p < sc.size() && sc[p] < y) ++p;
          if (p < sc.size() && sc[p] == y) bb.adj[(size_t)a * W + (p >> 6)] |= 1ULL << (p & 63);
        }
      }
      std::vector<uint64_t> P(W, 0);
      for (int a = 0; a < m; ++a) P[a >> 6] |= 1ULL << (a & 63);
      bb.expand(P, 1);
    }
    nodes += bb.nodes;
  }
  std::printf("graph=%s n=%u m=%llu degeneracy=%u omega=%d exact=%d load=%.2fs search=%.2fs nodes=%llu\n", argv[1], g.n,
              (unsigned long long)g.m, degen, best.load(), (int)!g_timeout, t_load, t.sec(), (unsigned long long)nodes.load());
  return 0;
}
