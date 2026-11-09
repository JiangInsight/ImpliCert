#include <array>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <random>
#include <functional>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include "cliques.hpp"
#include "absorber_claim.hpp"
#define NOMINMAX
#include <windows.h>
#include <psapi.h>

static constexpr double kBox = 300.0;
static inline double clampLam(double l, double l0) { return std::max(l0 - kBox, std::min(l0 + kBox, l)); }
static inline double weightOf(double l, double l0) {
  if (!std::isfinite(l) || !std::isfinite(l0)) return std::numeric_limits<double>::quiet_NaN();
  return std::exp(-std::max(-kBox, std::min(kBox, l - l0)));
}


static constexpr double kRoundSlack = 1e-9;


#ifndef KCDS_SOLVE_MAX_IT
#define KCDS_SOLVE_MAX_IT 100000
#endif
static constexpr int kSolveMaxIt = KCDS_SOLVE_MAX_IT, kBaselineMaxIt = 300;

enum class PrStop { Certified, IterationLimit, TimeLimit, NumericalFailure };
struct PrResult {
  PrStop stop;
  int iterations;
  bool certified() const { return stop == PrStop::Certified; }
  const char* name() const {
    switch (stop) {
      case PrStop::Certified: return "certified";
      case PrStop::IterationLimit: return "iteration_limit";
      case PrStop::TimeLimit: return "time_limit";
      default: return "numerical_failure";
    }
  }

  int log_iterations() const {
    return stop == PrStop::IterationLimit ? iterations + 1
         : (stop == PrStop::TimeLimit || stop == PrStop::NumericalFailure) ? -iterations : iterations;
  }
};

static bool certificate_closed(double ub, double lb, double factor) {
  if (!std::isfinite(ub) || !std::isfinite(lb) || !std::isfinite(factor) || ub < 0 || lb < 0 || factor < 1)
    return false;
  const double rhs = factor * lb;
  return std::isfinite(rhs) && ub <= rhs;
}
static void require_certificate(double ub, double lb, double factor, size_t witness_size, const char* phase) {
  if (witness_size && certificate_closed(ub, lb, factor)) return;
  char msg[256];
  std::snprintf(msg, sizeof msg, "uncertified phase=%s UB=%.17g LB=%.17g target=%.17g witness_size=%zu",
                phase, ub, lb, factor, witness_size);
  throw std::runtime_error(msg);
}
static double certificate_gap(double ub, double lb) {
  return lb > 0 ? ub / lb : ub == 0 ? 1.0 : std::numeric_limits<double>::infinity();
}
static void require_solve(const PrResult& result, const char* phase) {
  if (result.certified()) return;
  throw std::runtime_error(std::string("uncertified phase=") + phase + " reason=" + result.name() +
                           " iterations=" + std::to_string(result.iterations));
}


static constexpr int kStallRounds = 100;

static inline u64 ekey(u32 a, u32 b, u32 n) { return a < b ? (u64)a * n + b : (u64)b * n + a; }
static inline u64 mix64(u64 x) {
  x += 0x9E3779B97F4A7C15ull;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
  return x ^ (x >> 31);
}

struct Fenwick {
  std::vector<u64> t;
  u32 n = 0, top = 1;
  void build(const std::vector<u64>& a) {
    n = (u32)a.size();
    t.assign((size_t)n + 1, 0);
    for (u32 i = 1; i <= n; ++i) {
      t[i] += a[i - 1];
      const u32 j = i + (i & (0u - i));
      if (j <= n) t[j] += t[i];
    }
    while ((top << 1) <= n) top <<= 1;
  }
  void add(u32 i, long long d) {
    for (u32 j = i + 1; j <= n; j += j & (0u - j)) t[j] = (u64)((long long)t[j] + d);
  }
  u64 total() const {
    u64 s = 0;
    for (u32 j = n; j > 0; j -= j & (0u - j)) s += t[j];
    return s;
  }
  u32 find(u64 rank) const {
    u32 pos = 0;
    for (u32 step = top; step; step >>= 1)
      if (pos + step <= n && t[pos + step] <= rank) { pos += step; rank -= t[pos]; }
    return pos;
  }
};

struct DynGraph {
  u32 n = 0;
  std::vector<std::vector<u32>> adj;
  std::vector<std::vector<u32>> out;
  std::vector<u32> rank;
  void init(const Graph& g, bool exact_reserve = false) {
    n = g.n;
    std::vector<u32> order;
    core_decomposition(g, &order);
    rank.assign(n, 0);
    for (u32 i = 0; i < n; ++i) rank[order[i]] = i;
    adj.assign(n, {});
    out.assign(n, {});
    for (u32 v = 0; v < n; ++v) {
      adj[v].assign(g.nb(v), g.nb(v) + g.deg(v));
      if (exact_reserve) {
        size_t c = 0;
        for (u32 u : adj[v]) c += rank[u] > rank[v];
        out[v].reserve(c);
      }
      for (u32 u : adj[v]) if (rank[u] > rank[v]) out[v].push_back(u);
    }
  }
  static void ins_sorted(std::vector<u32>& a, u32 x) { a.insert(std::lower_bound(a.begin(), a.end(), x), x); }
  static void del_sorted(std::vector<u32>& a, u32 x) { a.erase(std::lower_bound(a.begin(), a.end(), x)); }


  u32 out_limit = 0xFFFFFFFFu;
  bool rerank_due = false;
  u32 max_out() const {
    u32 m = 0;
    for (const auto& o : out) m = std::max<u32>(m, (u32)o.size());
    return m;
  }
  void set_limit() {
    const u32 m = max_out();
    out_limit = std::max<u32>(2 * m, m + 16);
  }
  void rerank() {
    Graph g;
    g.n = n;
    g.off.assign((size_t)n + 1, 0);
    for (u32 v = 0; v < n; ++v) g.off[v + 1] = g.off[v] + adj[v].size();
    g.adj.resize(g.off[n]);
    for (u32 v = 0; v < n; ++v) std::copy(adj[v].begin(), adj[v].end(), g.adj.begin() + g.off[v]);
    g.m = g.off[n] / 2;
    std::vector<u32> order;
    core_decomposition(g, &order);
    for (u32 i = 0; i < n; ++i) rank[order[i]] = i;
#pragma omp parallel for schedule(dynamic, 1024)
    for (long long vv = 0; vv < (long long)n; ++vv) {
      const u32 v = (u32)vv;
      out[v].clear();
      for (u32 u : adj[v]) if (rank[u] > rank[v]) out[v].push_back(u);
    }
    set_limit();
    rerank_due = false;
  }
  void ins(u32 u, u32 v) {
    ins_sorted(adj[u], v);
    ins_sorted(adj[v], u);
    const u32 a = rank[u] < rank[v] ? u : v, b = a == u ? v : u;
    ins_sorted(out[a], b);
    if (out[a].size() > out_limit) rerank_due = true;
  }
  void del(u32 u, u32 v) {
    del_sorted(adj[u], v);
    del_sorted(adj[v], u);
    if (rank[u] < rank[v]) del_sorted(out[u], v); else del_sorted(out[v], u);
  }
};


struct LocalEnum {
  std::vector<u32> stamp, loc;
  u32 ver = 0;


  bool compact = false;
  std::vector<u32> hkey, hval;
  int hbits = 0;
  std::vector<std::vector<u32>> lout;
  std::vector<u32> glob;
  const std::vector<std::vector<u32>>* lop = nullptr;
  const std::vector<u32>* glp = nullptr;
  std::vector<unsigned char> lab;
  std::vector<std::vector<u32>> cand;
  std::vector<u32> stk, q;
  void init(u32 n, int k, bool compact_mode = false) {
    compact = compact_mode;
    if (!compact) {
      stamp.assign(n, 0);
      loc.assign(n, 0);
    }
    cand.resize(k + 2);
  }
  u32 hslot(u32 x) const { return (u32)(((u64)x * 0x9E3779B97F4A7C15ull) >> (64 - hbits)); }
  void hinsert(u32 x, u32 i) {
    const u32 mask = (1u << hbits) - 1;
    u32 p = hslot(x);
    while (hkey[p] != 0xFFFFFFFFu) p = (p + 1) & mask;
    hkey[p] = x;
    hval[p] = i;
  }
  u32 hfind(u32 x) const {
    const u32 mask = (1u << hbits) - 1;
    for (u32 p = hslot(x);; p = (p + 1) & mask) {
      if (hkey[p] == x) return hval[p];
      if (hkey[p] == 0xFFFFFFFFu) return 0xFFFFFFFFu;
    }
  }
  template <class CB>
  void rec(int level, int l, CB& cb) {
    const std::vector<u32>& Cl = cand[level];
    const std::vector<u32>& gl = *glp;
    if (l == 1) {
      for (u32 x : Cl) {
        q.clear();
        for (u32 s : stk) q.push_back(gl[s]);
        q.push_back(gl[x]);
        cb(q);
      }
      return;
    }
    const std::vector<std::vector<u32>>& lo = *lop;
    const unsigned char lv = (unsigned char)level, nl = (unsigned char)(level + 1);
    for (u32 x : Cl) {
      std::vector<u32>& nx = cand[level + 1];
      nx.clear();
      for (u32 y : lo[x]) if (lab[y] == lv) nx.push_back(y);
      if ((int)nx.size() + 1 < l) continue;
      for (u32 y : nx) lab[y] = nl;
      stk.push_back(x);
      rec(level + 1, l - 1, cb);
      stk.pop_back();
      for (u32 y : cand[level + 1]) lab[y] = lv;
    }
  }

  void build(const DynGraph& G, const std::vector<u32>& C, std::vector<std::vector<u32>>& lo, std::vector<u32>& gl) {
    const u32 m = (u32)C.size();
    gl = C;
    if (lo.size() < m) lo.resize(m);
    if (compact) {
      hbits = 4;
      while ((1u << hbits) < 2 * m + 2) ++hbits;
      if (hkey.size() < (1u << hbits)) {
        hkey.assign((size_t)1 << hbits, 0xFFFFFFFFu);
        hval.resize((size_t)1 << hbits);
      }
      for (u32 i = 0; i < m; ++i) hinsert(C[i], i);
      for (u32 i = 0; i < m; ++i) {
        std::vector<u32>& L = lo[i];
        L.clear();
        for (u32 y : G.out[C[i]]) {
          const u32 j = hfind(y);
          if (j != 0xFFFFFFFFu) L.push_back(j);
        }
      }
      for (u32 i = 0; i < m; ++i) {
        const u32 mask = (1u << hbits) - 1;
        for (u32 p = hslot(C[i]);; p = (p + 1) & mask)
          if (hkey[p] == C[i]) { hkey[p] = 0xFFFFFFFFu; break; }
      }
      return;
    }
    if (++ver == 0) { std::fill(stamp.begin(), stamp.end(), 0); ver = 1; }
    for (u32 i = 0; i < m; ++i) { stamp[C[i]] = ver; loc[C[i]] = i; }
    for (u32 i = 0; i < m; ++i) {
      std::vector<u32>& L = lo[i];
      L.clear();
      for (u32 y : G.out[C[i]]) if (stamp[y] == ver) L.push_back(loc[y]);
    }
  }

  template <class CB>
  void run_roots(const std::vector<std::vector<u32>>& lo, const std::vector<u32>& gl, int need, u32 off, u32 stride,
                 CB& cb) {
    const u32 m = (u32)gl.size();
    if (need == 1) {
      for (u32 i = off; i < m; i += stride) { q.assign(1, gl[i]); cb(q); }
      return;
    }
    lop = &lo;
    glp = &gl;
    lab.assign(m, 255);
    for (u32 i = off; i < m; i += stride) {
      if ((int)lo[i].size() + 1 < need) continue;
      cand[0] = lo[i];
      for (u32 y : cand[0]) lab[y] = 0;
      stk.assign(1, i);
      rec(0, need - 1, cb);
      for (u32 y : cand[0]) lab[y] = 255;
    }
  }
  template <class CB>
  void run(const DynGraph& G, const std::vector<u32>& C, int need, CB& cb) {
    if (need <= 0) { q.clear(); cb(q); return; }
    if (C.size() < (size_t)need) return;
    if (need == 1) {
      for (u32 x : C) { q.assign(1, x); cb(q); }
      return;
    }
    build(G, C, lout, glob);
    run_roots(lout, glob, need, 0, 1, cb);
  }
};

struct MaxTree {
  u32 sz = 1;
  std::vector<double> t;
  void build(const std::vector<double>& a) {
    while (sz < a.size()) sz <<= 1;
    t.assign(2 * sz, -1.0);
    for (size_t i = 0; i < a.size(); ++i) {
      if (!std::isfinite(a[i])) throw std::runtime_error("non-finite maintained load");
      t[sz + i] = a[i];
    }
    for (u32 i = sz - 1; i >= 1; --i) t[i] = std::max(t[2 * i], t[2 * i + 1]);
  }
  void set(u32 i, double v) {
    if (!std::isfinite(v)) throw std::runtime_error("non-finite maintained load");
    i += sz;
    t[i] = v;
    for (i >>= 1; i >= 1; i >>= 1) t[i] = std::max(t[2 * i], t[2 * i + 1]);
  }
  double max() const { return t[1]; }
  void above(double thr, std::vector<u32>& outv, u32 node = 1) const {
    if (t[node] <= thr) return;
    if (node >= sz) { outv.push_back(node - sz); return; }
    above(thr, outv, 2 * node);
    above(thr, outv, 2 * node + 1);
  }
};

struct ThreadState {
  LocalEnum le;
  std::vector<double> dr;
  std::vector<char> tflag;
  std::vector<u32> touched;
  std::vector<u32> hd;
  double dcS = 0;
  u64 cnt = 0;
  std::vector<u32> cand, q;

  std::vector<unsigned char> glab;
  std::vector<std::vector<u32>> gcand;
  std::vector<u32> gstack;
  void init(u32 n, int k, bool compact = false) {
    le.init(n, k, compact);
    dr.assign(n, 0.0);
    tflag.assign(n, 0);
    glab.assign(n, 255);
    gcand.resize(k + 1);
  }
  void add(u32 v, double d) {
    dr[v] += d;
    if (!tflag[v]) { tflag[v] = 1; touched.push_back(v); }
  }
};

struct Solver {
  int k = 4;
  double eps = 0.01, eta = 1.0;
  DynGraph G;
  std::vector<double> lam, w, r;
  double lam0 = 0;
  std::vector<char> inS;
  std::vector<u32> S;
  double cS = 0;
  MaxTree tree;
  std::vector<ThreadState> ts;
  std::vector<char> isChanged;
  std::vector<u32> changed;
  std::vector<int> xidx;
  std::vector<double> oldw;
  std::vector<u32> hpos, hstamp;
  u32 hver = 0;
  u64 cliques_touched = 0;
  double t_refresh = 0, t_relax = 0;
  int n_refresh = 0;

  double wOf(double l) const { return weightOf(l, lam0); }
  double rhoS() const { return S.empty() ? 0.0 : cS / (double)S.size(); }

  void merge() {
    for (auto& t : ts) {
      for (u32 v : t.touched) {
        r[v] += t.dr[v];
        t.dr[v] = 0;
        t.tflag[v] = 0;
        if (!isChanged[v]) { isChanged[v] = 1; changed.push_back(v); }
      }
      t.touched.clear();
      cS += t.dcS;
      t.dcS = 0;
      cliques_touched += t.cnt;
      t.cnt = 0;
    }
    for (u32 v : changed) { tree.set(v, r[v]); isChanged[v] = 0; }
    changed.clear();
  }


  std::vector<float> vwork;
  struct Task { int i, h; u32 off, stride; };
  struct SharedDAG { std::vector<std::vector<u32>> lo; std::vector<u32> gl; };
  std::vector<SharedDAG> sdag;
  std::vector<double> pass_cnt;
  bool split_heavy = std::getenv("KCDS_NOSPLIT") == nullptr;
  size_t n_split = 0;


  std::vector<Task> plan_tasks(const std::vector<u32>& X) {
    const int m = (int)X.size(), T = (int)ts.size();
    std::vector<Task> tasks;
    tasks.reserve(m + 8 * T);
    pass_cnt.assign(m, 0.0);
    std::vector<double> est(m);
    double tot = 0;
    for (int i = 0; i < m; ++i) {
      est[i] = vwork[X[i]] > 0 ? (double)vwork[X[i]] : (double)G.adj[X[i]].size();
      tot += est[i];
    }
    std::vector<int> ord(m);
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](int a, int b) { return est[a] > est[b]; });
    const double unit = tot / (4.0 * T);
    int nh = 0;
    if (split_heavy && T > 1)
      while (nh < m && est[ord[nh]] >= 2 * unit && vwork[X[ord[nh]]] > 1e4) nh++;
    if ((int)sdag.size() < nh) sdag.resize(nh);
#pragma omp parallel for schedule(dynamic, 1)
    for (int h = 0; h < nh; ++h) {
      ThreadState& t = ts[omp_get_thread_num()];
      const int i = ord[h];
      std::vector<u32> C;
      for (u32 u : G.adj[X[i]]) if (xidx[u] < 0 || xidx[u] > i) C.push_back(u);
      t.le.build(G, C, sdag[h].lo, sdag[h].gl);
    }
    for (int h = 0; h < nh; ++h) {
      const u32 c = (u32)std::max(2.0, std::min(2.0 * T, std::floor(est[ord[h]] / unit)));
      for (u32 o = 0; o < c; ++o) tasks.push_back({ord[h], h, o, c});
    }
    n_split += nh;
    for (int j = nh; j < m; ++j) tasks.push_back({ord[j], -1, 0, 1});
    return tasks;
  }

  template <class CB>
  void run_task(ThreadState& t, const Task& tk, const std::vector<u32>& X, CB& cb) {
    const u64 c0 = t.cnt;
    if (tk.h < 0) {
      t.cand.clear();
      for (u32 u : G.adj[X[tk.i]]) if (xidx[u] < 0 || xidx[u] > tk.i) t.cand.push_back(u);
      t.le.run(G, t.cand, k - 1, cb);
    } else {
      t.le.run_roots(sdag[tk.h].lo, sdag[tk.h].gl, k - 1, tk.off, tk.stride, cb);
    }
#pragma omp atomic
    pass_cnt[tk.i] += (double)(t.cnt - c0);
  }
  void record_work(const std::vector<u32>& X) {
    for (size_t i = 0; i < X.size(); ++i) vwork[X[i]] = (float)pass_cnt[i];
  }


  std::vector<char> bmark;
  void process_edges(const std::vector<std::pair<u32, u32>>& E, double sign) {
    if (E.empty()) return;
    const u32 n = G.n;
    std::unordered_map<u64, int> idx;
    idx.reserve(E.size() * 2);
    if (bmark.size() != n) bmark.assign(n, 0);
    const std::vector<char>& bv = bmark;
    for (int i = 0; i < (int)E.size(); ++i) {
      idx[ekey(E[i].first, E[i].second, n)] = i;
      bmark[E[i].first] = bmark[E[i].second] = 1;
    }
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < (int)E.size(); ++i) {
      ThreadState& t = ts[omp_get_thread_num()];
      const u32 u = E[i].first, v = E[i].second;
      t.cand.clear();
      std::set_intersection(G.adj[u].begin(), G.adj[u].end(), G.adj[v].begin(), G.adj[v].end(),
                            std::back_inserter(t.cand));
      auto cb = [&](const std::vector<u32>& rest) {
        t.q.assign({u, v});
        t.q.insert(t.q.end(), rest.begin(), rest.end());
        const auto& Q = t.q;

        for (size_t a = 0; a < Q.size(); ++a) {
          if (!bv[Q[a]]) continue;
          for (size_t b = a + 1; b < Q.size(); ++b) {
            if (!bv[Q[b]] || (a == 0 && b == 1)) continue;
            auto it = idx.find(ekey(Q[a], Q[b], n));
            if (it != idx.end() && it->second < i) return;
          }
        }
        double s = 0;
        bool allS = true;
        for (u32 x : Q) { s += w[x]; allS = allS && inS[x]; }
        for (u32 x : Q) t.add(x, sign * w[x] / s);
        if (allS) t.dcS += sign;
        t.cnt++;
        if (hdump) t.hd.insert(t.hd.end(), Q.begin(), Q.end());
      };
      t.le.run(G, t.cand, k - 2, cb);
    }
    merge();
    for (const auto& e : E) bmark[e.first] = bmark[e.second] = 0;
    if (hdump) {
      std::vector<u32>& dst = sign < 0 ? hd_del : hd_ins;
      for (auto& t : ts) { dst.insert(dst.end(), t.hd.begin(), t.hd.end()); t.hd.clear(); }
    }
  }
  bool hdump = false;
  std::vector<u32> hd_del, hd_ins;

  void apply_batch(const std::vector<std::pair<u32, u32>>& dels, const std::vector<std::pair<u32, u32>>& inss) {
    process_edges(dels, -1.0);
    for (auto [u, v] : dels) G.del(u, v);
    for (auto [u, v] : inss) G.ins(u, v);
    process_edges(inss, +1.0);
  }


  bool refreshS(double thr, size_t capH = 100000) {
    std::vector<u32> H;
    tree.above(thr, H);
    if (++hver == 0) { std::fill(hstamp.begin(), hstamp.end(), 0); hver = 1; }
    for (u32 v : H) hstamp[v] = hver;
    for (u32 v : S) if (hstamp[v] != hver) { hstamp[v] = hver; H.push_back(v); }
    std::sort(H.begin(), H.end(), [&](u32 a, u32 b) { return r[a] > r[b] || (r[a] == r[b] && a < b); });
    if (H.size() > capH) H.resize(capH);
    std::vector<u64> hist;
    prefix_hist(H, hist);
    double best = 0;
    size_t bs = 0;
    u64 acc = 0, bestc = 0;
    for (size_t i = 0; i < H.size(); ++i) {
      acc += hist[i];
      const double d = (double)acc / (double)(i + 1);
      if (d > best) { best = d; bs = i + 1; bestc = acc; }
    }
    if (best > rhoS() * (1 + 1e-12)) {
      for (u32 v : S) inS[v] = 0;
      S.assign(H.begin(), H.begin() + bs);
      for (u32 v : S) inS[v] = 1;
      cS = (double)bestc;
      return true;
    }
    return false;
  }


  void prefix_hist(const std::vector<u32>& H, std::vector<u64>& hist) {
    const u32 m = (u32)H.size();
    if (++hver == 0) { std::fill(hstamp.begin(), hstamp.end(), 0); hver = 1; }
    for (u32 i = 0; i < m; ++i) { hstamp[H[i]] = hver; hpos[H[i]] = i; }
    std::vector<std::vector<u32>> lo(m);
#pragma omp parallel for schedule(dynamic, 64)
    for (long long i = 0; i < (long long)m; ++i)
      for (u32 y : G.out[H[i]]) if (hstamp[y] == hver) lo[i].push_back(hpos[y]);
    const int T = (int)ts.size();
    std::vector<std::vector<u64>> hb(T, std::vector<u64>(m, 0));
#pragma omp parallel
    {
      std::vector<u64>& h = hb[omp_get_thread_num()];
      std::vector<unsigned char> lab(m, 255);
      std::vector<std::vector<u32>> cand(k + 1);
      std::vector<u32> stk;
      std::function<void(int, int)> rec = [&](int level, int l) {
        const std::vector<u32>& C = cand[level];
        if (l == 1) {
          u32 p0 = 0;
          for (u32 s : stk) p0 = std::max(p0, s);
          for (u32 x : C) h[std::max(p0, x)]++;
          return;
        }
        const unsigned char lv = (unsigned char)level, nl = (unsigned char)(level + 1);
        for (u32 x : C) {
          std::vector<u32>& nx = cand[level + 1];
          nx.clear();
          for (u32 y : lo[x]) if (lab[y] == lv) nx.push_back(y);
          if ((int)nx.size() + 1 < l) continue;
          for (u32 y : nx) lab[y] = nl;
          stk.push_back(x);
          rec(level + 1, l - 1);
          stk.pop_back();
          for (u32 y : cand[level + 1]) lab[y] = lv;
        }
      };
#pragma omp for schedule(dynamic, 1)
      for (long long i = 0; i < (long long)m; ++i) {
        if ((int)lo[i].size() + 1 < k) continue;
        cand[0] = lo[i];
        for (u32 y : cand[0]) lab[y] = 0;
        stk.assign(1, (u32)i);
        rec(0, k - 1);
        for (u32 y : cand[0]) lab[y] = 255;
      }
    }
    hist.assign(m, 0);
    for (auto& h : hb)
      for (u32 i = 0; i < m; ++i) hist[i] += h[i];
  }


  void dense_loads(const std::vector<double>& wv, std::vector<double>& out) {
    const u32 n = G.n;
#pragma omp parallel
    {
      ThreadState& t = ts[omp_get_thread_num()];
      std::vector<double>& buf = t.dr;
      auto leaf = [&](const std::vector<u32>& st, const std::vector<u32>& C) {
        double s0 = 0;
        for (u32 y : st) s0 += wv[y];
        double a = 0;
        for (u32 x : C) {
          const double inv = 1.0 / (s0 + wv[x]);
          a += inv;
          buf[x] += wv[x] * inv;
        }
        for (u32 y : st) buf[y] += wv[y] * a;
      };
#pragma omp for schedule(dynamic, 64)
      for (long long v = 0; v < (long long)n; ++v) {
        const u32 u = (u32)v;
        if ((int)G.out[u].size() + 1 < k) continue;
        t.gcand[0] = G.out[u];
        for (u32 y : t.gcand[0]) t.glab[y] = 0;
        t.gstack.assign(1, u);
        grec(t, 0, k - 1, leaf);
        for (u32 y : t.gcand[0]) t.glab[y] = 255;
      }
    }
#pragma omp parallel for schedule(static)
    for (long long v = 0; v < (long long)n; ++v) {
      double a = 0;
      for (auto& t : ts) { a += t.dr[v]; t.dr[v] = 0; }
      out[v] = a;
    }
  }

  double best_prefix_top(const std::vector<double>& rv, size_t capH, std::vector<u32>& Sout, u64& cnt) {
    const u32 n = G.n;
    std::vector<u32> H(n);
    std::iota(H.begin(), H.end(), 0u);
    auto cmp = [&](u32 a, u32 b) { return rv[a] > rv[b] || (rv[a] == rv[b] && a < b); };
    if (H.size() > capH) {
      std::nth_element(H.begin(), H.begin() + capH, H.end(), cmp);
      H.resize(capH);
    }
    std::sort(H.begin(), H.end(), cmp);
    std::vector<u64> hist;
    prefix_hist(H, hist);
    double best = 0;
    size_t bs = 0;
    u64 acc = 0;
    cnt = 0;
    for (size_t i = 0; i < H.size(); ++i) {
      acc += hist[i];
      const double d = (double)acc / (double)(i + 1);
      if (d > best) { best = d; bs = i + 1; cnt = acc; }
    }
    Sout.assign(H.begin(), H.begin() + bs);
    return best;
  }

  PrResult global_pr_c(double target, std::vector<double>& lamv, std::vector<double>& wv, std::vector<double>& rv,
                  double& LB, std::vector<u32>& Sout, u64& cnt, int maxIt, double* t_iter, double trigger = 0,
                  double time_cap = 0, bool rescue_prefix = false, size_t prefix_cap = 100000) {
    const u32 n = G.n;
    Timer t;
    if (maxIt < 1 || prefix_cap < 1 || n == 0 || !std::isfinite(target) || target < 1)
      throw std::invalid_argument("invalid compact PR budget, prefix cap, graph or target");
    PrResult result{PrStop::IterationLimit, maxIt};
    int it;

    std::vector<u32> Scur, bestS;
    u64 ccur = 0, best_cnt = 0;
    double best_lb = 0;
    auto keep_best = [&]() {
      Sout = bestS;
      cnt = best_cnt;
      LB = best_lb;
    };
    for (it = 1; it <= maxIt; ++it) {
#pragma omp parallel for schedule(static)
      for (long long v = 0; v < (long long)n; ++v) wv[v] = weightOf(lamv[v], lam0);
      dense_loads(wv, rv);
      double ub = 0;
      {
        std::vector<double> mx(omp_get_max_threads(), 0.0);
#pragma omp parallel
        {
          double m = 0;
#pragma omp for schedule(static) nowait
          for (long long v = 0; v < (long long)n; ++v)
            m = !std::isfinite(rv[v]) || rv[v] < 0 ? std::numeric_limits<double>::infinity() : std::max(m, rv[v]);
          mx[omp_get_thread_num()] = m;
        }
        for (double m : mx) ub = std::max(ub, m);
      }
      if (!std::isfinite(ub)) { result = {PrStop::NumericalFailure, it}; break; }
      double lb = best_prefix_top(rv, prefix_cap, Scur, ccur);


      if (prefix_cap < n && ((lb == 0 && ub > 0) ||
          (rescue_prefix && it == maxIt && !certificate_closed(ub, std::max(best_lb, lb), target)))) {
        prefix_cap = n;
        lb = best_prefix_top(rv, prefix_cap, Scur, ccur);
      }
      if (ub == 0 && bestS.empty()) bestS.assign(1, 0);
      if (lb > best_lb) {
        best_lb = lb;
        bestS = Scur;
        best_cnt = ccur;
      }
      if (it % 10 == 0)
        std::fprintf(stderr, "  global(c) it=%d UB=%.4f LB=%.4f gap=%.5f t=%.1fs\n", it, ub, best_lb, ub / best_lb,
                     t.sec());
      if (certificate_closed(ub, best_lb, target) || (it == 1 && certificate_closed(ub, best_lb, trigger))) {
        result = {PrStop::Certified, it}; break;
      }
      if (!(lb > 0) || !std::isfinite(lb)) { result = {PrStop::NumericalFailure, it}; break; }
      if (time_cap > 0 && t.sec() > time_cap) {
        keep_best();
        if (t_iter) *t_iter = t.sec();
        return {PrStop::TimeLimit, it};
      }

      if (it == maxIt) break;
#pragma omp parallel for schedule(static)
      for (long long v = 0; v < (long long)n; ++v)
        if (rv[v] > 0) lamv[v] = clampLam(lamv[v] + std::log(rv[v] / lb), lam0);
    }
    keep_best();
    if (t_iter) *t_iter = t.sec();
    return result;
  }

  u64 count_inside(const std::vector<u32>& Sv) {
    std::vector<u32> C = Sv;
    std::sort(C.begin(), C.end());
    SharedDAG d;
    ts[0].le.build(G, C, d.lo, d.gl);
    const int T = (int)ts.size();
    u64 total = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : total)
    for (int o = 0; o < T; ++o) {
      ThreadState& t = ts[omp_get_thread_num()];
      u64 c = 0;
      auto cb = [&](const std::vector<u32>&) { ++c; };
      t.le.run_roots(d.lo, d.gl, k, (u32)o, (u32)T, cb);
      total += c;
    }
    return total;
  }


  bool best_response = false;
  double damp = 1.0;
  static constexpr double kBin = 0.005, kOff = kBox + 10.0;
  static constexpr u32 kNB = (u32)(2 * kOff / kBin) + 16;
  std::vector<std::vector<double>> hcnt;
  std::vector<std::vector<u32>> hused;

  double best_lam(int tid, u32 v, double target) {
    ThreadState& t = ts[tid];
    std::vector<double>& h = hcnt[tid];
    std::vector<u32>& used = hused[tid];
    t.cand = G.adj[v];
    auto cb = [&](const std::vector<u32>& rest) {
      double s = 0;
      for (u32 x : rest) s += w[x];
      const u32 b = (u32)((std::log(s) + kOff) / kBin);
      if (h[b] == 0) used.push_back(b);
      h[b] += 1;
      t.cnt++;
    };
    t.le.run(G, t.cand, k - 1, cb);

    auto load = [&](double mu) {
      double acc = 0;
      for (u32 b : used) acc += h[b] / (1 + std::exp((b + 0.5) * kBin - kOff + mu));
      return acc;
    };
    double res = lam[v] - lam0, lo = -kBox, hi = kBox;
    if (!used.empty()) {
      if (load(hi) >= target) res = hi;
      else if (load(lo) <= target) res = lo;
      else {
        for (int it = 0; it < 60 && hi - lo > 1e-9; ++it) {
          const double mid = 0.5 * (lo + hi);
          if (load(mid) > target) lo = mid; else hi = mid;
        }
        res = 0.5 * (lo + hi);
      }
    }
    for (u32 b : used) h[b] = 0;
    used.clear();
    return lam0 + res;
  }


  double t_hess = 0, t_solve = 0;


  std::vector<double> absorb;
  std::vector<char> absflag;
  std::vector<std::vector<u32>> abslist;
  int q_abs = std::getenv("KCDS_QABS") ? std::atoi(std::getenv("KCDS_QABS")) : 256;
  size_t n_absorb = 0;
  bool abs_always = std::getenv("KCDS_ABS_TRIGGERED") == nullptr;


  std::vector<u32> gather_absorbers(double target) {
    std::vector<std::pair<double, u32>> cand;
    for (auto& l : abslist) {
      for (u32 u : l) {
        if (absflag[u]) {
          absflag[u] = 0;
          if (!inS[u] && r[u] < target * (1 - 0.5 * eps)) cand.push_back({absorb[u], u});
        }
        absorb[u] = 0;
      }
      l.clear();
    }
    std::sort(cand.begin(), cand.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
    if ((int)cand.size() > q_abs) cand.resize(q_abs);
    std::vector<u32> A;
    for (auto& c : cand) A.push_back(c.second);
    return A;
  }

  std::vector<double> newton_lams(const std::vector<u32>& X, const std::vector<double>& tv, double clip = 4.0) {
    const int m = (int)X.size();
    Timer th;
    const size_t mm = (size_t)m * m;
    const int T = (int)ts.size();

    const bool local = (size_t)T * mm <= ((size_t)32 << 20);


    constexpr int kInlineBlockSlots = 16;
    const int wide = std::min(k,m) > kInlineBlockSlots ? std::min(k,m) : 0;
    std::vector<int> wide_ids((size_t)T * wide);
    std::vector<double> wide_shares((size_t)T * wide);
    std::vector<double> H(mm, 0.0);
    std::vector<std::vector<double>> HL(local ? T : 0);
    for (auto& h : HL) h.assign(mm, 0.0);
    for (int i = 0; i < m; ++i) xidx[X[i]] = i;
    const std::vector<Task> tasks = plan_tasks(X);
#pragma omp parallel for schedule(dynamic, 1)
    for (int ti = 0; ti < (int)tasks.size(); ++ti) {
      const int tid = omp_get_thread_num();
      ThreadState& t = ts[tid];
      double* Ht = local ? HL[tid].data() : nullptr;
      const int i = tasks[ti].i;
      const u32 v = X[i];
      int small_ids[kInlineBlockSlots];
      double small_shares[kInlineBlockSlots];
      int* id = wide ? wide_ids.data() + (size_t)tid * wide : small_ids;
      double* xv = wide ? wide_shares.data() + (size_t)tid * wide : small_shares;
      auto cb = [&](const std::vector<u32>& rest) {
        double s = w[v];
        for (u32 x : rest) s += w[x];
        int c = 0;
        id[c] = i;
        xv[c++] = w[v] / s;
        for (u32 x : rest) if (xidx[x] >= 0) { id[c] = xidx[x]; xv[c++] = w[x] / s; }
        if (c < (int)rest.size() + 1) {
          double sS = inS[v] ? w[v] / s : 0.0;
          for (u32 x : rest) if (inS[x]) sS += w[x] / s;
          if (sS > 0)
            for (u32 x : rest) {
              if (xidx[x] >= 0) continue;
#pragma omp atomic
              absorb[x] += sS;
              if (claim_absorber(absflag[x])) abslist[tid].push_back(x);
            }
        }
        for (int a = 0; a < c; ++a) {
          const double dg = xv[a] * (1 - xv[a]);
          if (Ht) {
            Ht[(size_t)id[a] * m + id[a]] += dg;
            for (int b = a + 1; b < c; ++b) Ht[(size_t)id[a] * m + id[b]] -= xv[a] * xv[b];
            continue;
          }
#pragma omp atomic
          H[(size_t)id[a] * m + id[a]] += dg;
          for (int b = a + 1; b < c; ++b) {
            const double p = xv[a] * xv[b];
#pragma omp atomic
            H[(size_t)id[a] * m + id[b]] -= p;
#pragma omp atomic
            H[(size_t)id[b] * m + id[a]] -= p;
          }
        }
        t.cnt++;
      };
      run_task(t, tasks[ti], X, cb);
    }
    record_work(X);
    if (local) {

#pragma omp parallel for schedule(static)
      for (int i = 0; i < m; ++i)
        for (int j = i; j < m; ++j) {
          double acc = 0;
          for (int t = 0; t < T; ++t) acc += HL[t][(size_t)i * m + j] + (i != j ? HL[t][(size_t)j * m + i] : 0.0);
          H[(size_t)i * m + j] = H[(size_t)j * m + i] = acc;
        }
    }
    for (u32 v : X) xidx[v] = -1;
    t_hess += th.sec();
    Timer tsv;


    std::vector<int> sat;
    for (int i = 0; i < m; ++i)
      if (H[(size_t)i * m + i] < sat_frac * std::max(r[X[i]], 1.0)) sat.push_back(i);
    std::vector<char> isSat(m, 0);
    for (int i : sat) {
      isSat[i] = 1;
      for (int j = 0; j < m; ++j) H[(size_t)i * m + j] = H[(size_t)j * m + i] = 0.0;
      H[(size_t)i * m + i] = 1.0;
    }
    double dmax = 0;
    for (int i = 0; i < m; ++i) if (!isSat[i]) dmax = std::max(dmax, H[(size_t)i * m + i]);


    std::vector<double> diag(m), b(m), x(m), res(m), z(m), p(m), Ap(m), dinv(m);
    for (int i = 0; i < m; ++i) { diag[i] = H[(size_t)i * m + i]; b[i] = isSat[i] ? 0.0 : r[X[i]] - tv[i]; }
    double delta = 1e-6 * std::max(dmax, 1.0);
    auto matvec = [&](const std::vector<double>& in, std::vector<double>& out) {
#pragma omp parallel for schedule(static)
      for (int i = 0; i < m; ++i) {
        const double* row = &H[(size_t)i * m];
        double acc = delta * in[i];
        for (int j = 0; j < m; ++j) acc += row[j] * in[j];
        out[i] = acc;
      }
    };
    double bn = 0;
    for (int i = 0; i < m; ++i) bn += b[i] * b[i];
    int attempt = 0, cg_it = 0;
    double xm = 0;
    for (attempt = 0; attempt < 12; ++attempt) {
      std::fill(x.begin(), x.end(), 0.0);
      for (int i = 0; i < m; ++i) dinv[i] = 1.0 / (diag[i] + delta);
      res = b;
      for (int i = 0; i < m; ++i) { z[i] = dinv[i] * res[i]; p[i] = z[i]; }
      double rz = 0;
      for (int i = 0; i < m; ++i) rz += res[i] * z[i];
      for (int it = 0; it < 200 && bn > 0; ++it, ++cg_it) {
        matvec(p, Ap);
        double pAp = 0;
        for (int i = 0; i < m; ++i) pAp += p[i] * Ap[i];
        if (pAp <= 0) break;
        const double alpha = rz / pAp;
        double rn = 0;
        for (int i = 0; i < m; ++i) { x[i] += alpha * p[i]; res[i] -= alpha * Ap[i]; rn += res[i] * res[i]; }
        if (rn <= 1e-8 * bn) break;
        double rz2 = 0;
        for (int i = 0; i < m; ++i) { z[i] = dinv[i] * res[i]; rz2 += res[i] * z[i]; }
        for (int i = 0; i < m; ++i) p[i] = z[i] + (rz2 / rz) * p[i];
        rz = rz2;
      }
      xm = 0;
      for (int i = 0; i < m; ++i) xm = std::max(xm, std::fabs(x[i]));
      if (xm <= clip) break;
      delta *= 10;
    }
    if (verbose) {
      int nlo = 0, nhi = 0;
      for (int i = 0; i < m; ++i) {
        const double mu = lam[X[i]] - lam0;
        nlo += mu <= -kBox + 1e-9;
        nhi += mu >= kBox - 1e-9;
      }
      std::fprintf(stderr,
                   "      newton m=%d sat=%zu attempts=%d delta/dmax=%.1e cg_it=%d max|d|=%.3f |b|=%.1f box_lo=%d box_hi=%d\n",
                   m, sat.size(), attempt, delta / std::max(dmax, 1.0), cg_it, xm, std::sqrt(bn), nlo, nhi);
    }
    std::vector<double> nl(m);
    for (int i = 0; i < m; ++i) nl[i] = lam[X[i]] + newton_damp * std::max(-clip, std::min(clip, x[i]));
#pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < (int)sat.size(); ++j) {
      const u32 v = X[sat[j]];
      nl[sat[j]] = lam[v] + sat_damp * (best_lam(omp_get_thread_num(), v, tv[sat[j]]) - lam[v]);
    }
    t_solve += tsv.sec();
    return nl;
  }

  double sat_frac = std::getenv("KCDS_SAT") ? std::atof(std::getenv("KCDS_SAT")) : 0.1;
  bool uniform_targets = std::getenv("KCDS_UNIFORM") != nullptr;
  double newton_damp = std::getenv("KCDS_NDAMP") ? std::atof(std::getenv("KCDS_NDAMP")) : 1.0;
  double sat_damp = std::getenv("KCDS_SDAMP") ? std::atof(std::getenv("KCDS_SDAMP")) : 1.0;

  void relax(const std::vector<u32>& X, double target, const std::vector<double>* given = nullptr) {
    oldw.resize(X.size());
    std::vector<double> nl(X.size());
    if (given) nl = *given;
    else if (best_response) {
#pragma omp parallel for schedule(dynamic, 1)
      for (int i = 0; i < (int)X.size(); ++i) nl[i] = best_lam(omp_get_thread_num(), X[i], target);
    }
    for (size_t i = 0; i < X.size(); ++i) {
      const u32 v = X[i];
      oldw[i] = w[v];
      xidx[v] = (int)i;
      const double step = given ? nl[i] - lam[v]
                          : best_response ? damp * (nl[i] - lam[v])
                                          : eta * std::log(std::max(r[v], 1e-300) / target);
      lam[v] = clampLam(lam[v] + step, lam0);
      w[v] = wOf(lam[v]);
    }
    const std::vector<Task> tasks = plan_tasks(X);
#pragma omp parallel for schedule(dynamic, 1)
    for (int ti = 0; ti < (int)tasks.size(); ++ti) {
      ThreadState& t = ts[omp_get_thread_num()];
      const int i = tasks[ti].i;
      const u32 v = X[i];
      auto cb = [&](const std::vector<u32>& rest) {
        double so = oldw[i], sn = w[v];
        for (u32 x : rest) { sn += w[x]; so += xidx[x] >= 0 ? oldw[xidx[x]] : w[x]; }
        t.add(v, w[v] / sn - oldw[i] / so);
        for (u32 x : rest) t.add(x, w[x] / sn - (xidx[x] >= 0 ? oldw[xidx[x]] : w[x]) / so);
        t.cnt++;
      };
      run_task(t, tasks[ti], X, cb);
    }
    record_work(X);
    for (u32 v : X) xidx[v] = -1;
    merge();
  }


  template <class Leaf>
  void grec(ThreadState& t, int level, int l, Leaf& leaf) {
    const std::vector<u32>& C = t.gcand[level];
    if (l == 1) { leaf(t.gstack, C); return; }
    const unsigned char lv = (unsigned char)level, nl = (unsigned char)(level + 1);
    for (u32 x : C) {
      std::vector<u32>& nx = t.gcand[level + 1];
      nx.clear();
      for (u32 y : G.out[x]) if (t.glab[y] == lv) nx.push_back(y);
      if ((int)nx.size() + 1 < l) continue;
      for (u32 y : nx) t.glab[y] = nl;
      t.gstack.push_back(x);
      grec(t, level + 1, l - 1, leaf);
      t.gstack.pop_back();
      for (u32 y : t.gcand[level + 1]) t.glab[y] = lv;
    }
  }

  std::vector<char> inR;
  std::vector<u32> rstamp;
  u32 rver = 0;


  void relax2(const std::vector<u32>& R, double target) {
    oldw.resize(R.size());
    for (size_t i = 0; i < R.size(); ++i) {
      const u32 v = R[i];
      oldw[i] = w[v];
      xidx[v] = (int)i;
      inR[v] = 1;
      const double step = std::max(-2.0, std::min(2.0, std::log(std::max(r[v], 1e-300) / target)));
      lam[v] = clampLam(lam[v] + eta * step, lam0);
      w[v] = wOf(lam[v]);
    }
    if (++rver == 0) { std::fill(rstamp.begin(), rstamp.end(), 0); rver = 1; }
    std::vector<u32> roots;
    for (u32 v : R) {
      if (rstamp[v] != rver) { rstamp[v] = rver; roots.push_back(v); }
      for (u32 u : G.adj[v])
        if (G.rank[u] < G.rank[v] && rstamp[u] != rver) { rstamp[u] = rver; roots.push_back(u); }
    }
#pragma omp parallel for schedule(dynamic, 16)
    for (long long ri = 0; ri < (long long)roots.size(); ++ri) {
      ThreadState& t = ts[omp_get_thread_num()];
      const u32 m = roots[ri];
      if ((int)G.out[m].size() + 1 < k) continue;
      t.gcand[0] = G.out[m];
      for (u32 y : t.gcand[0]) t.glab[y] = 0;
      t.gstack.assign(1, m);
      auto leaf = [&](const std::vector<u32>& st, const std::vector<u32>& C) {
        bool stR = false;
        double so0 = 0, sn0 = 0;
        for (u32 y : st) {
          stR = stR || inR[y];
          sn0 += w[y];
          so0 += xidx[y] >= 0 ? oldw[xidx[y]] : w[y];
        }
        for (u32 x : C) {
          if (!stR && !inR[x]) continue;
          const double wo_x = xidx[x] >= 0 ? oldw[xidx[x]] : w[x];
          const double sn = sn0 + w[x], so = so0 + wo_x;
          for (u32 y : st) {
            if (!inR[y] && sn == so) continue;
            t.add(y, w[y] / sn - (xidx[y] >= 0 ? oldw[xidx[y]] : w[y]) / so);
          }
          t.add(x, w[x] / sn - wo_x / so);
          t.cnt++;
        }
      };
      grec(t, 0, k - 1, leaf);
      for (u32 y : t.gcand[0]) t.glab[y] = 255;
    }
    for (u32 v : R) { xidx[v] = -1; inR[v] = 0; }
    merge();
  }

  bool two_sided = true;
  int strategy = 2;
  int q = 64;
  bool verbose = std::getenv("KCDS_VERBOSE") != nullptr;
  bool use_newton = false, hybrid = false;
  int newton_max = 3000;
  u64 newton_cost = 0;
  int n_newton = 0;


  double budget_s = 0;


  int repair(int maxRounds = 300) {
    double lb = rhoS();
    if (tree.max() <= (1 + eps) * (1 - kRoundSlack) * lb) return 0;
    Timer tbud;
    const double goal = 1 + 0.5 * eps;
    u64 spent_targeted = 0;
    int streak = 0;
    for (int round = 1; round <= maxRounds; ++round) {
      if (budget_s > 0 && round > 1 && tbud.sec() > budget_s) return -2;
      if (round == 1 || round % 10 == 0) {
        Timer tt;
        if (refreshS(lb)) lb = rhoS();
        t_refresh += tt.sec();
        n_refresh++;
      }
      if (tree.max() <= goal * lb) return round;
      Timer trx;
      std::vector<u32> X;
      tree.above(goal * lb, X);
      if (strategy == 2) {

        std::sort(X.begin(), X.end(), [&](u32 a, u32 b) { return r[a] > r[b]; });
        if ((int)X.size() > q) X.resize(q);
        std::vector<u32> under;
        for (u32 v : S) if (r[v] < lb * (1 - 0.5 * eps)) under.push_back(v);
        std::sort(under.begin(), under.end(), [&](u32 a, u32 b) { return r[a] < r[b]; });
        if (verbose) {
          const u32 top = X.empty() ? 0 : X[0];
          size_t nover = 0, overS = 0;
          {
            std::vector<u32> all;
            tree.above(goal * lb, all);
            nover = all.size();
            for (u32 v : all) overS += inS[v];
          }
          double sumS = 0, muS = 0;
          for (u32 v : S) { sumS += r[v]; muS += lam[v] - lam0; }
          std::fprintf(stderr,
                       "    round=%d UB=%.2f LB=%.2f over=%zu (inS %zu) under=%zu |S|=%zu top: r=%.2f inS=%d mu=%.2f"
                       " | avgS/LB=%.5f crossS/cS=%.5f meanmuS=%.2f\n",
                       round, tree.max(), lb, nover, overS, under.size(), S.size(), r[top], (int)inS[top],
                       lam[top] - lam0, sumS / S.size() / lb, (sumS - cS) / cS, muS / S.size());
        }
        if ((int)under.size() > q) under.resize(q);
        std::vector<u32> NX;


        const bool newton_now =
            use_newton && (!hybrid || (newton_cost == 0 ? round > 5
                                                        : 2 * spent_targeted >= newton_cost || streak >= kStallRounds));
        if (newton_now) {

          tree.above(goal * lb, NX);
          if (++hver == 0) { std::fill(hstamp.begin(), hstamp.end(), 0); hver = 1; }
          for (u32 v : NX) hstamp[v] = hver;
          for (u32 v : S) if (hstamp[v] != hver) { hstamp[v] = hver; NX.push_back(v); }
          if ((int)NX.size() > newton_max) NX.clear();
        }
        const u64 c0 = cliques_touched;
        if (!NX.empty()) {


          const double tgt = lb * (1 + 0.25 * eps);
          double sumS = 0;
          for (u32 v : S) sumS += r[v];
          const double tS = uniform_targets ? tgt : std::max(tgt, sumS / std::max<size_t>(S.size(), 1));
          std::vector<double> tv(NX.size());
          for (size_t i = 0; i < NX.size(); ++i) tv[i] = inS[NX[i]] ? tS : tgt;
          const std::vector<double> nl = newton_lams(NX, tv);
          relax(NX, tgt, &nl);


          const std::vector<u32> A = gather_absorbers(tgt);
          if ((tS > tgt || abs_always) && !A.empty()) { relax(A, tgt); n_absorb += A.size(); }
          newton_cost = cliques_touched - c0;
          n_newton++;
          streak = 0;
        } else if (best_response) {


          if (!X.empty()) relax(X, lb);
          under.clear();
          for (u32 v : S) if (r[v] < lb * (1 - 0.5 * eps)) under.push_back(v);
          std::sort(under.begin(), under.end(), [&](u32 a, u32 b) { return r[a] < r[b]; });
          if ((int)under.size() > q) under.resize(q);
          if (!under.empty()) relax(under, lb);
        } else {
          X.insert(X.end(), under.begin(), under.end());
          relax(X, lb);
          spent_targeted += cliques_touched - c0;
          ++streak;
        }
      } else if (two_sided) {

        if (++hver == 0) { std::fill(hstamp.begin(), hstamp.end(), 0); hver = 1; }
        for (u32 v : X) hstamp[v] = hver;
        for (u32 v : S) if (hstamp[v] != hver) { hstamp[v] = hver; X.push_back(v); }
        relax2(X, lb);
      } else {
        relax(X, lb);
      }
      t_relax += trx.sec();
      lb = rhoS();
      if (tree.max() <= goal * lb) return round;
    }
    return -1;
  }
};


struct Stream {
  u32 n = 0;
  u64 M = 0;
  std::vector<u32> u, v;
  std::vector<int8_t> op;
};


static Stream load_events(const char* path) {
  Stream s;
  FILE* f = std::fopen(path, "rb");
  u64 h[2];
  if (!f || std::fread(h, 8, 2, f) != 2) { std::fprintf(stderr, "cannot read %s\n", path); std::exit(1); }
  s.n = (u32)h[0];
  s.M = h[1];
  s.u.resize(s.M);
  s.v.resize(s.M);
  s.op.resize(s.M);
  if (std::fread(s.u.data(), 4, s.M, f) != s.M || std::fread(s.v.data(), 4, s.M, f) != s.M ||
      std::fread(s.op.data(), 1, s.M, f) != s.M) {
    std::fprintf(stderr, "truncated %s\n", path);
    std::exit(1);
  }
  std::fclose(f);
  return s;
}
static Stream load_stream(const char* path) {
  Stream s;
  FILE* f = std::fopen(path, "rb");
  u64 h[2];
  if (!f || std::fread(h, 8, 2, f) != 2) { std::fprintf(stderr, "cannot read %s\n", path); std::exit(1); }
  s.n = (u32)h[0];
  s.M = h[1];
  s.u.resize(s.M);
  s.v.resize(s.M);
  if (std::fread(s.u.data(), 4, s.M, f) != s.M || std::fread(s.v.data(), 4, s.M, f) != s.M) {
    std::fprintf(stderr, "truncated %s\n", path);
    std::exit(1);
  }
  std::fclose(f);
  return s;
}


static constexpr u32 kNone = 0xFFFFFFFFu;
static constexpr size_t kPrefixCap = 100000;
struct GTask { u32 root, off, stride; };


template <class Leaf>
static void run_gtask(const DAG& d, int k, const GTask& t, std::vector<unsigned char>& lab,
                      std::vector<std::vector<u32>>& cand, std::vector<u32>& stack, Leaf& leaf) {
  const u32 u = t.root, du = d.outdeg(u);
  stack.assign(1, u);
  if (k == 2) {
    std::vector<u32>& sub = cand[1];
    sub.clear();
    for (u32 i = t.off; i < du; i += t.stride) sub.push_back(d.out(u)[i]);
    leaf(stack, sub);
    return;
  }
  cand[0].assign(d.out(u), d.out(u) + du);
  for (u32 y : cand[0]) lab[y] = 0;
  for (u32 i = t.off; i < du; i += t.stride) {
    const u32 x = cand[0][i];
    std::vector<u32>& next = cand[1];
    next.clear();
    for (u64 e = d.off[x]; e < d.off[x + 1]; ++e) if (lab[d.adj[e]] == 0) next.push_back(d.adj[e]);
    if (next.size() + 2 < (size_t)k) continue;
    for (u32 y : next) lab[y] = 1;
    stack.push_back(x);
    kc_rec(d, lab, cand, stack, 1, k - 2, leaf);
    stack.pop_back();
    for (u32 y : cand[1]) lab[y] = 0;
  }
  for (u32 y : cand[0]) lab[y] = 255;
}


static void gpr_loads(const DAG& d, int k, const std::vector<GTask>& tasks, size_t nhead, const std::vector<double>& w,
                      std::vector<std::vector<double>>& rb, std::vector<double>& r, std::vector<u64>* work) {
#pragma omp parallel
  {
    std::vector<double>& buf = rb[omp_get_thread_num()];
    u64 cnt = 0;
    auto leaf = [&](const std::vector<u32>& st, const std::vector<u32>& C) {
      double s0 = 0;
      for (u32 y : st) s0 += w[y];
      double a = 0;
      for (u32 x : C) {
        const double inv = 1.0 / (s0 + w[x]);
        a += inv;
        buf[x] += w[x] * inv;
      }
      for (u32 y : st) buf[y] += w[y] * a;
      cnt += C.size();
    };
    std::vector<unsigned char> lab(d.n, 255);
    std::vector<std::vector<u32>> cand(k + 1);
    std::vector<u32> stack;
    auto run = [&](long long i) {
      const u64 c0 = cnt;
      run_gtask(d, k, tasks[i], lab, cand, stack, leaf);
      if (work) (*work)[i] = cnt - c0;
    };
#pragma omp for schedule(dynamic, 1) nowait
    for (long long i = 0; i < (long long)nhead; ++i) run(i);
#pragma omp for schedule(dynamic, 64)
    for (long long i = (long long)nhead; i < (long long)tasks.size(); ++i) run(i);
  }
  reduce_into(rb, r);
}


static void gpr_split(const DAG& d, std::vector<GTask>& tasks, const std::vector<u64>& work, int T, size_t& nhead) {
  u64 total = 0;
  for (u64 x : work) total += x;
  const u64 lim = std::max<u64>(1, total / (4 * (u64)T)), piece = std::max<u64>(1, total / (8 * (u64)T));
  std::vector<std::pair<u64, GTask>> out;
  out.reserve(tasks.size());
  for (size_t i = 0; i < tasks.size(); ++i) {
    const GTask& t = tasks[i];
    if (work[i] == 0) continue;
    const u32 du = d.outdeg(t.root);
    if (work[i] > lim && t.stride == 1 && du > 1) {
      const u32 parts = (u32)std::min<u64>(du, (work[i] + piece - 1) / piece);
      for (u32 p = 0; p < parts; ++p) out.push_back({work[i] / parts, GTask{t.root, p, parts}});
    } else {
      out.push_back({work[i], t});
    }
  }
  std::stable_sort(out.begin(), out.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
  tasks.clear();
  for (auto& p : out) tasks.push_back(p.second);
  nhead = std::min<size_t>(tasks.size(), 64 * (size_t)T);
}


template <class Leaf>
static void lo_rec(const std::vector<std::vector<u32>>& lo, std::vector<unsigned char>& lab,
                   std::vector<std::vector<u32>>& cand, std::vector<u32>& stack, int level, int l, Leaf& leaf) {
  const std::vector<u32>& C = cand[level];
  if (l == 1) { leaf(stack, C); return; }
  const unsigned char lv = (unsigned char)level, nl = (unsigned char)(level + 1);
  for (u32 x : C) {
    std::vector<u32>& next = cand[level + 1];
    next.clear();
    for (u32 y : lo[x]) if (lab[y] == lv) next.push_back(y);
    if (next.size() + 1 < (size_t)l) continue;
    for (u32 y : next) lab[y] = nl;
    stack.push_back(x);
    lo_rec(lo, lab, cand, stack, level + 1, l - 1, leaf);
    stack.pop_back();
    for (u32 y : cand[level + 1]) lab[y] = lv;
  }
}


static double gpr_prefix(const DAG& d, int k, const std::vector<double>& r, std::vector<u32>& H, std::vector<u32>& pos,
                         u32& sz, u64& cnt, size_t prefix_cap = kPrefixCap) {
  const u32 n = d.n;
  const int T = omp_get_max_threads();
  std::vector<std::vector<u32>> part(T);
#pragma omp parallel
  {
    std::vector<u32>& p = part[omp_get_thread_num()];
#pragma omp for schedule(static)
    for (long long v = 0; v < (long long)n; ++v) if (r[v] > 0) p.push_back((u32)v);
  }
  H.clear();
  for (auto& p : part) H.insert(H.end(), p.begin(), p.end());
  auto cmp = [&](u32 a, u32 b) { return r[a] > r[b] || (r[a] == r[b] && a < b); };
  if (H.size() > prefix_cap) {
    std::nth_element(H.begin(), H.begin() + prefix_cap, H.end(), cmp);
    H.resize(prefix_cap);
  }
  std::sort(H.begin(), H.end(), cmp);
  const u32 m = (u32)H.size();
  for (u32 i = 0; i < m; ++i) pos[H[i]] = i;
  std::vector<std::vector<u32>> lo(m);
#pragma omp parallel for schedule(dynamic, 256)
  for (long long i = 0; i < (long long)m; ++i)
    for (u64 e = d.off[H[i]]; e < d.off[H[i] + 1]; ++e) {
      const u32 p = pos[d.adj[e]];
      if (p != kNone) lo[i].push_back(p);
    }
  for (u32 i = 0; i < m; ++i) pos[H[i]] = kNone;
  std::vector<u32> roots;
  for (u32 i = 0; i < m; ++i) if (lo[i].size() + 1 >= (size_t)k) roots.push_back(i);
  std::stable_sort(roots.begin(), roots.end(), [&](u32 a, u32 b) { return lo[a].size() > lo[b].size(); });
  const size_t nhead = std::min<size_t>(roots.size(), 64 * (size_t)T);
  std::vector<std::vector<u64>> hb(T, std::vector<u64>(m, 0));
#pragma omp parallel
  {
    std::vector<u64>& h = hb[omp_get_thread_num()];
    auto leaf = [&](const std::vector<u32>& st, const std::vector<u32>& C) {
      u32 p0 = 0;
      for (u32 s : st) p0 = std::max(p0, s);
      for (u32 x : C) h[std::max(p0, x)]++;
    };
    std::vector<unsigned char> lab(m, 255);
    std::vector<std::vector<u32>> cand(k + 1);
    std::vector<u32> stack;
    auto run = [&](long long j) {
      const u32 i = roots[j];
      cand[0] = lo[i];
      for (u32 y : cand[0]) lab[y] = 0;
      stack.assign(1, i);
      lo_rec(lo, lab, cand, stack, 0, k - 1, leaf);
      for (u32 y : cand[0]) lab[y] = 255;
    };
#pragma omp for schedule(dynamic, 1) nowait
    for (long long j = 0; j < (long long)nhead; ++j) run(j);
#pragma omp for schedule(dynamic, 16)
    for (long long j = (long long)nhead; j < (long long)roots.size(); ++j) run(j);
  }
  std::vector<u64> hist(m, 0);
#pragma omp parallel for schedule(static)
  for (long long i = 0; i < (long long)m; ++i) {
    u64 s = 0;
    for (auto& h : hb) s += h[i];
    hist[i] = s;
  }
  double best = 0;
  u64 acc = 0;
  sz = 0;
  cnt = 0;
  for (u32 i = 0; i < m; ++i) {
    acc += hist[i];
    const double dens = (double)acc / (double)(i + 1);
    if (dens > best) { best = dens; sz = i + 1; cnt = acc; }
  }
  return best;
}


static PrResult global_pr(const Graph& g, int k, double target, std::vector<double>& lam, double lam0,
                     std::vector<double>& w, std::vector<double>& r, double& LB, std::vector<u32>& Sout, int maxIt,
                     double* t_prep, double* t_iter, double trigger = 0, double time_cap = 0, u64* cnt_out = nullptr,
                     bool rescue_prefix = false, size_t prefix_cap = kPrefixCap) {
  if (maxIt < 1 || prefix_cap < 1 || g.n == 0 || !std::isfinite(target) || target < 1)
    throw std::invalid_argument("invalid PR budget, prefix cap, graph or target");
  PrResult result{PrStop::IterationLimit, maxIt};
  Timer ttot;
  Timer t;
  const int T = omp_get_max_threads();
  std::vector<u32> order;
  const std::vector<u32> core = core_decomposition(g, &order);
  std::vector<u32> nid(g.n, kNone), old;
  for (u32 v : order) if (core[v] + 1 >= (u32)k) { nid[v] = (u32)old.size(); old.push_back(v); }
  const u32 n2 = (u32)old.size();
  DAG d;
  d.n = n2;
  d.off.assign((size_t)n2 + 1, 0);
#pragma omp parallel for schedule(dynamic, 1024)
  for (long long a = 0; a < (long long)n2; ++a) {
    const u32 v = old[a];
    u64 c = 0;
    for (u64 e = g.off[v]; e < g.off[v + 1]; ++e) {
      const u32 b = nid[g.adj[e]];
      if (b != kNone && b > (u32)a) ++c;
    }
    d.off[a + 1] = c;
  }
  for (u32 a = 0; a < n2; ++a) d.off[a + 1] += d.off[a];
  d.adj.resize(d.off[n2]);
#pragma omp parallel for schedule(dynamic, 1024)
  for (long long a = 0; a < (long long)n2; ++a) {
    const u32 v = old[a];
    u64 p = d.off[a];
    for (u64 e = g.off[v]; e < g.off[v + 1]; ++e) {
      const u32 b = nid[g.adj[e]];
      if (b != kNone && b > (u32)a) d.adj[p++] = b;
    }
  }
  std::vector<GTask> tasks;
  for (u32 a = 0; a < n2; ++a) if (d.outdeg(a) + 1 >= (u32)k) tasks.push_back(GTask{a, 0, 1});
  std::stable_sort(tasks.begin(), tasks.end(),
                   [&](const GTask& x, const GTask& y) { return d.outdeg(x.root) > d.outdeg(y.root); });
  size_t nhead = std::min<size_t>(tasks.size(), 64 * (size_t)T);
  std::vector<u64> work(tasks.size(), 0);
  std::vector<std::vector<double>> rb(T, std::vector<double>(n2, 0.0));
  std::vector<double> wl(n2), rl(n2, 0.0);
  std::vector<u32> pos(n2, kNone), H;
  if (t_prep) *t_prep = t.sec();
  t.reset();
  double lb = 0, best_lb = 0;
  u32 sz = 0;
  u64 cnt = 0, best_cnt = 0;
  std::vector<u32> bestS;
  auto finish = [&]() {
#pragma omp parallel for schedule(static)
    for (long long v = 0; v < (long long)g.n; ++v) { r[v] = 0; w[v] = weightOf(lam[v], lam0); }
#pragma omp parallel for schedule(static)
    for (long long a = 0; a < (long long)n2; ++a) r[old[a]] = rl[a];
    Sout = bestS;
    LB = best_lb;
    if (cnt_out) *cnt_out = best_cnt;
    if (t_iter) *t_iter = t.sec();
  };
  int it;
  for (it = 1; it <= maxIt; ++it) {
#pragma omp parallel for schedule(static)
    for (long long a = 0; a < (long long)n2; ++a) wl[a] = weightOf(lam[old[a]], lam0);
    gpr_loads(d, k, tasks, nhead, wl, rb, rl, it == 1 ? &work : nullptr);
    double ub = 0;
    {
      std::vector<double> mx(T, 0.0);
#pragma omp parallel
      {
        double m = 0;
#pragma omp for schedule(static) nowait
        for (long long a = 0; a < (long long)n2; ++a)
          m = !std::isfinite(rl[a]) || rl[a] < 0 ? std::numeric_limits<double>::infinity() : std::max(m, rl[a]);
        mx[omp_get_thread_num()] = m;
      }
      for (double m : mx) ub = std::max(ub, m);
    }
    if (!std::isfinite(ub)) { result = {PrStop::NumericalFailure, it}; break; }
    lb = gpr_prefix(d, k, rl, H, pos, sz, cnt, prefix_cap);
    if (prefix_cap < n2 && ((lb == 0 && ub > 0) ||
        (rescue_prefix && it == maxIt && !certificate_closed(ub, std::max(best_lb, lb), target)))) {
      prefix_cap = n2;
      lb = gpr_prefix(d, k, rl, H, pos, sz, cnt, prefix_cap);
    }
    if (ub == 0 && bestS.empty()) bestS.assign(1, 0);
    if (lb > best_lb) {
      best_lb = lb;
      best_cnt = cnt;
      bestS.resize(sz);
      for (u32 i = 0; i < sz; ++i) bestS[i] = old[H[i]];
    }
    if (it % 10 == 0)
      std::fprintf(stderr, "  global it=%d UB=%.4f LB=%.4f gap=%.5f t=%.1fs\n", it, ub, best_lb, ub / best_lb, t.sec());
    if (certificate_closed(ub, best_lb, target) || (it == 1 && certificate_closed(ub, best_lb, trigger))) {
      result = {PrStop::Certified, it}; break;
    }
    if (!(lb > 0) || !std::isfinite(lb)) { result = {PrStop::NumericalFailure, it}; break; }
    if (time_cap > 0 && ttot.sec() > time_cap) {
      finish();
      return {PrStop::TimeLimit, it};
    }
    if (it == maxIt) break;
    if (it == 1) gpr_split(d, tasks, work, T, nhead);

#pragma omp parallel for schedule(static)
    for (long long a = 0; a < (long long)n2; ++a)
      if (rl[a] > 0) {
        const u32 v = old[a];
        lam[v] = clampLam(lam[v] + std::log(rl[a] / lb), lam0);
      }
  }
  finish();
  return result;
}


struct HDump {
  FILE* f = nullptr;
  FILE* q = nullptr;
  int k = 0;
  u64 next_id = 0, steps = 0;
  using Key = std::array<u32, 6>;
  struct KeyHash {
    size_t operator()(const Key& a) const {
      u64 h = 0x9E3779B97F4A7C15ull;
      for (u32 x : a) h = mix64(h ^ x);
      return (size_t)h;
    }
  };
  std::unordered_map<Key, u64, KeyHash> ids;
  std::vector<char> fbuf;
  void open(const std::string& prefix, u32 n, int kk) {
    k = kk;
    if (k < 2 || k > 6) { std::fprintf(stderr, "KCDS_HDUMP supports 2 <= k <= 6\n"); std::exit(1); }
    f = std::fopen((prefix + ".updates").c_str(), "wb");
    q = std::fopen((prefix + ".queries.tsv").c_str(), "wb");
    if (!f || !q) { std::fprintf(stderr, "KCDS_HDUMP: cannot open %s.*\n", prefix.c_str()); std::exit(1); }
    fbuf.resize((size_t)1 << 26);
    std::setvbuf(f, fbuf.data(), _IOFBF, fbuf.size());
    std::fprintf(f, "%u\n", n);
  }
  void write(const std::vector<u32>& flat, char op) {
    for (size_t i = 0; i + (size_t)k <= flat.size(); i += (size_t)k) {
      Key a;
      a.fill(0xFFFFFFFFu);
      for (int j = 0; j < k; ++j) a[j] = flat[i + j];
      std::sort(a.begin(), a.begin() + k);
      u64 id;
      if (op == '+') {
        id = next_id++;
        if (!ids.emplace(a, id).second) { std::fprintf(stderr, "KCDS_HDUMP: clique inserted twice\n"); std::exit(1); }
      } else {
        auto it = ids.find(a);
        if (it == ids.end()) { std::fprintf(stderr, "KCDS_HDUMP: deleted clique not active\n"); std::exit(1); }
        id = it->second;
        ids.erase(it);
      }
      std::fprintf(f, "%c %llu %d", op, (unsigned long long)id, k);
      for (int j = 0; j < k; ++j) std::fprintf(f, " %u", a[j]);
      std::fputc('\n', f);
      ++steps;
    }
  }
  void query(const std::string& label) {
    std::fprintf(q, "%llu\t%s\n", (unsigned long long)(steps ? steps - 1 : 0), label.c_str());
  }
  void close() {
    std::fclose(f);
    std::fclose(q);
  }
};

#include "barrier.hpp"
static int run_driver(int argc, char** argv) {
  Timer trace_clock;
  const std::string arm = std::getenv("TRACE_ARM");
  const bool global_arm = arm == "cold" || arm == "warm";
  if (argc < 8) {
    std::fprintf(stderr, "usage: %s graph.bin k eps batch nbatches mode seed [check_every] [hold_frac]\n", argv[0]);
    return 1;
  }


  char* k_end = nullptr;
  const long requested_k = std::strtol(argv[2], &k_end, 10);
  if (k_end == argv[2] || *k_end != '\0' || requested_k < 2 || requested_k > 257)
    throw std::invalid_argument("clique size must be an integer in [2,257]");
  const int k = static_cast<int>(requested_k);
  const double eps = std::atof(argv[3]);
  const int batch = std::atoi(argv[4]);
  const int nb = std::atoi(argv[5]);
  const std::string mode = argv[6];
  const unsigned seed = (unsigned)std::atoi(argv[7]);
  if (!std::isfinite(eps) || eps <= 0 || 1 + 0.25 * eps == 1 || batch < 0 || nb < 0)
    throw std::invalid_argument("invalid eps or workload size");
  const int check_every = argc > 8 ? std::atoi(argv[8]) : 5;
  const double hold_frac = argc > 9 ? std::atof(argv[9]) : 0.1;
  std::mt19937_64 rng(seed);
  const int T = omp_get_max_threads();

  Solver s;
  s.k = k;
  s.eps = eps;
  u32 n;
  std::vector<std::pair<u32, u32>> pool, cur;


  const bool stream = mode == "stream";


  const bool events = mode == "events", streamlike = stream || events;
  std::unordered_set<u64> dlinks;


  bool compact = std::getenv("KCDS_COMPACT") != nullptr;
  Stream st;
  size_t tail = 0, head = 0;
  std::unordered_map<u64, u32> mult;
  auto apply_event = [&](size_t i) {
    const u32 a = st.u[i], b = st.v[i];
    const u64 dk = (u64)a * st.n + b, uk = ekey(a, b, st.n);
    if (st.op[i] > 0) {
      if (dlinks.insert(dk).second) mult[uk]++;
    } else if (dlinks.erase(dk)) {
      if (--mult[uk] == 0) mult.erase(uk);
    }
  };
  {
    if (events) {
      st = load_events(argv[1]);
      n = st.n;
      const size_t need = (size_t)nb * batch;
      head = st.M > need ? st.M - need : 0;
      for (size_t i = 0; i < head; ++i) apply_event(i);
      cur.reserve(mult.size());
      for (const auto& [key, c] : mult) cur.push_back({(u32)(key / n), (u32)(key % n)});
      std::sort(cur.begin(), cur.end());
      std::printf("events total=%llu initial=%zu edges=%zu links=%zu\n", (unsigned long long)st.M, head, cur.size(),
                  dlinks.size());
    } else if (stream) {
      st = load_stream(argv[1]);
      n = st.n;
      const size_t W = (size_t)hold_frac, need = W + (size_t)nb * batch;
      tail = st.M > need ? st.M - need : 0;
      head = std::min<size_t>(st.M, tail + W);
      mult.reserve(2 * W);
      for (size_t i = tail; i < head; ++i)
        if (mult[ekey(st.u[i], st.v[i], n)]++ == 0) cur.push_back({st.u[i], st.v[i]});
      std::printf("stream events=%llu window=[%zu,%zu) distinct_pairs=%zu\n", (unsigned long long)st.M, tail, head,
                  cur.size());
    } else {
      if (!compact) {
        FILE* f = std::fopen(argv[1], "rb");
        u64 h[2] = {0, 0};
        if (f) { if (std::fread(h, 8, 2, f) != 2) h[0] = 0; std::fclose(f); }
        compact = h[0] > 30000000ull;
      }
      if (compact && mode != "replay") { std::fprintf(stderr, "compact mode supports replay and stream only\n"); return 1; }
      Graph g0 = load_bin(argv[1]);
      n = g0.n;
      if (compact) {


        const u64 thr = (u64)(hold_frac * 18446744073709551615.0), salt = mix64(seed);
        std::vector<u64> noff((size_t)n + 1);
        u64 wpos = 0;
        for (u32 v = 0; v < n; ++v) {
          const u64 s0 = g0.off[v], s1 = g0.off[v + 1];
          noff[v] = wpos;
          for (u64 e = s0; e < s1; ++e) {
            const u32 u = g0.adj[e];
            const u32 a = std::min(u, v), b = std::max(u, v);
            if (mix64(((u64)a << 32 | b) ^ salt) < thr) {
              if (v < u) pool.push_back({v, u});
              continue;
            }
            g0.adj[wpos++] = u;
          }
        }
        noff[n] = wpos;
        g0.off.swap(noff);
        g0.adj.resize(wpos);
        g0.m = wpos / 2;
        std::shuffle(pool.begin(), pool.end(), rng);
        s.G.init(g0, true);
      } else {
        for (u32 v = 0; v < n; ++v)
          for (u32 u : std::vector<u32>(g0.nb(v), g0.nb(v) + g0.deg(v))) if (u > v) cur.push_back({v, u});
        std::shuffle(cur.begin(), cur.end(), rng);
        const size_t hold = (size_t)(cur.size() * hold_frac);
        pool.assign(cur.end() - hold, cur.end());
        cur.resize(cur.size() - hold);
      }
    }
    if (!(compact && !streamlike)) {
      std::vector<std::vector<u32>> lists(n);
      for (auto [u, v] : cur) { lists[u].push_back(v); lists[v].push_back(u); }
      for (auto& l : lists) std::sort(l.begin(), l.end());
      Graph g1 = csr_from_lists(lists);
      s.G.init(g1);
    }
  }
  if (n == 0) throw std::invalid_argument("a non-empty vertex universe is required");
  if (check_every < 1) throw std::invalid_argument("check_every must be positive");
  std::printf("mode=%s compact=%d\n", mode.c_str(), (int)compact);


  const bool rerank_on = !compact && std::getenv("KCDS_NORERANK") == nullptr;
  if (rerank_on) s.G.set_limit();
  int nrerank = 0;


  Fenwick fen;
  if (compact && !streamlike) {
    std::vector<u64> deg(n);
    for (u32 v = 0; v < n; ++v) deg[v] = s.G.adj[v].size();
    fen.build(deg);
  }
  std::unordered_map<u64, size_t> where;
  if (!streamlike && !compact) {
    where.reserve(cur.size() * 2);
    for (size_t i = 0; i < cur.size(); ++i) where[ekey(cur[i].first, cur[i].second, n)] = i;
  }


  std::vector<std::pair<u32, u32>> plant;
  if (mode == "plant") {
    if (compact) { std::fprintf(stderr, "plant mode needs the edge index (no compact mode)\n"); return 1; }
    const size_t need = (size_t)(batch - batch / 2) * (size_t)nb;
    u32 p = 2;
    while ((size_t)p * (p - 1) / 2 < need + need / 10) ++p;
    std::vector<u32> P;
    std::unordered_set<u32> inP;
    while (P.size() < p) {
      const u32 v = (u32)(rng() % n);
      if (inP.insert(v).second) P.push_back(v);
    }
    for (size_t i = 0; i < P.size(); ++i)
      for (size_t j = i + 1; j < P.size(); ++j) {
        const u32 a = std::min(P[i], P[j]), b = std::max(P[i], P[j]);
        if (!where.count(ekey(a, b, n))) plant.push_back({a, b});
      }
    std::shuffle(plant.begin(), plant.end(), rng);
    std::printf("plant |P|=%u missing_pairs=%zu\n", p, plant.size());
  }

  s.lam.assign(n, 0.0);
  s.w.assign(n, 1.0);
  s.r.assign(n, 0.0);
  s.inS.assign(n, 0);
  s.isChanged.assign(n, 0);
  s.xidx.assign(n, -1);
  s.hpos.assign(n, 0);
  s.hstamp.assign(n, 0);
  s.inR.assign(n, 0);
  s.rstamp.assign(n, 0);
  {
    const std::string st = argc > 10 ? argv[10] : "targeted";
    s.strategy = st == "one" ? 0 : (st == "two" ? 1 : 2);
    s.two_sided = s.strategy == 1;
    s.best_response = st == "br";
    s.use_newton = st == "newton" || st == "hybrid";
    s.hybrid = st == "hybrid";
    if (argc > 11) s.q = std::atoi(argv[11]);
    if (argc > 12) s.damp = std::atof(argv[12]);
  }
  s.ts.resize(T);
  for (auto& t : s.ts) t.init(n, k, compact);
  s.hcnt.assign(T, std::vector<double>(Solver::kNB, 0.0));
  s.hused.assign(T, {});
  s.absorb.assign(n, 0.0);
  s.absflag.assign(n, 0);
  s.abslist.assign(T, {});
  s.vwork.assign(n, 0.0f);


  const double init_gap = std::getenv("KCDS_INIT_GAP") ? std::atof(std::getenv("KCDS_INIT_GAP")) : 0.0;

  auto cold_solve = [&](const char* tag) {
    Timer tc;
    const double target = std::string(tag) == "cold_batch" ? 1 + eps : 1 + 0.5 * eps;
    char solve_line[2048] = {}, finish_line[1024] = {};
    std::fill(s.lam.begin(), s.lam.end(), s.lam0);
    std::fill(s.w.begin(), s.w.end(), 1.0);
    for (u32 v : s.S) s.inS[v] = 0;
    s.S.clear();
    if (compact) {
      double LB = 0, ti = 0;
      u64 cnt = 0, m2 = 0;
      for (u32 v = 0; v < n; ++v) m2 += s.G.adj[v].size();
      const PrResult solve = s.global_pr_c(init_gap > 0 ? init_gap : target, s.lam, s.w, s.r, LB, s.S, cnt,
                                   kSolveMaxIt, &ti, 0, 0, true);
      require_solve(solve, tag);
      const int it = solve.iterations;
      std::snprintf(solve_line, sizeof solve_line, "%s graph=%s n=%u m=%llu k=%d eps=%.3f threads=%d iters=%d prep=%.2fs iter=%.2fs LB=%.4f |S|=%zu\n",
                  tag, argv[1], n, (unsigned long long)(m2 / 2), k, eps, T, it, 0.0, ti, LB, s.S.size());
      for (u32 v : s.S) s.inS[v] = 1;
      s.cS = (double)cnt;
    } else {
      Graph gs = csr_from_lists(s.G.adj);
      double LB = 0, tp = 0, ti = 0;
      u64 cnt = 0;
      const PrResult solve = global_pr(gs, k, init_gap > 0 ? init_gap : target, s.lam, s.lam0, s.w, s.r, LB, s.S,
                               kSolveMaxIt, &tp, &ti, 0, 0, &cnt, true);
      require_solve(solve, tag);
      const int it = solve.iterations;
      std::snprintf(solve_line, sizeof solve_line, "%s graph=%s n=%u m=%llu k=%d eps=%.3f threads=%d iters=%d prep=%.2fs iter=%.2fs LB=%.4f |S|=%zu\n",
                  tag, argv[1], n, (unsigned long long)gs.m, k, eps, T, it, tp, ti, LB, s.S.size());
      for (u32 v : s.S) s.inS[v] = 1;
      s.cS = (double)cnt;
    }
    s.tree.build(s.r);
    if (init_gap > 0) {
      Timer tf;
      const double b0 = s.budget_s;
      s.budget_s = 0;
      const u64 c0 = s.cliques_touched;
      const int rr = s.repair(1000);
      s.budget_s = b0;
      std::snprintf(finish_line, sizeof finish_line, "%s_finish rounds=%d time=%.2fs cliques=%llu UB=%.4f LB=%.4f gap=%.5f |S|=%zu\n", tag, rr, tf.sec(),
                  (unsigned long long)(s.cliques_touched - c0), s.tree.max(), s.rhoS(), certificate_gap(s.tree.max(), s.rhoS()),
                  s.S.size());
    }
    require_certificate(s.tree.max(), s.rhoS(), target, s.S.size(), tag);
    std::fputs(solve_line, stdout);
    std::fputs(finish_line, stdout);
    s.t_refresh = s.t_relax = s.t_hess = s.t_solve = 0;
    s.n_refresh = 0;
    s.n_newton = 0;
    s.n_absorb = 0;
    s.n_split = 0;
    return tc.sec();
  };


  auto warm_solve = [&](double cap, double target) {
    Timer tw;
    for (u32 v : s.S) s.inS[v] = 0;
    double LB = 0, tp = 0, ti = 0;
    PrResult solve{PrStop::NumericalFailure, 0};
    if (compact) {
      u64 cnt = 0;
      solve = s.global_pr_c(target, s.lam, s.w, s.r, LB, s.S, cnt, 300, &ti, 0, cap);
      s.cS = (double)cnt;
    } else {
      Graph gs = csr_from_lists(s.G.adj);
      u64 cnt = 0;
      solve = global_pr(gs, k, target, s.lam, s.lam0, s.w, s.r, LB, s.S, kBaselineMaxIt, &tp, &ti, 0, cap, &cnt);
      s.cS = (double)cnt;
    }
    for (u32 v : s.S) s.inS[v] = 1;
    s.tree.build(s.r);
    const int it = solve.log_iterations();
    const bool ok = solve.certified() && certificate_closed(s.tree.max(), s.rhoS(), target);
    std::printf("  warm_fallback iters=%d prep=%.2fs iter=%.2fs time=%.2fs UB=%.4f LB=%.4f closed=%d\n", it, tp, ti,
                tw.sec(), s.tree.max(), s.rhoS(), (int)ok);
    return ok;
  };
  const char* hdump_prefix = std::getenv("KCDS_HDUMP");
  HDump hd;
  double t_cold = 0;
  if (hdump_prefix) {
    s.hdump = true;
    s.tree.build(s.r);
    hd.open(hdump_prefix, n, k);
    Graph g0 = csr_from_lists(s.G.adj);
    std::vector<u32> order0;
    core_decomposition(g0, &order0);
    DAG d0 = orient(g0, order0);
    std::vector<std::vector<u32>> parts(T);
    for_all_cliques(d0, k, [&](int tid) {
      std::vector<u32>* p = &parts[tid];
      return [p](const std::vector<u32>& st, const std::vector<u32>& C) {
        for (u32 x : C) {
          p->insert(p->end(), st.begin(), st.end());
          p->push_back(x);
        }
      };
    });
    for (auto& p : parts) {
      hd.write(p, '+');
      std::vector<u32>().swap(p);
    }
    hd.query("init");
    std::printf("hdump init cliques=%llu\n", (unsigned long long)hd.steps);
  } else {
    t_cold = cold_solve("init");
  }
  checkpoint_barrier(s, 0, trace_clock, t_cold, 0, 0, 0);
  s.cliques_touched = 0;


  const bool warm_first = std::getenv("KCDS_WARMFB") != nullptr && std::getenv("KCDS_COLDFB") == nullptr;

  const double fb = std::getenv("KCDS_FALLBACK") ? std::atof(std::getenv("KCDS_FALLBACK")) : 1.0;
  s.budget_s = fb > 0 ? fb * t_cold : 0;
  std::printf("init_total time=%.2fs repair_budget=%.2fs\n", t_cold, s.budget_s);


  auto drift = [&]() {
    u32 maxout = 0;
    double lmin = 1e300, lmax = -1e300;
    for (u32 v = 0; v < n; ++v) {
      maxout = std::max<u32>(maxout, (u32)s.G.out[v].size());
      if (s.r[v] > 0) { lmin = std::min(lmin, s.lam[v]); lmax = std::max(lmax, s.lam[v]); }
    }
    char buf[96];
    std::snprintf(buf, sizeof buf, "maxout=%u lam_range=%.3f", maxout, lmax >= lmin ? lmax - lmin : 0.0);
    return std::string(buf);
  };
  std::printf("init_drift %s\n", drift().c_str());

  double tot = 0;
  int nrepairs = 0;

  const int last = std::getenv("KCDS_MAX_BATCH") ? std::min(nb, std::atoi(std::getenv("KCDS_MAX_BATCH"))) : nb;


  const bool nocheck = std::getenv("KCDS_NOCHECK") != nullptr;
  for (int b = 1; b <= last; ++b) {
    Timer batch_clock;
    std::vector<std::pair<u32, u32>> dels, inss;
    const char* del_frac = std::getenv("KCDS_DEL_FRAC");
    const int nd = streamlike ? 0
                   : del_frac ? std::max(0, std::min(batch, (int)std::lround(batch * std::atof(del_frac))))
                              : batch / 2,
              ni = streamlike ? 0 : batch - nd;
    if (events) {

      std::unordered_map<u64, char> before;
      for (int i = 0; i < batch && head < st.M; ++i, ++head) {
        const u64 key = ekey(st.u[head], st.v[head], n);
        if (!before.count(key)) before[key] = mult.count(key) ? 1 : 0;
        apply_event(head);
      }
      for (const auto& [key, was] : before) {
        const bool now = mult.count(key) > 0;
        const std::pair<u32, u32> e{(u32)(key / n), (u32)(key % n)};
        if (was && !now) dels.push_back(e);
        if (!was && now) inss.push_back(e);
      }
      std::sort(dels.begin(), dels.end());
      std::sort(inss.begin(), inss.end());
    }
    if (stream) {

      std::unordered_map<u64, std::pair<u32, u32>> seen;
      std::unordered_map<u64, char> before;
      auto touch = [&](size_t i) {
        const u64 key = ekey(st.u[i], st.v[i], n);
        if (seen.emplace(key, std::make_pair(st.u[i], st.v[i])).second) before[key] = mult.count(key) ? 1 : 0;
        return key;
      };
      for (int i = 0; i < batch && head < st.M; ++i, ++head) mult[touch(head)]++;
      for (int i = 0; i < batch && tail < head; ++i, ++tail) {
        const u64 key = touch(tail);
        if (--mult[key] == 0) mult.erase(key);
      }
      for (auto& [key, e] : seen) {
        const bool now = mult.count(key) > 0;
        if (before[key] && !now) dels.push_back(e);
        if (!before[key] && now) inss.push_back(e);
      }
    }
    if (compact && !streamlike) {

      std::unordered_set<u64> chosen;
      const u64 tot2 = fen.total();
      for (int i = 0; i < nd && tot2 > 0; ++i) {
        for (int tries = 0; tries < 100; ++tries) {
          const u32 v = fen.find(rng() % tot2);
          const auto& a = s.G.adj[v];
          if (a.empty()) continue;
          const u32 u = a[rng() % a.size()];
          if (chosen.insert(ekey(u, v, n)).second) { dels.push_back({std::min(u, v), std::max(u, v)}); break; }
        }
      }
    }
    for (int i = 0; i < nd && !cur.empty() && !compact; ++i) {
      size_t id;
      if ((mode == "core" || mode == "plant") && s.S.size() > 1) {
        bool found = false;
        std::pair<u32, u32> e{0, 0};
        for (int tries = 0; tries < 50 && !found; ++tries) {
          const u32 u = s.S[rng() % s.S.size()];
          const auto& a = s.G.adj[u];
          if (a.empty()) continue;
          const u32 v = a[rng() % a.size()];
          if (s.inS[v]) { e = {std::min(u, v), std::max(u, v)}; found = true; }
        }
        if (!found) continue;
        auto it = where.find(ekey(e.first, e.second, n));
        if (it == where.end()) continue;
        id = it->second;
      } else {
        id = rng() % cur.size();
      }
      auto e = cur[id];
      where.erase(ekey(e.first, e.second, n));
      if (id != cur.size() - 1) { cur[id] = cur.back(); where[ekey(cur[id].first, cur[id].second, n)] = id; }
      cur.pop_back();
      dels.push_back(e);
    }
    std::vector<std::pair<u32, u32>>& src = mode == "plant" ? plant : pool;
    for (int i = 0; i < ni && !src.empty(); ++i) {
      auto e = src.back();
      src.pop_back();
      inss.push_back(e);
      if (compact) continue;
      where[ekey(e.first, e.second, n)] = cur.size();
      cur.push_back(e);
    }
    s.cliques_touched = 0;
    Timer tu;
    if (global_arm) {
      for (auto [u, v] : dels) {
        DynGraph::del_sorted(s.G.adj[u], v);
        DynGraph::del_sorted(s.G.adj[v], u);
      }
      for (auto [u, v] : inss) {
        DynGraph::ins_sorted(s.G.adj[u], v);
        DynGraph::ins_sorted(s.G.adj[v], u);
      }
      if (arm != "warm" || !warm_solve(s.budget_s, 1 + eps)) cold_solve("cold_batch");
      checkpoint_barrier(s, b, batch_clock, tu.sec(), 0, dels.size(), inss.size());
      continue;
    }
    s.apply_batch(dels, inss);
    if (s.G.rerank_due) {
      s.G.rerank();
      ++nrerank;
      std::printf("  rerank batch=%d out_limit=%u\n", b, s.G.out_limit);
    }
    const double t_upd = tu.sec();
    if (compact && !streamlike) {
      for (auto [u, v] : dels) { fen.add(u, -1); fen.add(v, -1); }
      for (auto [u, v] : inss) { fen.add(u, +1); fen.add(v, +1); }
    }
    const u64 c_upd = s.cliques_touched;
    s.cliques_touched = 0;
    if (hdump_prefix) {
      const size_t nd_c = s.hd_del.size() / k, ni_c = s.hd_ins.size() / k;
      hd.write(s.hd_del, '-');
      hd.write(s.hd_ins, '+');
      s.hd_del.clear();
      s.hd_ins.clear();
      hd.query("b" + std::to_string(b));
      std::printf("hdump batch=%d dels=%zu ins=%zu cliques_del=%zu cliques_ins=%zu steps=%llu\n", b, dels.size(),
                  inss.size(), nd_c, ni_c, (unsigned long long)hd.steps);
      std::fflush(stdout);
      continue;
    }
    const bool check = !nocheck && (b % check_every == 0 || b == last);
    std::vector<double> lam_pre;
    if (check) lam_pre = s.lam;
    Timer tr;
    int rounds = s.repair();
    if (rounds < 0 && fb > 0) {
      std::printf("  fallback batch=%d reason=%s repair_spent=%.3fs\n", b, rounds == -2 ? "budget" : "cap", tr.sec());
      if (!(warm_first && warm_solve(s.budget_s, 1 + 0.5 * eps))) cold_solve("fallback");
      rounds = -2;
    }
    if (rounds < 0 && fb <= 0) throw std::runtime_error("uncertified repair; fallback disabled");
    require_certificate(s.tree.max(), s.rhoS(), 1 + eps, s.S.size(), "batch");
    const double t_rep = tr.sec();
    if (rounds != 0) nrepairs++;
    if (rounds != 0)
      std::printf("  repair breakdown: refresh=%.3fs (%d calls) relax=%.3fs (hessian=%.3fs solve=%.3fs) newton_rounds=%d "
                  "absorbed=%zu split=%zu\n",
                  s.t_refresh, s.n_refresh, s.t_relax, s.t_hess, s.t_solve, s.n_newton, s.n_absorb, s.n_split);
    s.t_refresh = s.t_relax = s.t_hess = s.t_solve = 0; s.n_refresh = 0; s.n_newton = 0; s.n_absorb = 0;
    s.n_split = 0;
    tot += t_upd + t_rep;
    const double ub = s.tree.max(), lb = s.rhoS();
    char batch_line[1024];
    std::snprintf(batch_line, sizeof batch_line, "batch=%d dels=%zu ins=%zu upd=%.4fs cliques_upd=%llu repair=%.4fs rounds=%d cliques_rep=%llu UB=%.4f "
                "LB=%.4f gap=%.5f |S|=%zu\n",
                b, dels.size(), inss.size(), t_upd, (unsigned long long)c_upd, t_rep, rounds,
                (unsigned long long)s.cliques_touched, ub, lb, certificate_gap(ub, lb), s.S.size());
    if (check && b <= last && compact) {

      std::vector<double> rex(n);
      s.dense_loads(s.w, rex);
      double maxerr = 0, ubx = 0;
      for (u32 v = 0; v < n; ++v) {
        if (!std::isfinite(rex[v])) throw std::runtime_error("non-finite recounted load");
        maxerr = std::max(maxerr, std::fabs(rex[v] - s.r[v]) / std::max(1.0, rex[v]));
        ubx = std::max(ubx, rex[v]);
      }
      std::vector<double>().swap(rex);
      const double exactS = (double)s.count_inside(s.S);
      std::vector<double> lamw = lam_pre, ww(n), rw(n);
      double LBw = 0, LBc = 0, iw = 0, ic = 0;
      std::vector<u32> Sw, Sc;
      u64 cw = 0, cc = 0;
      const PrResult sw = s.global_pr_c(1 + 0.5 * eps, lamw, ww, rw, LBw, Sw, cw, kBaselineMaxIt, &iw, 1 + eps);
      std::vector<double>().swap(lamw);
      std::vector<double>().swap(ww);
      std::vector<double>().swap(rw);
      std::vector<double> lamc(n, s.lam0), wc(n, 1.0), rc(n);
      const PrResult sc = s.global_pr_c(1 + eps, lamc, wc, rc, LBc, Sc, cc, kBaselineMaxIt, &ic);
      const int itw = sw.log_iterations(), itc = sc.log_iterations();
      std::printf("baseline_status batch=%d warm=%s cold=%s\n", b, sw.name(), sc.name());
      if (exactS != s.cS) throw std::runtime_error("witness recount mismatch");
      require_certificate(ubx, exactS / s.S.size(), 1 + eps, s.S.size(), "recount");


      std::printf("check batch=%d load_relerr=%.2e cS=%.0f exact=%.0f | warm_global iters=%d prep=%.2fs iter=%.2fs "
                  "LB=%.4f | cold_global iters=%d iter=%.2fs LB=%.4f | ours avg/batch=%.4fs repairs=%d/%d | "
                  "exact_ub=%.6f size=%zu cold_prep=%.2fs %s cliques=0\n",
                  b, maxerr, s.cS, exactS, itw, 0.0, iw, LBw, itc, ic, LBc, tot / b, nrepairs, b, ubx, s.S.size(), 0.0,
                  drift().c_str());
      std::fflush(stdout);
    }
    if (check && b <= last && !compact) {
      Graph gs = csr_from_lists(s.G.adj);
      if (const char* dump = std::getenv("KCDS_DUMP")) {

        const std::string path = std::string(dump) + "_b" + std::to_string(b) + ".txt";
        FILE* f = std::fopen(path.c_str(), "w");
        if (f) {
          std::fprintf(f, "%u %llu\n", gs.n, (unsigned long long)gs.m);
          for (u32 v = 0; v < gs.n; ++v)
            for (u64 e = gs.off[v]; e < gs.off[v + 1]; ++e)
              if (gs.adj[e] > v) std::fprintf(f, "%u %u\n", v, gs.adj[e]);
          std::fclose(f);
          std::printf("dumped %s\n", path.c_str());
        }
      }
      std::vector<u32> order;
      core_decomposition(gs, &order);
      DAG d = orient(gs, order);
      std::vector<std::vector<double>> rb(T, std::vector<double>(n, 0.0));
      std::vector<double> rex(n);
      potential_loads(d, k, s.w, rb, rex);
      double maxerr = 0, ubx = 0;
      for (u32 v = 0; v < n; ++v) {
        if (!std::isfinite(rex[v])) throw std::runtime_error("non-finite recounted load");
        maxerr = std::max(maxerr, std::fabs(rex[v] - s.r[v]) / std::max(1.0, rex[v]));
        ubx = std::max(ubx, rex[v]);
      }
      std::vector<std::vector<double>> cnt(T, std::vector<double>(2, 0.0));
      for_all_cliques(d, k, [&](int tid) {
        double* c = &cnt[tid][0];
        const std::vector<char>* in = &s.inS;
        return [c, in](const std::vector<u32>& st, const std::vector<u32>& C) {
          c[1] += (double)C.size();
          for (u32 y : st) if (!(*in)[y]) return;
          for (u32 x : C) if ((*in)[x]) c[0] += 1.0;
        };
      });
      double exactS = 0, allC = 0;
      for (auto& c : cnt) { exactS += c[0]; allC += c[1]; }


      std::vector<double> lamw = lam_pre, ww(n), rw(n), lamc(n, 0.0), wc(n, 1.0), rc(n);
      double LBw, LBc, pw, iw, pc, ic;
      std::vector<u32> Sw, Sc;
      const PrResult sw = global_pr(gs, k, 1 + 0.5 * eps, lamw, s.lam0, ww, rw, LBw, Sw, kBaselineMaxIt, &pw, &iw, 1 + eps);
      const PrResult sc = global_pr(gs, k, 1 + eps, lamc, 0.0, wc, rc, LBc, Sc, kBaselineMaxIt, &pc, &ic);
      const int itw = sw.log_iterations(), itc = sc.log_iterations();
      std::printf("baseline_status batch=%d warm=%s cold=%s\n", b, sw.name(), sc.name());
      if (exactS != s.cS) throw std::runtime_error("witness recount mismatch");
      require_certificate(ubx, exactS / s.S.size(), 1 + eps, s.S.size(), "recount");

      std::printf("check batch=%d load_relerr=%.2e cS=%.0f exact=%.0f | warm_global iters=%d prep=%.2fs iter=%.2fs "
                  "LB=%.4f | cold_global iters=%d iter=%.2fs LB=%.4f | ours avg/batch=%.4fs repairs=%d/%d | "
                  "exact_ub=%.6f size=%zu cold_prep=%.2fs %s cliques=%.0f\n",
                  b, maxerr, s.cS, exactS, itw, pw, iw, LBw, itc, ic, LBc, tot / b, nrepairs, b, ubx, s.S.size(), pc,
                  drift().c_str(), allC);
      std::fflush(stdout);
    }

    checkpoint_barrier(s, b, batch_clock, t_upd + t_rep, rounds, dels.size(), inss.size());
    std::fputs(batch_line, stdout);
    std::fflush(stdout);
  }
  std::printf("reranks total=%d\n", nrerank);
  if (hdump_prefix) {
    hd.close();
    std::printf("hdump done updates=%llu active=%zu\n", (unsigned long long)hd.steps, hd.ids.size());
    return 0;
  }
  {


    double graph_b = 0, vec_b = 0, thr_b = 0;
    for (u32 v = 0; v < n; ++v)
      graph_b += 4.0 * (s.G.adj[v].capacity() + s.G.out[v].capacity()) + 2 * sizeof(std::vector<u32>);
    graph_b += 4.0 * s.G.rank.capacity();
    vec_b = 8.0 * (s.lam.capacity() + s.w.capacity() + s.r.capacity() + s.tree.t.capacity() + s.absorb.capacity()) +
            4.0 * (s.S.capacity() + s.hpos.capacity() + s.hstamp.capacity() + s.rstamp.capacity() + s.xidx.capacity()) +
            1.0 * (s.inS.capacity() + s.isChanged.capacity() + s.inR.capacity() + s.absflag.capacity() +
                   s.bmark.capacity());
    if (compact) vec_b += 4.0 * s.vwork.capacity() + 8.0 * fen.t.capacity();
    for (auto& t : s.ts)
      thr_b += 8.0 * t.dr.capacity() + 1.0 * (t.tflag.capacity() + t.glab.capacity()) +
               4.0 * (t.le.stamp.capacity() + t.le.loc.capacity() + t.le.hkey.capacity() + t.le.hval.capacity());
    for (auto& h : s.hcnt) thr_b += 8.0 * h.capacity();
    PROCESS_MEMORY_COUNTERS pmc{};
    GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
    std::printf("memory state_graph=%.1fMB state_vectors=%.1fMB state_thread_scratch=%.1fMB process_peak=%.1fMB\n",
                graph_b / 1048576, vec_b / 1048576, thr_b / 1048576, pmc.PeakWorkingSetSize / 1048576.0);
  }
  return 0;
}

int main(int argc, char** argv) {
  try { return run_driver(argc, argv); }
  catch (const std::exception& e) {
    std::fflush(stdout);
    std::fprintf(stderr, "KCDS_FAILURE: %s\n", e.what());
    return 2;
  }
}
