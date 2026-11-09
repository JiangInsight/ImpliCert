#pragma once
#include <omp.h>
#include "graph.hpp"


template <class Leaf>
static void kc_rec(const DAG& d, std::vector<unsigned char>& lab, std::vector<std::vector<u32>>& cand,
                   std::vector<u32>& stack, int level, int l, Leaf& leaf) {
  const std::vector<u32>& C = cand[level];
  if (l == 1) { leaf(stack, C); return; }
  const unsigned char lv = static_cast<unsigned char>(level), nl = static_cast<unsigned char>(level + 1);
  for (u32 x : C) {
    std::vector<u32>& next = cand[level + 1];
    next.clear();
    for (u64 e = d.off[x]; e < d.off[x + 1]; ++e)
      if (lab[d.adj[e]] == lv) next.push_back(d.adj[e]);
    if (next.size() + 1 < static_cast<size_t>(l)) continue;
    for (u32 y : next) lab[y] = nl;
    stack.push_back(x);
    kc_rec(d, lab, cand, stack, level + 1, l - 1, leaf);
    stack.pop_back();
    for (u32 y : cand[level + 1]) lab[y] = lv;
  }
}


template <class MakeLeaf>
static void for_all_cliques(const DAG& d, int k, MakeLeaf&& make_leaf) {
#pragma omp parallel
  {
    auto leaf = make_leaf(omp_get_thread_num());
    std::vector<unsigned char> lab(d.n, 255);
    std::vector<std::vector<u32>> cand(k + 1);
    std::vector<u32> stack;
#pragma omp for schedule(dynamic, 64)
    for (long long v = 0; v < (long long)d.n; ++v) {
      const u32 u = static_cast<u32>(v);
      if (d.outdeg(u) + 1 < static_cast<u32>(k)) continue;
      cand[0].assign(d.out(u), d.out(u) + d.outdeg(u));
      for (u32 w : cand[0]) lab[w] = 0;
      stack.assign(1, u);
      kc_rec(d, lab, cand, stack, 0, k - 1, leaf);
      for (u32 w : cand[0]) lab[w] = 255;
    }
  }
}


template <class T>
static void reduce_into(std::vector<std::vector<T>>& buf, std::vector<T>& out) {
  const long long n = (long long)out.size();
#pragma omp parallel for schedule(static)
  for (long long i = 0; i < n; ++i) {
    T s = 0;
    for (auto& b : buf) { s += b[i]; b[i] = 0; }
    out[i] = s;
  }
}


static void potential_loads(const DAG& d, int k, const std::vector<double>& w, std::vector<std::vector<double>>& rb,
                            std::vector<double>& r) {
  for_all_cliques(d, k, [&](int tid) {
    std::vector<double>* buf = &rb[tid];
    return [buf, &w](const std::vector<u32>& st, const std::vector<u32>& C) {
      double s0 = 0;
      for (u32 y : st) s0 += w[y];
      double a = 0;
      for (u32 x : C) {
        const double inv = 1.0 / (s0 + w[x]);
        a += inv;
        (*buf)[x] += w[x] * inv;
      }
      for (u32 y : st) (*buf)[y] += w[y] * a;
    };
  });
  reduce_into(rb, r);
}


static std::pair<double, u32> best_prefix(const DAG& d, int k, const std::vector<double>& score,
                                          std::vector<std::vector<u64>>& hb, std::vector<u32>* order_out = nullptr) {
  const u32 n = d.n;
  std::vector<u32> ord(n);
  for (u32 i = 0; i < n; ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(), [&](u32 a, u32 b) { return score[a] > score[b] || (score[a] == score[b] && a < b); });
  std::vector<u32> pos(n);
  for (u32 i = 0; i < n; ++i) pos[ord[i]] = i;
  for_all_cliques(d, k, [&](int tid) {
    std::vector<u64>* h = &hb[tid];
    return [h, &pos](const std::vector<u32>& st, const std::vector<u32>& C) {
      u32 p0 = 0;
      for (u32 y : st) p0 = std::max(p0, pos[y]);
      for (u32 x : C) (*h)[std::max(p0, pos[x])]++;
    };
  });
  std::vector<u64> hist(n);
  reduce_into(hb, hist);
  double best = 0;
  u32 bs = 0;
  u64 acc = 0;
  for (u32 i = 0; i < n; ++i) {
    acc += hist[i];
    double dens = (double)acc / (double)(i + 1);
    if (dens > best) { best = dens; bs = i + 1; }
  }
  if (order_out) *order_out = std::move(ord);
  return {best, bs};
}


static Graph csr_from_lists(const std::vector<std::vector<u32>>& adj) {
  Graph g;
  g.n = static_cast<u32>(adj.size());
  g.off.assign(g.n + 1, 0);
  for (u32 v = 0; v < g.n; ++v) g.off[v + 1] = g.off[v] + adj[v].size();
  g.adj.resize(g.off[g.n]);
  for (u32 v = 0; v < g.n; ++v) std::copy(adj[v].begin(), adj[v].end(), g.adj.begin() + g.off[v]);
  g.m = g.off[g.n] / 2;
  return g;
}
