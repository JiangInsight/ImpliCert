#define _CRT_SECURE_NO_WARNINGS
#pragma once
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

using u32 = uint32_t;
using u64 = uint64_t;

struct Timer {
  std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
  double sec() const { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); }
  void reset() { t0 = std::chrono::steady_clock::now(); }
};


struct Graph {
  u32 n = 0;
  u64 m = 0;
  std::vector<u64> off;
  std::vector<u32> adj;
  u32 deg(u32 v) const { return static_cast<u32>(off[v + 1] - off[v]); }
  const u32* nb(u32 v) const { return adj.data() + off[v]; }
};

inline Graph load_bin(const std::string& path) {
  FILE* f = std::fopen(path.c_str(), "rb");
  if (!f) throw std::runtime_error("cannot open " + path);
  u64 hdr[2];
  if (std::fread(hdr, sizeof(u64), 2, f) != 2) throw std::runtime_error("bad header");
  Graph g;
  g.n = static_cast<u32>(hdr[0]);
  g.m = hdr[1];
  g.off.resize(g.n + 1);
  g.adj.resize(2 * g.m);
  if (std::fread(g.off.data(), sizeof(u64), g.n + 1, f) != g.n + 1) throw std::runtime_error("bad offsets");
  if (std::fread(g.adj.data(), sizeof(u32), 2 * g.m, f) != 2 * g.m) throw std::runtime_error("bad adj");
  std::fclose(f);
  return g;
}


inline std::vector<u32> core_decomposition(const Graph& g, std::vector<u32>* order = nullptr) {
  const u32 n = g.n;
  std::vector<u32> deg(n), pos(n), vert(n);
  u32 md = 0;
  for (u32 v = 0; v < n; ++v) { deg[v] = g.deg(v); md = std::max(md, deg[v]); }
  std::vector<u32> bin(md + 1, 0);
  for (u32 v = 0; v < n; ++v) bin[deg[v]]++;
  u32 start = 0;
  for (u32 d = 0; d <= md; ++d) { u32 num = bin[d]; bin[d] = start; start += num; }
  for (u32 v = 0; v < n; ++v) { pos[v] = bin[deg[v]]; vert[pos[v]] = v; bin[deg[v]]++; }
  for (u32 d = md; d >= 1; --d) bin[d] = bin[d - 1];
  bin[0] = 0;
  for (u32 i = 0; i < n; ++i) {
    u32 v = vert[i];
    for (u64 e = g.off[v]; e < g.off[v + 1]; ++e) {
      u32 u = g.adj[e];
      if (deg[u] > deg[v]) {
        u32 du = deg[u], pu = pos[u], pw = bin[du], w = vert[pw];
        if (u != w) { pos[u] = pw; vert[pu] = w; pos[w] = pu; vert[pw] = u; }
        bin[du]++;
        deg[u]--;
      }
    }
  }
  if (order) *order = vert;
  return deg;
}


struct DAG {
  u32 n = 0;
  std::vector<u64> off;
  std::vector<u32> adj;
  u32 outdeg(u32 v) const { return static_cast<u32>(off[v + 1] - off[v]); }
  const u32* out(u32 v) const { return adj.data() + off[v]; }
};

inline DAG orient(const Graph& g, const std::vector<u32>& order) {
  std::vector<u32> rank(g.n);
  for (u32 i = 0; i < g.n; ++i) rank[order[i]] = i;
  DAG d;
  d.n = g.n;
  d.off.assign(g.n + 1, 0);
  for (u32 v = 0; v < g.n; ++v)
    for (u64 e = g.off[v]; e < g.off[v + 1]; ++e)
      if (rank[g.adj[e]] > rank[v]) d.off[v + 1]++;
  for (u32 v = 0; v < g.n; ++v) d.off[v + 1] += d.off[v];
  d.adj.resize(d.off[g.n]);
  std::vector<u64> pos(d.off.begin(), d.off.end() - 1);
  for (u32 v = 0; v < g.n; ++v)
    for (u64 e = g.off[v]; e < g.off[v + 1]; ++e)
      if (rank[g.adj[e]] > rank[v]) d.adj[pos[v]++] = g.adj[e];
  return d;
}
