#define main guarded_solver_cli
#include "implicert.cpp"
#undef main

#include <bit>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace {
constexpr u32 vertex_count = 5;
constexpr int clique_size = 4;
constexpr int epsilon_denominator = 100;

Solver initialize() {
  std::vector<std::vector<u32>> adjacency(vertex_count);
  for (u32 u = 0; u < vertex_count; ++u)
    for (u32 v = 0; v < vertex_count; ++v)
      if (u != v) adjacency[u].push_back(v);
  Graph graph = csr_from_lists(adjacency);

  Solver state;
  state.k = clique_size;
  state.eps = 1.0 / epsilon_denominator;
  state.G.init(graph);
  state.lam.assign(vertex_count, 0.0);
  state.w.assign(vertex_count, 1.0);
  state.r.assign(vertex_count, 0.0);
  state.inS.assign(vertex_count, 0);
  state.isChanged.assign(vertex_count, 0);
  state.xidx.assign(vertex_count, -1);
  state.hpos.assign(vertex_count, 0);
  state.hstamp.assign(vertex_count, 0);
  state.inR.assign(vertex_count, 0);
  state.rstamp.assign(vertex_count, 0);
  state.strategy = 2;
  state.two_sided = false;
  state.hybrid = state.use_newton = true;
  state.ts.resize(1);
  state.ts.front().init(vertex_count, clique_size);
  state.hcnt.assign(1, std::vector<double>(Solver::kNB, 0.0));
  state.hused.resize(1);
  state.absorb.assign(vertex_count, 0.0);
  state.absflag.assign(vertex_count, 0);
  state.abslist.resize(1);
  state.vwork.assign(vertex_count, 0.0f);

  double lb = 0, prep = 0, iter = 0;
  u64 witness_cliques = 0;
  auto status = global_pr(graph, state.k, 1 + state.eps / 2, state.lam, state.lam0,
                          state.w, state.r, lb, state.S, kSolveMaxIt, &prep, &iter,
                          0, 0, &witness_cliques, true);
  require_solve(status, "example initialization");
  for (u32 v : state.S) state.inS[v] = 1;
  state.cS = static_cast<double>(witness_cliques);
  state.tree.build(state.r);
  return state;
}

void save(const Solver& s, const char* path) {
  require_certificate(s.tree.max(), s.rhoS(), 1 + s.eps, s.S.size(), "example checkpoint");
  std::ofstream out(path);
  out.exceptions(std::ios::badbit | std::ios::failbit);
  out << "{\n  \"schema\": \"implicit-kclique-certificate-v1\",\n  \"n\": " << s.G.n
      << ", \"k\": " << s.k << ", \"epsilon\": \"1/" << epsilon_denominator << "\",\n  \"edges\": [";
  bool first = true;
  for (u32 u = 0; u < s.G.n; ++u) {
    for (u32 v : s.G.adj[u]) {
      if (u >= v) continue;
      if (!first) out << ", ";
      first = false;
      out << '[' << u << ", " << v << ']';
    }
  }
  out << "],\n  \"witness\": [";
  for (size_t i = 0; i < s.S.size(); ++i) {
    if (i) out << ", ";
    out << s.S[i];
  }
  out << "],\n  \"weight_bits\": [";
  for (size_t i = 0; i < s.w.size(); ++i) {
    if (i) out << ", ";
    out << '"' << std::hex << std::setfill('0') << std::setw(16)
        << std::bit_cast<uint64_t>(s.w[i]) << std::dec << '"';
  }


  constexpr int64_t scale = 1000000000;
  const auto ub_grid = static_cast<int64_t>(std::ceil(s.tree.max() * scale));
  out << "],\n  \"witness_count\": " << static_cast<u64>(s.cS)
      << ",\n  \"lower\": \"" << static_cast<u64>(s.cS) << '/' << s.S.size()
      << "\", \"upper\": \"" << ub_grid << '/' << scale << "\"\n}\n";
}
}

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "usage: export_checkpoint before.json after.json\n";
    return 2;
  }
  try {
    omp_set_num_threads(1);
    Solver s = initialize();
    save(s, argv[1]);
    s.apply_batch({{0, 1}}, {});
    const int rounds = s.repair();
    if (rounds < 0) throw std::runtime_error("example repair did not finish");
    save(s, argv[2]);
    std::cout << "Deleted edge (0,1); repair rounds: " << rounds << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "Export failed: " << error.what() << '\n';
    return 1;
  }
}
