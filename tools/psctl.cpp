#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>
#include "graph.h"
#include "sct.h"
#include "pava.h"

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    constexpr unsigned k = 4;
    constexpr unsigned max_passes = 300;
    graph* g = build_graph(argv[1], k);
    auto start = std::chrono::steady_clock::now();
    std::vector<double> loads(g->n), counts(g->n);
    double best_density = 0;
    for (unsigned t = 1; t <= max_passes; ++t) {
        build_iter(g, k, loads.data(), false);
        auto bounds = pava_iter(g, loads.data(), counts.data(), k, t);
        double cliques = std::accumulate(counts.begin(), counts.end(), 0.0);
        double mass = std::accumulate(loads.begin(), loads.end(), 0.0) / t;
        if (std::abs(cliques - mass) > 1e-7 * std::max(1.0, cliques)) return 5;
        if (t == 1) std::printf("CLIQUES %.0f\n", cliques);
        best_density = std::max(best_density, bounds.second);
        double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        std::printf("PASS %u %.17g %.17g %.9f\n", t, best_density, bounds.first, seconds);
        std::fflush(stdout);
        if (!std::isfinite(bounds.first) || !std::isfinite(bounds.second) ||
            best_density < 0 || bounds.first + 1e-8 < best_density) return 3;
        if (best_density > 0 && bounds.first <= 1.01 * best_density) return 0;
    }
    return 4;
}
