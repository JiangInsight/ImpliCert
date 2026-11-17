#pragma once
#include <numeric>

namespace guard_audit {
constexpr u64 exact_limit = u64(1) << 53;
constexpr u64 too_large = exact_limit + 1;
struct alignas(64) Box {
  u64 writes = 0, lower = 0, upper = 0, clipped = 0, nan_proposals = 0, inf_proposals = 0, invalid_stored = 0;
};
static std::vector<Box> boxes;
static u64 largest_graph_bound = 0, largest_witness = 0;

inline void record_potential(double proposal, double center, double result, double radius) {
  auto& b = boxes[omp_get_thread_num()];
  ++b.writes;
  b.lower += result == center-radius;
  b.upper += result == center+radius;
  b.clipped += proposal<center-radius || proposal>center+radius;
  b.nan_proposals += std::isnan(proposal);
  b.inf_proposals += std::isinf(proposal);
  b.invalid_stored += !std::isfinite(result);
}
inline void report(int batch) {
  Box b;
  for (const auto& x : boxes) {
    b.writes+=x.writes; b.lower+=x.lower; b.upper+=x.upper;
    b.clipped+=x.clipped; b.nan_proposals+=x.nan_proposals;
    b.inf_proposals+=x.inf_proposals; b.invalid_stored+=x.invalid_stored;
  }
  std::printf("audit_state batch=%d writes=%llu box_lo=%llu box_hi=%llu clipped=%llu nan_proposals=%llu inf_proposals=%llu invalid_stored=%llu max_graph_clique_bound=%llu max_witness=%llu\n",
              batch,b.writes,b.lower,b.upper,b.clipped,b.nan_proposals,b.inf_proposals,b.invalid_stored,largest_graph_bound,largest_witness);
  if (b.invalid_stored) throw std::runtime_error("nonfinite_stored_potential");
}
template<class State> void stored_box(const State& state,int batch,double radius) {
  u64 at_lower=0,at_upper=0,bad=0;
  for (double v:state.lam) {
    at_lower+=v==state.lam0-radius;
    at_upper+=v==state.lam0+radius;
    bad+=!std::isfinite(v) || v<state.lam0-radius || v>state.lam0+radius;
  }
  std::printf("stored_box batch=%d vertices=%zu lower=%llu upper=%llu invalid=%llu\n",
              batch,state.lam.size(),at_lower,at_upper,bad);
  if (bad) throw std::runtime_error("invalid_stored_potential");
}
inline double exact_count(u64 value) {
  if (value>exact_limit) throw std::runtime_error("witness_integer_range");
  largest_witness=std::max(largest_witness,value);
  return static_cast<double>(value);
}
inline void observe_witness(double value) {
  if (!std::isfinite(value) || !(value>=0 && value<=double(exact_limit)) || std::trunc(value)!=value)
    throw std::runtime_error("witness_integer_range");
  largest_witness=std::max(largest_witness,static_cast<u64>(value));
}
inline u64 choose_bound(u64 n,u64 k) {
  if (n<k) return 0;
  k=std::min(k,n-k);
  u64 v=1;
  for (u64 j=1; j<=k; ++j) {
    u64 a=n-j+1, b=j, g=std::gcd(a,b); a/=g; b/=g;
    v/=b;  
    if (v>exact_limit/a) return too_large;
    v*=a;
  }
  return v;
}

template<class Out> void graph_range(const Out& out,int k) {
  std::vector<u64> table;
  u64 total=0;
  for (const auto& neighbours:out) {
    const size_t d=neighbours.size();
    while (d>=table.size()) table.push_back(choose_bound(table.size(),k-1));
    if (table[d]>exact_limit-total) throw std::runtime_error("graph_clique_bound_exceeds_exact_range");
    total+=table[d];
  }
  largest_graph_bound=std::max(total,largest_graph_bound);
}
}
