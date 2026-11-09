#pragma once
#include <bit>
#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string_view>
#include <limits>
#include <cmath>

template<class State>
void export_certificate_if_requested(const State& s, int batch) {
  const char* dst=std::getenv("KCDS_CERT_FILE");
  if (!dst) return;
  const char* chosen=std::getenv("KCDS_CERT_BATCH");
  if (!chosen) throw std::runtime_error("KCDS_CERT_BATCH is required with KCDS_CERT_FILE");
  std::string_view token(chosen);
  int requested=0;
  const auto [end,error]=std::from_chars(token.data(),token.data()+token.size(),requested);
  if (error!=std::errc{} || end!=token.data()+token.size() || requested<1)
    throw std::runtime_error("invalid KCDS_CERT_BATCH");
  if (batch!=requested) return;
  static_assert(std::endian::native==std::endian::little && sizeof(double)==8
                && std::numeric_limits<double>::is_iec559);
  const std::filesystem::path target(dst), partial(target.string()+".partial");
  if (std::filesystem::exists(target) || std::filesystem::exists(partial))
    throw std::runtime_error("certificate path already exists");
  const uint32_t n=s.G.n, k=s.k, witness_size=static_cast<uint32_t>(s.S.size());
  if (s.S.empty() || s.S.size()>n || s.w.size()!=n || !std::isfinite(s.cS)
      || s.cS<0 || s.cS>0x1p53 || s.cS!=std::floor(s.cS))
    throw std::runtime_error("invalid certificate state");
  uint64_t arcs=0;
  for (const auto& nbr:s.G.adj) arcs+=nbr.size();
  const uint64_t m=arcs/2, count=static_cast<uint64_t>(s.cS);
  double raw_ub=s.tree.max(),epsilon=s.eps;
  std::ofstream f; f.exceptions(std::ios::badbit|std::ios::failbit);
  f.open(partial,std::ios::binary);
  auto one=[&](const auto& value){f.write(reinterpret_cast<const char*>(&value),sizeof(value));};
  f.write("ICERT001",8);one(n);one(k);one(m);one(witness_size);one(count);one(raw_ub);one(epsilon);
  uint64_t written=0;
  for (uint32_t u=0;u<n;++u) for (auto v:s.G.adj[u]) if (u<v) {
    uint32_t end_node=v;one(u);one(end_node);++written;
  }
  if (arcs%2 || written!=m) throw std::runtime_error("graph edge count mismatch");
  for (double w:s.w) one(w);
  for (uint32_t v:s.S) one(v);
  f.close();std::filesystem::rename(partial,target);
}
