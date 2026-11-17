#define _CRT_SECURE_NO_WARNINGS
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <vector>

int main(int argc, char** argv) {
  if (argc < 3) { std::fprintf(stderr, "usage: prep_events <konect out file> <out.ebin>\n"); return 1; }
  FILE* f = std::fopen(argv[1], "rb");
  if (!f) { std::fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
  std::vector<uint32_t> U, V;
  std::vector<int8_t> OP;
  std::vector<int64_t> T;
  uint64_t lines = 0, loops = 0, bad = 0, zero_w = 0;
  uint32_t maxid = 0;
  bool has_t = true;
  static char buf[1 << 16];
  while (std::fgets(buf, sizeof buf, f)) {
    ++lines;
    const char* p = buf;
    while (*p == ' ' || *p == '\t') ++p;
    if (*p == '%' || *p == '\n' || *p == '\r' || *p == 0) continue;
    char* e;
    const unsigned long long a = std::strtoull(p, &e, 10);
    if (e == p) { ++bad; continue; }
    p = e;
    const unsigned long long b = std::strtoull(p, &e, 10);
    if (e == p) { ++bad; continue; }
    p = e;
    const double w = std::strtod(p, &e);
    if (e == p) { ++bad; continue; }
    p = e;
    const long long t = std::strtoll(p, &e, 10);
    if (e == p) has_t = false;
    if (a == 0 || b == 0 || a > 0xFFFFFFFEull || b > 0xFFFFFFFEull) { ++bad; continue; }
    if (a == b) { ++loops; continue; }
    if (w == 0) { ++zero_w; continue; }
    U.push_back((uint32_t)(a - 1));
    V.push_back((uint32_t)(b - 1));
    OP.push_back(w > 0 ? 1 : -1);
    T.push_back(e == p ? 0 : t);
    maxid = std::max<uint32_t>(maxid, (uint32_t)std::max(a, b));
  }
  std::fclose(f);
  const uint64_t M = U.size(), n = maxid;
  std::vector<uint64_t> idx(M);
  std::iota(idx.begin(), idx.end(), 0ull);
  if (has_t) std::stable_sort(idx.begin(), idx.end(), [&](uint64_t x, uint64_t y) { return T[x] < T[y]; });
  FILE* o = std::fopen(argv[2], "wb");
  if (!o) { std::fprintf(stderr, "cannot write %s\n", argv[2]); return 1; }
  std::fwrite(&n, 8, 1, o);
  std::fwrite(&M, 8, 1, o);
  auto put = [&](auto& src, size_t sz) {
    std::vector<char> tmp(M * sz);
    for (uint64_t i = 0; i < M; ++i) std::memcpy(&tmp[i * sz], &src[idx[i]], sz);
    std::fwrite(tmp.data(), 1, tmp.size(), o);
  };
  put(U, 4);
  put(V, 4);
  put(OP, 1);
  put(T, 8);
  std::fclose(o);
  uint64_t adds = 0;
  for (int8_t x : OP) adds += x > 0;
  std::printf("lines=%llu events=%llu adds=%llu removes=%llu self_loops=%llu zero_weight=%llu bad=%llu n=%llu "
              "timestamps=%d first_t=%lld last_t=%lld\n",
              (unsigned long long)lines, (unsigned long long)M, (unsigned long long)adds,
              (unsigned long long)(M - adds), (unsigned long long)loops, (unsigned long long)zero_w,
              (unsigned long long)bad, (unsigned long long)n, (int)has_t, M ? (long long)T[idx[0]] : 0LL,
              M ? (long long)T[idx[M - 1]] : 0LL);
  return 0;
}
