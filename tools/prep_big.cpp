#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <io.h>
#include <fcntl.h>
#include <omp.h>

using u32 = uint32_t;
using u64 = uint64_t;

int main(int argc, char** argv) {
  if (argc < 2) { std::fprintf(stderr, "usage: gzip -dc in.txt.gz | prep_big out.bin\n"); return 1; }
  _setmode(_fileno(stdin), _O_BINARY);
  std::vector<u32> eu, ev;
  eu.reserve(1u << 20);
  ev.reserve(1u << 20);
  const size_t BUF = 1 << 24;
  std::vector<char> buf(BUF + 1);
  size_t carry = 0;
  bool comment = false, in_num = false;
  u64 cur = 0, maxid = 0;
  int field = 0;
  u64 first = 0;
  u64 lines = 0;
  for (;;) {
    const size_t got = std::fread(buf.data(), 1, BUF, stdin);
    if (got == 0) break;
    for (size_t i = 0; i < got; ++i) {
      const char c = buf[i];
      if (c == '\n') {
        if (in_num) {
          if (field == 1) { eu.push_back((u32)first); ev.push_back((u32)cur); maxid = std::max(maxid, cur); }
        }
        in_num = false; comment = false; field = 0; cur = 0; ++lines;
        continue;
      }
      if (comment) continue;
      if (c == '#' && field == 0 && !in_num) { comment = true; continue; }
      if (c >= '0' && c <= '9') {
        cur = cur * 10 + (u64)(c - '0');
        in_num = true;
      } else if (in_num) {   
        if (field == 0) { first = cur; maxid = std::max(maxid, cur); field = 1; }
        else if (field == 1) { eu.push_back((u32)first); ev.push_back((u32)cur); maxid = std::max(maxid, cur); field = 2; }
        in_num = false; cur = 0;
      }
    }
    (void)carry;
  }
  if (in_num && field == 1) { eu.push_back((u32)first); ev.push_back((u32)cur); maxid = std::max(maxid, cur); }
  if (maxid >= 0xFFFFFFFFull) { std::fprintf(stderr, "ids exceed 32 bits\n"); return 1; }
  const u64 P = eu.size();
  std::fprintf(stderr, "read %llu pairs (%llu lines), max id %llu\n", (unsigned long long)P,
               (unsigned long long)lines, (unsigned long long)maxid);
   
  std::vector<u32> lab((size_t)maxid + 1, 0xFFFFFFFFu);
  for (u64 i = 0; i < P; ++i)
    if (eu[i] != ev[i]) { lab[eu[i]] = 0; lab[ev[i]] = 0; }
  u32 n = 0;
  for (u64 x = 0; x <= maxid; ++x) if (lab[x] == 0) lab[x] = n++;
  std::vector<u64> off((size_t)n + 1, 0);
  for (u64 i = 0; i < P; ++i)
    if (eu[i] != ev[i]) { off[lab[eu[i]] + 1]++; off[lab[ev[i]] + 1]++; }
  for (u32 v = 0; v < n; ++v) off[v + 1] += off[v];
  std::vector<u32> adj(off[n]);
  {
    std::vector<u64> pos(off.begin(), off.end() - 1);
    for (u64 i = 0; i < P; ++i) {
      if (eu[i] == ev[i]) continue;
      const u32 a = lab[eu[i]], b = lab[ev[i]];
      adj[pos[a]++] = b;
      adj[pos[b]++] = a;
    }
  }
  std::vector<u32>().swap(eu);
  std::vector<u32>().swap(ev);
  std::vector<u32>().swap(lab);
   
  std::vector<u64> newdeg((size_t)n, 0);
#pragma omp parallel for schedule(dynamic, 4096)
  for (long long v = 0; v < (long long)n; ++v) {
    u32* s = adj.data() + off[v];
    u32* e = adj.data() + off[v + 1];
    std::sort(s, e);
    newdeg[v] = (u64)(std::unique(s, e) - s);
  }
  u64 w = 0;
  std::vector<u64> noff((size_t)n + 1, 0);
  for (u32 v = 0; v < n; ++v) {
    const u64 s = off[v];
    noff[v] = w;
    for (u64 j = 0; j < newdeg[v]; ++j) adj[w++] = adj[s + j];
  }
  noff[n] = w;
  const u64 m = w / 2;
  FILE* f = std::fopen(argv[1], "wb");
  if (!f) { std::fprintf(stderr, "cannot write %s\n", argv[1]); return 1; }
  const u64 hdr[2] = {n, m};
  std::fwrite(hdr, sizeof(u64), 2, f);
  std::fwrite(noff.data(), sizeof(u64), (size_t)n + 1, f);
  std::fwrite(adj.data(), sizeof(u32), (size_t)w, f);
  std::fclose(f);
  u64 maxdeg = 0;
  for (u32 v = 0; v < n; ++v) maxdeg = std::max(maxdeg, noff[v + 1] - noff[v]);
  std::fprintf(stderr, "%s: n=%u m=%llu maxdeg=%llu\n", argv[1], n, (unsigned long long)m,
               (unsigned long long)maxdeg);
  return 0;
}
