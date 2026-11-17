# Baseline provenance

| Executable | Public source | Fixed revision |
|---|---|---|
| CCAS | https://github.com/forxenn/ccas | `c1a8032ec4f415958e78e1d9c3a49207dbda9cf4` |
| KCCA | https://github.com/forxenn/KCCA | `de420e1badc878989443bec75a95eca8d6e277b1` |
| KClist++ | https://github.com/btsun/kClistpp | `2704e72c8fa26a4be1adaa4ee8128456811ade99` |
| CAP, Udshp, CQ | https://github.com/chenjb1997/cap-artifact | `6bd103bf86e7c40b32f122efc6b8bb599b868e70` |
| PSCTL reconstruction | CCAS source above | same CCAS revision |

The external checkouts retain their original licenses and provenance. KClist++ is MIT licensed, copyright Bintao Sun. CCAS/KCCA do not provide a blanket license in the pinned checkout. The CAP artifact's `baselines/PROVENANCE.md` identifies its Udshp source as the Bera et al. author artifact and notes that no license file was supplied for those files. This repository does not assign new rights to third-party code. Sources are fetched from upstream rather than republished under an assumed license.

CQ is the CAP artifact's reconstruction of Chekuri–Quanrud, not the original authors' implementation. In particular, its explicit scheduling of changed endpoints after deletions must not be presented as measurement of the original worst-case guarantee. Its default copying rule can be very expensive, including on small inputs. The runner preserves timeouts; the release functional checks use three updates at epsilon=0.01/0.1/0.5; these checks do not claim to reproduce full-workload timing results.

CCAS and KCCA use their documented GCC optimization flags; CCAS additionally includes `cassert` for current GCC compatibility. KClist++ uses MSVC's LLVM OpenMP. The KClist++ adapter changes three OpenMP loop counters for MSVC and adds the paper's per-iteration quality reporting and separate iteration/report timers. `benchmarks/kclist_adapter.json` records the exact transformation and source/result hashes. It does not change clique enumeration or Frank–Wolfe updates. Relative output names avoid the upstream fixed-size filename buffer.

The CAP and CQ adapters write their already-stopped cumulative update/query timers at scheduled queries; algorithm updates are unchanged. Udshp uses the artifact's supplied service timer. zlib 1.3.1 is fetched from the [official release](https://github.com/madler/zlib/releases/tag/v1.3.1), checked against SHA256 `9a93b2b7dfdac77ceba5a558a580e74667dd6fede4585b91eefb60f03b72df23`, and compiled with the same GCC distribution. Its upstream license is preserved in the extracted source.

PSCTL is reconstructed from the exposed `build_iter(...,false)` and `pava_iter` routines in CCAS. It is not a separately released PSCTL executable. The matched reconstruction restores callback stack ordering, uses a consistent load/ID tie order in prefix charging, releases traversal scratch, initializes the binomial table once, and evaluates large binomials without indexing past its fixed table. The public driver checks clique-count/mass consistency and reports its own bounds at every pass. The earlier unmatched counting experiment is not a valid baseline result. `results/reference/review3/static/` contains the matched measurements; see the paper's distinction between standalone certification and externally matched quality.
