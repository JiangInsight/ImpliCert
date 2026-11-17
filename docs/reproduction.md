# Reproduction

The public core and experiment drivers target x64 Windows, matching the measured implementation. Install Python 3.12 or newer, Git, Visual Studio C++ Build Tools with LLVM OpenMP, and a WinLibs GCC distribution with C++20 and `__int128`. The paper used MSVC 19.51 and GCC 16.2. `--cxx` selects the compiler explicitly; none of the scripts needs the original research directory. The strict verifier is compiled with `/fp:strict`. Large instances need substantial RAM and disk; the original host had 63 GB RAM. Runs default to a 32 GB process limit, adjustable with `--memory-mb`.

```powershell
python -m pip install -r requirements.txt
python tools/build.py all
python tools/baselines.py all --cxx C:/toolchains/mingw64/bin/g++.exe
python tools/audit_reference.py
python tools/validate_release.py --baselines --real-data
```

`validate_release.py` is a correctness and integration check, not a replacement for the paper's timing experiments. It exercises the core contracts, three strict arms, triggered PR, all component ablations, orientation reuse, guarded diagnostics, fixed adversarial policies, historical versions, controlled repair, static and temporal conversion, and all seven external baseline executables. It downloads the complete Enron graph and checks the converted file against the original input hash before three short numerical/strict trajectories. CQ's functional checks use three updates at epsilon 0.01, 0.1 and 0.5; they do not replace full-workload timings. Each stage must return a successful status. All commands, output logs, input/binary hashes and acceptance records remain under `runs/`.

## Data

`benchmarks/datasets.json` identifies the 15 main graph inputs and the separate Friendster stress input. Sources are [SNAP](https://snap.stanford.edu/data/) and [KONECT](http://konect.cc/networks/link-dynamic-frwiki/). Downloads retain the upstream archive; data ownership and original terms remain with the providers.

```powershell
python tools/data.py fetch dblp wiki enron
python tools/data.py fetch friendster --native
python tools/data.py fetch frwiki --source C:/datasets/download.tsv.link-dynamic-frwiki.tar.bz2
```

KONECT availability varies; `--source` accepts the original archive obtained from the provider or an existing research copy. It is not a substitute graph. `fetch all` requests every complete input, including Friendster, and can require many gigabytes of downloads. `--directory` selects another dataset directory; pass it to `reproduce.py --data` as well. Stored data are ignored by Git. The repository's small CSV/JSON inputs are sufficient for the functional checks.

Static conversion removes self loops, merges undirected duplicates, relabels IDs in increasing order, and writes sorted CSR. Temporal conversion removes self loops, performs a stable timestamp sort, preserves repeated contacts, and writes the two endpoint arrays and timestamps. Events conversion preserves directed existence updates for the driver's undirected union state. It sorts by time only when every event has a timestamp. Native `prep_big` supplies the same CSR format without pandas for Friendster. SHA256 sidecars record every converted input. Five recorded input hashes are enforced automatically.

The formats use little-endian integers: CSR has `uint64 n,m`, `uint64 offsets[n+1]`, `uint32 adjacency[2m]`; temporal has `uint64 n,M`, `uint32 u[M],v[M]`, `int64 t[M]`; events additionally stores `int8 op[M]` before timestamps. Nested vertex samples use the original fixed permutation and seed 47, retain `round(p*n)` vertices, and relabel the induced graph in increasing original order.

## Main runs and strict service

```powershell
python reproduce.py list main
python reproduce.py run main --data datasets --output runs/main-guarded --snapshots
python reproduce.py run main --version v7 --data datasets --output runs/main-v7
python reproduce.py run main --version v7b --data datasets --output runs/main-v7b
python reproduce.py run strict --data datasets --output runs/strict
python reproduce.py run triggered --data datasets --output runs/triggered
python tools/summarize.py runs/strict
```

`list` prints the exact plan without running it. Select inputs with `--datasets dblp wiki`; `--case` filters case names. The default main plan uses 24 workers, k=4, epsilon=0.01, seed 1 and 20 batches. Static replay uses 1,000 edge updates and 10% holdout; wiki-talk uses 10,000 contacts and a 2M-contact window; StackOverflow uses a 10M-contact window. frwiki replays the final 200,000 observed updates. Checkpoint PR baselines are called at batches 10 and 20, with Orkut checked at 20. The same executable computes maintenance, cold checkpoints and warm checkpoints.

Strict service runs maintained/warm/cold arms on dblp, wiki-talk, Amazon, Enron and web-Google, with three repetitions and rotated order. It exports and verifies the initial state and every answer before publication. Its supported contract is k=4 and exact rational epsilon=1/100. wiki-talk uses 100,000 contacts per batch in this experiment. Triggered PR is a separate instrumented experiment: it retains incremental updates and hysteresis and runs local/warm/cold repair only when triggered. Global repair changes future potentials, so trigger counts need not remain identical. Its iteration budget is not the 300-pass checkpoint cap.

`--graph`, `--batch`, `--batches`, `--hold`, `--threads` and `--check` permit explicit diagnostic overrides; those values are saved in the run identity and do not reproduce the original timing setting. A completed run is reused only when command, input hash, executable hash, environment and limits agree. Failed or timed-out records are preserved; use a new output directory to retry. Runs execute serially. For timing measurements, use an otherwise idle host. `--affinity` accepts a Windows hexadecimal affinity mask; choose a performance core for single-worker baselines instead of assuming identical CPU numbering.

## External comparisons

The pinned sources and adaptations are documented in [baselines.md](baselines.md). Build-time downloads remain under ignored `external/`; licenses and provenance files are retained there. The build fails on a modified or differently pinned checkout.

```powershell
python tools/compare.py static runs/main-guarded/004_dblp_k4_replay_hybrid_s1_r1/snapshot_b20.txt --output runs/static-dblp
python reproduce.py run main --datasets dblp --hypergraph --output runs/export-dblp
python tools/compare.py dynamic runs/export-dblp/004_dblp_k4_replay_hybrid_s1_r1/hypergraph --output runs/dynamic-dblp-e001 --epsilon .01
python tools/compare.py dynamic runs/export-dblp/004_dblp_k4_replay_hybrid_s1_r1/hypergraph --output runs/dynamic-dblp-e01 --epsilon .1
python tools/summarize.py runs/dynamic-dblp-e01
```

Snapshot export and hyperedge export are separate runs. Hyperedge export deliberately bypasses the solver and is labelled `hypergraph-export-only`; it does not measure maintenance. Dynamic baselines receive exactly the generated clique insertion/deletion IDs and the query schedule, including zero-update batches. Enumeration is not charged to these baselines, as in the paper. CAP/CQ emit cumulative service timers; the summarizer differences consecutive queries. Udshp emits per-update service times and per-query extraction times. Serialization and preprocessing are not silently charged to the paper's service-time metric.

Static comparison passes the same snapshot to CCAS, KCCA, KClist++ and the matched PSCTL reconstruction. CCAS receives `1-1/(1+epsilon)`, KCCA/KClist++ receive `1+epsilon`. Maximum-clique preprocessing has its own record. KClist++ performs three runs by default and reports after every iteration; the reported solver time is cumulative iterations plus the final report. Whole-process wall time remains separately recorded. PSCTL supports the published k=4, 1% setting with a 300-pass cap. A cap or timeout is a recorded outcome, never a zero latency or an invented missing value.

## Ablations and additional experiments

All plans are generated by `tools/plans.py` and run through the same entry point:

| Suite | Coverage |
|---|---|
| `ablations` | full/noSAT/noABS/uniform/noSplit/noFallback, wiki/StackOverflow/skitter, three runs |
| `work` | separately instrumented clique visits, including fallback and prefix scans |
| `strategies` | hybrid, Newton and targeted repair on five inputs |
| `adversarial` | core deletions, planted regions, seeds and cold/warm fallback |
| `sensitivity` | k=3/5/6/7, epsilon=0.005/0.02/0.05, batch sizes, workers and seeds |
| `long` | wiki 200, StackOverflow/frwiki 1,000 batches, rerank and fixed orientation |
| `memory` | standard and compact state, including Friendster |
| `orientation` | rebuild versus reuse in the instrumented control |
| `policies` | a shared recorded adversarial trace, potential-box and PCG/q sweeps |
| `guarded` | guarded recount diagnostics on main/adversarial/high-k/small-epsilon cases |
| `friend` | ten controlled compact-state stress batches |
| `scaling` | nested 20/40/60/80/100% LiveJournal and Orkut vertex samples |

For example, `python reproduce.py run work --output runs/work` followed by `python tools/summarize.py runs/work` produces per-run means, nearest-rank p95, maxima, standard deviations and clique-work totals. `repeat-summary.json` summarizes variation of run means across seeds/repetitions. It does not treat ordinary `cliques_rep` as total work: full fallback and prefix counts require instrumented `work` records. Historical instrumentation predates the final atomic-claim fix; its identity remains distinct from `guarded` and `trace`.
