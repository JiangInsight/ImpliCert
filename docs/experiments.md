# Experiment-to-artifact mapping

The archive and rerunnable versions are deliberately identified separately. `results/reference/` contains historical records, not outputs of a new timing campaign. `tools/audit_reference.py` verifies every main log against the recorded SHA256 and produces a compact summary. `benchmarks/source_manifest.json` records both research-source hashes and the packaged source hashes after removing comments and adjusting local includes. Current public build records include source/header and executable hashes.

| Paper evidence | Source/entry point | Data and recorded output |
|---|---|---|
| Main numerical maintenance and checkpoint PR | archived `versions/v7/implicert-v7.exe`; `main --version v7`, corrected `guarded` and source-buildable `v7b` reruns | 15 input/settings rows in `results/reference/version_map.csv`; unchanged raw logs in `results/reference/main/` |
| Original versus guarded comparison | `versions/guarded_audit`, `guarded` suite; `v7b` for the matched original arm | `results/reference/round2/` and version map, with matching batch signatures |
| Triggered local/warm/cold global PR | `versions/triggered`, `triggered` suite | wiki-talk, skitter, LiveJournal; archived round2 policy records |
| Corrected per-answer strict service | `src/trace.cpp`, `src/barrier.hpp`, native outward checker and exact rational Python decision; `strict` suite | five complete graph traces, three repetitions, archived `data/strict-service.csv`, `data/strict-extension-runs.csv`, and reference `icde_review`/`strict_extension` |
| SAT/absorber/split/target/fallback and strategy ablations | `ablations`, `strategies`, separately instrumented `work` | original numerical records versus rerun outputs remain separate; work counters include fallback |
| Orientation and potential-box/PCG policy diagnostics | `versions/control`, `versions/policy`; `orientation`, `policies` | `results/reference/round2/`, `results/reference/review3/policies/` |
| Guarded count ranges, saturation and adversarial updates | `versions/guarded_audit`, contract/absorber microtests, `guarded` | archived guarded and saturation CSVs; current release tests |
| Friendster controlled repair | `versions/friend`, `friend`; compact mode in `memory` | `results/reference/review3/friend/`; full Friendster must be downloaded separately |
| CCAS/KCCA/KClist++ and matched PSCTL | `tools/baselines.py`, `tools/compare.py static` | solver snapshots from `--snapshots`; archived matched PSCTL CSVs |
| CAP/Udshp/CQ dynamic comparison | same pinned build, `tools/compare.py dynamic` | exported `.updates` and `.queries.tsv`, epsilon 0.01/0.1, one worker |
| k, epsilon, size, workers, update size, long traces and seeds | `sensitivity`, `scaling`, `memory`, `long`, `adversarial` | explicit plans and recorded run identities |

## Version boundary

The original main-table v7 binary is distributed in `versions/v7/implicert-v7.exe` with its SHA256 identity and runs through `main --version v7`, but its exact corresponding source snapshot was not recovered. The research directory, Overleaf history, source/supplement packages and reproduce archives were searched. The oldest recoverable later source is v7b, described by the original queue as v7 plus a stagnation rule. v7c adds the deletion-fraction workload switch. Removing that rule and naming the result v7 would not establish binary/source correspondence, so this artifact does not do so. The unchanged original logs and binary hash remain available in the version map. The runnable v7b plan is a historical-version rerun, not a claim of exact original-v7 reproduction. The original executable permits binary-level reruns; source-level reproduction of v7 remains unavailable.

`src/implicert.cpp` is the guarded release with atomic absorber ownership. The historical `triggered`, `control`, `profile`, `policy` and guarded-audit diagnostic snapshots preserve their original identities and can predate that final fix. The `friend` stress driver applies the atomic ownership helper. The current strict `trace` uses the corrected release plus a synchronous independent publication barrier. No experiment from the earlier snapshots is silently relabelled as a final-release measurement.

Binary64 maintenance records contain numerical closure/recounts. They do not themselves provide a rigorous floating-point upper bound. The `check` line's historical name `exact_ub` means a fresh binary64 recount, not directed-rounding certification. Only successfully accepted strict checkpoints are labelled independently verified. Failed states, timeouts, potential/count-range rejection and incomplete output are retained as failures, never folded into successful timing statistics.

## Timing definitions

The numerical main table averages all batch update-plus-repair times, including batches that do not trigger a repair. Its cold speedup divides the mean cold checkpoint time by that all-batch mean. Checkpoint warm/cold PR has its own cap and stopping rule and differs from triggered global repair. Strict service includes graph/witness/weight export and independent verification for every answer. Dynamic materialized-clique baselines are not charged for graph-to-hyperedge enumeration. Static solver timers exclude input/core preparation and maximum-clique search where the original comparison excluded them. Whole-process wall time is stored independently for auditing.

The provided full plans enable new measurements; they are not newly completed runs. The release validation is intentionally small, with a complete Enron input but only three batches and one worker, plus deterministic fixtures. Friendster preprocessing/repair, the large temporal traces, 24-worker timing, three-run full strict results and one-hour static limits must be rerun explicitly on the intended measurement host. No GPU is required.
