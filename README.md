# ImpliCert

Implementation of *Implicit Certificates for Batch-Dynamic k-Clique Densest Subgraphs Without Materializing Cliques*.

<img src="figures/clique-deletion.png" width="650" alt="Deleting edge (0,1) from a five-vertex clique leaves two four-cliques. Repair closes the load gap.">

## Verify a deletion

With Python 3.10 or later, verify both sides of the deletion in the figure:

```sh
python verify_certificate.py data/before.json data/after.json
```

The checker enumerates the four-cliques independently. The witness density and maximum load agree at 1 before deletion and at 2/5 after deletion.

## Batch replay

On x64 Windows with Visual Studio C++ Build Tools and LLVM OpenMP:

```powershell
.\build.cmd
python run_strict_example.py
```

The example runs incremental maintenance, warm global repair and cold global repair on the same trace. Mean service time over three runs includes checkpoint export and independent verification after every update. [dblp and wiki-talk timings](data/strict-service.csv) · [Amazon, Enron and web-Google timings](data/strict-extension-runs.csv)

| Trace | ImpliCert (s) | Warm global PR (s) | Cold global PR (s) |
|---|---:|---:|---:|
| dblp | **15.06** | 17.70 | 18.91 |
| wiki-talk | 30.29 | 45.86 | **29.06** |
| Amazon | **16.74** | 24.75 | 23.26 |
| Enron | **14.89** | 25.96 | 21.53 |
| web-Google | **39.10** | 46.76 | 74.75 |

<details>
<summary>Numerical batch latency and repair frequency</summary>

<img src="figures/batch-outcomes.png" width="620" alt="Batch latency, cold checkpoint costs and triggered repairs across fifteen numerical workloads.">

</details>
