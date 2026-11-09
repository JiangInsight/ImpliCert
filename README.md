# ImpliCert

*Implicit Certificates for Batch-Dynamic k-Clique Densest Subgraphs Without Materializing Cliques*

Deleting one edge from a five-vertex clique destroys three of its five four-cliques. ImpliCert recovers their shares from vertex potentials. A local repair then brings the maximum load down to the witness density of 2/5.

<p align="center"><img src="figures/clique-deletion.png" width="660" alt="Deleting edge (0,1) opens a gap between density and maximum load. Updating the weights closes it."></p>

### Time through independent verification

Mean seconds for all 20 updates of a trace, over three runs. Every batch waits for export and an outward-rounded recount.

| Trace | ImpliCert | Warm global PR | Cold global PR |
|---|---:|---:|---:|
| dblp | **15.06** | 17.70 | 18.91 |
| wiki-talk | 30.29 | 45.86 | **29.06** |

<details>
<summary>Figure 2 · Numerical batch latency, cold checkpoint solves and repair triggers</summary>

<img src="figures/batch-outcomes.png" width="610" alt="Numerical batch latency and cold checkpoint costs across fifteen workloads, alongside the updates that triggered repair.">

</details>

Check the pictured deletion with Python 3.10+:

~~~sh
python verify_certificate.py data/before.json data/after.json
~~~

For the native strict replay on x64 Windows, run `build.cmd` with Visual Studio C++ Build Tools and LLVM OpenMP installed, then `python run_strict_example.py`.
