# Implicit Certificates for Batch-Dynamic k-Clique Densest Subgraphs Without Materializing Cliques

ImpliCert maintains a dense vertex set under edge insertions and deletions. Vertex potentials recover each affected clique’s shares. A strict four-clique query checks the exported witness and rounds the loads outward before returning.

## Run

64-bit Windows, Python 3.10+ and Visual Studio C++ Build Tools with LLVM OpenMP.

~~~powershell
.\build.cmd
python verify_certificate.py data/before.json data/after.json
python run_strict_example.py
~~~

[Data](data) includes the deletion checkpoints and all 378 verification records from the dblp and wiki-talk comparison.
