# Runbook — Binary Alloy Cluster Optimization

## 1. Build

```bash
# VS 2022 Community default path:
"C:/Program Files/Microsoft Visual Studio/2022/Community/MSBuild/Current/Bin/MSBuild.exe" Cluster.sln -p:Configuration=Release
```

Output: `Release/Cluster.exe`

## 2. Configure

Edit `Cluster/data/config.txt`. Minimal example:

```ini
totalAtoms=23
numElementA=13
numElementB=10
elementA=Pt
elementB=Cu

potentialType=Gupta
potentialFile=data/gupta_PtCu.txt

algorithm=CDE
cde.populationSize=30
cde.maxGenerations=100
cde.useLocalSearch=true
cde.useMultiPopulation=true

numRuns=3
runAllCompositions=false

outputDirectory=results
verbose=true
```

### Potential ↔ parameter file mapping

| potentialType | parameter file |
|---|---|
| `Gupta` | `data/gupta_XxYy.txt` |
| `FinnisSinclair` | `data/fs_XxYy.txt` |
| `SuttonChen` | `data/sc_XxYy.txt` |

## 3. Run

### Single composition (from `Cluster/` directory)

```bash
cd Cluster
echo "" | ../Release/Cluster.exe -config data/config.txt
```

### Scan all compositions

Set `runAllCompositions=true` in config, then run. Each (nA, nB) pair is optimized sequentially.

### Command-line overrides

```
-atoms 38                  # override totalAtoms
-composition 19            # override numElementA
-elements Pt Co            # override elementA/B
-algorithm SaNSDE          # CDE | SaNSDE | PSO
-potential FinnisSinclair  # Gupta | FinnisSinclair | SuttonChen
-generations 300           # override maxGenerations
-runs 5                    # override numRuns
-allcompositions           # set runAllCompositions=true
```

## 4. Multithreading

CDE and SaNSDE support parallel population evaluation via `std::thread`.

### Enable

```ini
cde.useThreading=true      # for CDE
sansde.useThreading=true   # for SaNSDE
```

### How it works

- Each algorithm spawns 3 independent populations (sub-populations).
- When `useThreading=true`, populations run in parallel threads, one per population.
- A `std::mutex` protects the global-best individual when populations exchange or update it.
- `RandomGenerator` uses `thread_local` Mersenne Twister engines — no lock contention on RNG.

### Benchmark configs

Two pre-made configs for A/B testing:

| file | threading | purpose |
|---|---|---|
| `data/bench_thread.txt` | on | measure threaded throughput |
| `data/bench_no_thread.txt` | off | baseline comparison |

Run both and compare wall-clock time:

```bash
cd Cluster
time echo "" | ../Release/Cluster.exe -config data/bench_thread.txt
time echo "" | ../Release/Cluster.exe -config data/bench_no_thread.txt
```

### When to use threading

- **ON** — multi-core CPU, large population (≥30), FinnisSinclair/SuttonChen potentials (more expensive per evaluation).
- **OFF** — small systems (<20 atoms), small population, debugging (easier to trace single-threaded).

## 5. Output

Results land in `results/{PotentialType}/{CompositionName}/`:

```
results/Gupta/Pt13Cu10/
  1_best_energy_per_generation.txt   # convergence trace, run 1
  2_best_energy_per_generation.txt   # convergence trace, run 2
  3_best_energy_per_generation.txt   # convergence trace, run 3
  best.txt                           # best structure (Diamond format)
  best.xyz                           # best structure (XYZ format)
  energy.txt                         # per-generation population energies
  historical_best.txt                # cross-session best (Diamond)
  historical_best.xyz                # cross-session best (XYZ)
```

Each potential type gets its own subdirectory — SuttonChen / Gupta / FinnisSinclair results never overwrite each other.

## 6. Backup

```bash
bash backup.sh
# → ClusterOpt_20260525_HHMMSS.zip
```

Uses `git ls-files` to dynamically discover all source files. Automatically excludes `results/`, build artifacts, and VS caches. No hardcoded file lists — new source files are captured automatically.

## 7. Restore

```bash
unzip ClusterOpt_YYYYMMDD_HHMMSS.zip -d /path/to/restore/
```

## 8. Git Push

```bash
git add -A
git commit -m "your message"
git push origin master
```
