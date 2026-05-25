# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build

### VS IDE

- Open `Cluster.sln` in Visual Studio 2019+ (C++17, x64, Release)
- Build → Build Solution (`Ctrl+Shift+B`), output at `x64/Release/BinaryAlloyOptimizer.exe`

### Command line (MSBuild)

Replace `<VS_PATH>` with your VS install path (e.g. `C:/Program Files/Microsoft Visual Studio/2022/Community`):

```bash
"<VS_PATH>/MSBuild/Current/Bin/MSBuild.exe" Cluster.sln -p:Configuration=Release
```

Output: `Release/Cluster.exe` (or `x64/Release/` if platform explicitly set).

## Architecture

### Strategy pattern for interatomic potentials

`PotentialBase` (`include/PotentialBase.h`) defines the interface. Three implementations:

- **GuptaPotential** — `r0, A, xi, p, q` parameters, best for FCC/noble metals (Pt, Au, Pd, Cu, Ag, Ni, Co)
- **FinnisSinclairPotential** — `cutoff, c0..c2, d, c, beta` parameters, for BCC transition metals (Fe, Cr, Mo, W)
- **SuttonChenPotential** — `epsilon, c, a, n, m` parameters, alternative for FCC metals

Each potential reads 3 lines from a parameter file (A-A, B-B, A-B interactions). The base class declares pure virtual `calculateEnergy()` and `calculateEnergyWithForces()` — the latter is what L-BFGS calls.

### Three optimization algorithms with similar structure

All operate on `BinaryAlloyCluster` (3N coordinates + N type integers):

- **CDE** — Collaborative Differential Evolution. Multi-population scheme (`CDE_Population`) with 5 mutation strategies (`RAND1`, `RAND2`, `BEST1`, `BEST2`, `RAND_TO_BEST1`). Populations exchange best individuals periodically.
- **SaNSDE** — Self-adaptive DE with Neighborhood Search. 3 parallel populations. Self-adapts F, CR, and mutation strategy selection via a success memory deque. Supports `neighborhoodSearch()` for local refinement.
- **PSO** — Particle Swarm Optimization. Velocity-based, inertia weight decays linearly. More experimental, less used.

### L-BFGS local optimization

`NELbfgs` (`include/NELbfgs.h`) wraps the C `liblbfgs` library. It provides static C callbacks (`evaluate`, `progress`) that bridge to the C++ potential interface. `LocalOptimizer` (`include/LocalOptimizer.h`) is a thin facade. Algorithms inject local search after crossover/selection steps, controlled by `localSearchFrequency`.

### Configuration flow

1. `Configuration::loadFromFile()` parses `data/config.txt` (ini-style `key=value`)
2. Command-line args in `main.cpp` override config values
3. `Configuration::SystemConfig` bundles all settings — algorithm params, potential selection, run control
4. `runAllCompositions=true` iterates all (nA, nB) pairs for a given `totalAtoms`

### Output structure

Results are nested by potential type to avoid overwriting (e.g. `results/Gupta/Pt7Cu6/`):

`ResultManager` handles file I/O per run:
- `best.xyz` / `best.txt` — best structure for a composition
- `*_best_energy_per_generation.txt` — convergence trace per run
- `historical_best.*` — cross-run best, persisted and loaded between sessions
- `all_compositions_summary.txt` — aggregate when scanning compositions
- XYZ format for visualization (VESTA/Ovito), Diamond format for internal use

### Key data classes

- **`BinaryAlloyCluster`** — flat `vector<double>` coordinates (3N) + `vector<int>` atom types. Knows its element symbols. Has `randomInitialize(boxSize)` for random starting points.
- **`RandomGenerator`** — thread-local Mersenne Twister with uniform, normal, cauchy distributions

### Thread safety

`CDE` and `SaNSDE` use `std::mutex` for updating the global best individual across parallel populations. `RandomGenerator` uses `thread_local` engines.
