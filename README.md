# Binary Alloy Cluster Global Optimization

Evolutionary algorithms (CDE, SaNSDE) for optimizing binary alloy nanoparticle structures with Gupta, Finnis-Sinclair, and Sutton-Chen interatomic potentials. L-BFGS local refinement included.

## Quick Start

Open `Cluster.sln` in Visual Studio 2019+ and build as Release.

### Configure

Edit `Cluster/data/config.txt`:

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

### Run

```bash
cd Cluster
echo "" | ../Release/Cluster.exe -config data/config.txt
```

Command-line overrides: `-atoms`, `-composition`, `-elements`, `-algorithm`, `-potential`, `-generations`, `-runs`, `-allcompositions`.

## Multithreading

CDE and SaNSDE support parallel population evaluation. Three sub-populations run in parallel threads, with `std::mutex` guarding the global best. `RandomGenerator` uses `thread_local` engines — no RNG lock contention.

Enable in config:
```ini
cde.useThreading=true      # or sansde.useThreading=true
```

## Potentials

| Type | Config Key | Parameters | Best For |
|------|-----------|------------|----------|
| **Gupta** | `Gupta` | r0, A, xi, p, q | FCC/noble metals (Pt, Au, Pd, Cu, Ag, Ni, Co) |
| **Finnis-Sinclair** | `FinnisSinclair` | cutoff, c0, c1, c2, d, c, beta | BCC transition metals (Fe, Cr, Mo, W) |
| **Sutton-Chen** | `SuttonChen` | epsilon, c, a, n, m | FCC metals alternative |

Parameter files in `Cluster/data/` follow a 3-line format (A-A, B-B, A-B). Available pairs: AgPd, AuAg, CuAu, NiAl, PtCo, PtCu, PtNi, PtPd, plus pure elements (CuCu, FeFe, PtPt for FS).

## Output Structure

Results are nested by potential type to prevent overwriting:

```
results/
  Gupta/Pt13Cu10/
    best.xyz                          # Best structure (XYZ)
    best.txt                          # Best structure (Diamond)
    historical_best.txt               # Cross-session best
    1_best_energy_per_generation.txt  # Run 1 convergence
    energy.txt                        # Population energies
  SuttonChen/...
  FinnisSinclair/...
```

## Parameter File Format

One file per element pair. Line 1: element names. Lines 2–4: A-A, B-B, A-B parameters (CSV).

Gupta example (`data/gupta_PtCu.txt`):
```
Pt Cu
2.7747, 0.2975, 2.695, 10.612, 4.004
2.556, 0.0894, 1.224, 10.96, 2.278
2.66535, 0.163, 1.816, 10.786, 3.141
```

If A-B parameters are unavailable, use mixing rules: `r0_AB = (r0_AA + r0_BB) / 2`, `A_AB = sqrt(A_AA * A_BB)`, etc.

## Visualization

Open `best.xyz` with [VESTA](https://jp-minerals.org/vesta/), [Ovito](https://www.ovito.org/), or [ASE](https://wiki.fysik.dtu.dk/ase/).

## Backup

```bash
bash backup.sh   # creates ClusterOpt_YYYYMMDD_HHMMSS.zip
```

Uses `git ls-files` for dynamic file discovery — new source files are automatically captured. Excludes build artifacts.

## Configuration Reference

| Key | Default | Description |
|-----|---------|-------------|
| `totalAtoms` | 38 | Total atoms in cluster |
| `elementA` / `elementB` | Pt / Co | Element symbols |
| `potentialType` | Gupta | Gupta, FinnisSinclair, or SuttonChen |
| `potentialFile` | data/gupta_PtCo.txt | Path to parameter file |
| `algorithm` | CDE | CDE, SaNSDE, or PSO |
| `cde.populationSize` | 40 | Sub-population size |
| `cde.maxGenerations` | 200 | Max generations |
| `cde.exchangeInterval` | 20 | Inter-population exchange frequency |
| `cde.useLocalSearch` | true | Enable L-BFGS refinement |
| `cde.useMultiPopulation` | true | Use 3 sub-populations |
| `cde.useThreading` | false | Parallel population evaluation |
| `sansde.populationSize` | 40 | Population size |
| `sansde.maxGenerations` | 200 | Max generations |
| `sansde.useLocalSearch` | true | Enable L-BFGS refinement |
| `sansde.useThreading` | false | Parallel population evaluation |
| `numRuns` | 1 | Independent runs per composition |
| `runAllCompositions` | false | Scan all (nA, nB) pairs |
| `outputDirectory` | results | Output path |
| `verbose` | true | Detailed console output |

## Version History

### v3.0 (Current)
- Multithreading support (parallel populations with thread-local RNG)
- Results auto-nested by potential type to prevent overwriting

### v2.0
- Added Finnis-Sinclair and Sutton-Chen potentials
- Runtime potential selection via config/command-line

### v1.0
- Initial release: Gupta potential, CDE and SaNSDE algorithms

## References

1. **Gupta**: Cleri & Rosato (1993), *Phys. Rev. B*, 48, 22-33
2. **Finnis-Sinclair**: Finnis & Sinclair (1984), *Phil. Mag. A*, 50, 45-55
3. **Sutton-Chen**: Sutton & Chen (1990), *Phil. Mag. Lett.*, 61, 139-146

## License

MIT
