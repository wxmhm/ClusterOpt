# Binary Alloy Cluster Global Optimization

A C++ implementation of evolutionary algorithms for optimizing binary alloy nanoparticle structures with **multiple interatomic potential support**.

## ✨ Key Features

- **Multiple Interatomic Potentials:**
  - **Gupta Potential**: Optimized for metallic clusters (Au, Pt, Pd, Cu, Ag, Ni, Co)
  - **Finnis-Sinclair Potential**: Ideal for transition metals (Fe, Cr, Mo, W)
  - **Sutton-Chen Potential**: Alternative for FCC metals

- **Three Optimization Algorithms:**
  - **CDE**: Improved Differential Evolution with multiple mutation strategies
  - **SaNSDE**: Self-adaptive Differential Evolution with Neighborhood Search
  - **PSO**: Particle Swarm Optimization

- **Advanced Features:**
  - L-BFGS local optimization for structure refinement
  - Flexible configuration system
  - Multi-run capability with statistical analysis
  - Historical best tracking across runs
  - Composition screening mode

## 📋 Requirements

- Windows 10 or later
- Visual Studio 2019 or later (with C++17 support)
- liblbfgs library (included in lib/ folder)

## 🚀 Quick Start

### 1. Build the Project

1. Open `BinaryAlloyOptimizer.sln` in Visual Studio 2019+
2. Select **Release** configuration and **x64** platform
3. Build → Build Solution (or press `Ctrl+Shift+B`)
4. Executable will be in `x64/Release/BinaryAlloyOptimizer.exe`

### 2. Prepare Configuration File

Edit `data/config.txt`:

```ini
# Basic settings
totalAtoms=38
numElementA=19
numElementB=19
elementA=Pt
elementB=Co

# Select potential type
potentialType=Gupta
potentialFile=data/gupta_PtCo.txt

# Algorithm settings
algorithm=CDE
CDE.populationSize=40
CDE.maxGenerations=200
CDE.useLocalSearch=true

# Run control
numRuns=10
runAllCompositions=false
```

### 3. Prepare Potential Parameters

#### For Gupta Potential

Create/edit `data/gupta_PtCo.txt`:

```
Pt Co
2.6975,0.2975,1.64249,10.612,4.004
2.4985,0.1765,1.18997,10.960,3.814
2.5980,0.2370,1.41623,10.786,3.909
```

Format: Three lines with `r0, A, xi, p, q` for:
- Line 1: Pt-Pt interactions
- Line 2: Co-Co interactions
- Line 3: Pt-Co interactions

#### For Finnis-Sinclair Potential

Create `data/fs_FeCr.txt`:

```
Fe Cr
5.0,1.0,-0.5,0.1,4.5,1.0
5.0,1.0,-0.5,0.1,4.5,1.0
5.0,1.0,-0.5,0.1,4.5,1.0
```

Format: `cutoff, c0, c1, c2, d, c` for A-A, B-B, A-B interactions

#### For Sutton-Chen Potential

Create `data/sc_AuCu.txt`:

```
Au Cu
1.0,39.432,3.61,10.0,8.0
1.0,39.432,3.61,10.0,8.0
1.0,39.432,3.61,10.0,8.0
```

Format: `epsilon, c, a, n, m` for A-A, B-B, A-B interactions

### 4. Run the Program

#### Method 1: Double-click (Recommended)

1. Double-click `BinaryAlloyOptimizer.exe`
2. It automatically loads `data/config.txt`
3. Results saved to `results/` folder

#### Method 2: Command Line with Options

```bash
# Basic usage
BinaryAlloyOptimizer.exe -atoms 38 -composition 19

# Specify potential
BinaryAlloyOptimizer.exe -potential Gupta -atoms 38 -composition 19

# Full control
BinaryAlloyOptimizer.exe -potential FinnisSinclair ^
    -atoms 38 -composition 19 ^
    -elements Fe Cr ^
    -runs 10 -generations 200
```

**Command Line Options:**
- `-config <file>` - Configuration file (default: data/config.txt)
- `-potential <type>` - Gupta, FinnisSinclair, or SuttonChen
- `-potentialFile <file>` - Path to potential parameter file
- `-atoms <n>` - Total number of atoms
- `-composition <nA>` - Number of element A atoms
- `-elements <A> <B>` - Element symbols
- `-algorithm <type>` - CDE, SaNSDE, or PSO
- `-generations <n>` - Maximum generations
- `-runs <n>` - Number of independent runs
- `-allcompositions` - Run all possible compositions
- `-help` - Show help message

## 📊 Example Cases

### Case 1: Pt-Co Core-Shell (Gupta Potential)

**Configuration** (`data/config.txt`):
```ini
totalAtoms=38
numElementA=19
numElementB=19
elementA=Pt
elementB=Co
potentialType=Gupta
potentialFile=data/gupta_PtCo.txt
algorithm=CDE
CDE.populationSize=40
CDE.maxGenerations=200
CDE.useLocalSearch=true
numRuns=50
```

**Expected Output**: Core-shell structure with Pt shell, Co core

### Case 2: Fe-Cr BCC Cluster (Finnis-Sinclair)

```ini
totalAtoms=19
numElementA=13
numElementB=6
elementA=Fe
elementB=Cr
potentialType=FinnisSinclair
potentialFile=data/fs_FeCr.txt
algorithm=CDE
CDE.maxGenerations=300
numRuns=10
```

### Case 3: Au-Cu Alloy (Sutton-Chen)

```ini
totalAtoms=55
numElementA=38
numElementB=17
elementA=Au
elementB=Cu
potentialType=SuttonChen
potentialFile=data/sc_AuCu.txt
algorithm=SaNSDE
sansde.maxGenerations=500
numRuns=20
```

### Case 4: Composition Screening

```ini
totalAtoms=38
elementA=Pt
elementB=Ni
potentialType=Gupta
potentialFile=data/gupta_PtNi.txt
runAllCompositions=true
numRuns=5
```

Optimizes all compositions: Pt38, Pt37Ni1, ..., Ni38

### Case 5: Comparing Potentials

```bash
# Test Gupta
program.exe -potential Gupta -atoms 38 -composition 19 ^
    -elements Au Cu -potentialFile data/gupta_AuCu.txt

# Test Sutton-Chen
program.exe -potential SuttonChen -atoms 38 -composition 19 ^
    -elements Au Cu -potentialFile data/sc_AuCu.txt
```

## 🗂️ Project Structure

```
BinaryAlloyOptimizer/
├── data/
│   ├── config.txt                    # Main configuration
│   ├── gupta_PtCo.txt               # Gupta parameters
│   ├── fs_FeCr.txt                  # Finnis-Sinclair (optional)
│   ├── sc_AuCu.txt                  # Sutton-Chen (optional)
│   └── initial_structures/          # Initial structures (optional)
│
├── include/                          # Header files
│   ├── PotentialBase.h              # ← NEW: Base class
│   ├── GuptaPotential.h
│   ├── FinnisSinclairPotential.h    # ← NEW
│   ├── SuttonChenPotential.h        # ← NEW
│   ├── Configuration.h               # ← Modified
│   └── ...
│
├── src/                              # Source files
│   ├── main.cpp                     # ← Modified
│   ├── GuptaPotential.cpp
│   ├── FinnisSinclairPotential.cpp  # ← NEW
│   ├── SuttonChenPotential.cpp      # ← NEW
│   ├── Configuration.cpp             # ← Modified
│   └── ...
│
├── lib/
│   ├── lbfgs.h                      # L-BFGS library
│   └── lbfgs.lib
│
└── results/                          # Output directory
    └── Pt19Co19/
        ├── best.xyz
        ├── best.txt
        ├── historical_best.txt
        └── ...
```

## 📤 Output Files

For each composition (e.g., `results/Pt19Co19/`):

| File | Description |
|------|-------------|
| `best.xyz` | Best structure (XYZ format for visualization) |
| `best.txt` | Best structure (Diamond format) |
| `historical_best.txt` | All-time best across all runs |
| `N_best_energy_per_generation.txt` | Energy convergence for run N |
| `energy.txt` | Population energies |
| `all_compositions_summary.txt` | Summary (if runAllCompositions=true) |

## 🔬 Interatomic Potentials GuCDE

### Which Potential to Use?

| Potential | Best For | Elements | Structure | Accuracy |
|-----------|----------|----------|-----------|----------|
| **Gupta** | Noble & late transition metals | Au, Pt, Pd, Cu, Ag, Ni, Co | FCC, HCP | ★★★★★ |
| **Finnis-Sinclair** | BCC transition metals | Fe, Cr, Mo, W, V, Nb, Ta | BCC | ★★★★★ |
| **Sutton-Chen** | FCC metals (simpler) | Au, Pt, Pd, Cu, Ag, Ni | FCC | ★★★★☆ |

### Common Metal Systems

**Gupta Potential** (Recommended):
- Pt-Co, Pt-Ni (catalysts)
- Au-Cu, Au-Ag (plasmonic)
- Pd-Cu (hydrogen storage)

**Finnis-Sinclair**:
- Fe-Cr (stainless steel)
- Mo-W (high temperature)

**Sutton-Chen**:
- Au-Cu (alternative to Gupta)
- Ag-Cu (alloys)

### Adding New Metal Parameters

#### Step 1: Find Parameters

Search literature (Physical Review B, J. Chem. Phys., etc.)

#### Step 2: Create Parameter File

For Gupta (`data/gupta_Metal1_Metal2.txt`):
```
Metal1 Metal2
r0_11, A_11, xi_11, p_11, q_11     # Metal1-Metal1
r0_22, A_22, xi_22, p_22, q_22     # Metal2-Metal2
r0_12, A_12, xi_12, p_12, q_12     # Metal1-Metal2
```

Use mixing rules for A-B interactions if not available:
- `r0_AB = (r0_AA + r0_BB) / 2`
- `A_AB = sqrt(A_AA * A_BB)`
- `xi_AB = sqrt(xi_AA * xi_BB)`
- `p_AB = (p_AA + p_BB) / 2`
- `q_AB = (q_AA + q_BB) / 2`

#### Step 3: Update Config

```ini
elementA=Metal1
elementB=Metal2
potentialType=Gupta
potentialFile=data/gupta_Metal1_Metal2.txt
```

## 🎨 Visualization

### Recommended Tools

| Tool | Platform | Free | Best For |
|------|----------|------|----------|
| **VESTA** | Win/Mac/Linux | ✅ | Crystals, publication |
| **Ovito** | Win/Mac/Linux | ✅ (Basic) | Analysis, animations |
| **VMD** | Win/Mac/Linux | ✅ | Molecular dynamics |
| **ASE** | Python | ✅ | Quick viewing |
| **Avogadro** | Win/Mac/Linux | ✅ | Molecular editing |

### Quick Visualization

```bash
# With Ovito
ovito results/Pt19Co19/best.xyz

# With ASE (Python)
ase gui results/Pt19Co19/best.xyz

# With VESTA
# Just drag and drop best.xyz onto VESTA window
```

## ⚙️ Configuration Reference

### Potential Selection

```ini
potentialType=Gupta              # Gupta, FinnisSinclair, SuttonChen
potentialFile=data/gupta_PtCo.txt
```

### General

```ini
totalAtoms=38
numElementA=19
numElementB=19
elementA=Pt
elementB=Co
algorithm=CDE                    # CDE, SaNSDE, PSO
```

### CDE Parameters

```ini
CDE.populationSize=40            # Population per sub-population
CDE.maxGenerations=200           # Max generations
CDE.exchangeInterval=20          # Population exchange interval
CDE.useLocalSearch=true          # Enable L-BFGS
CDE.localSearchFrequency=1       # Local search every N evals
CDE.useMultiPopulation=true      # Multiple sub-populations
```

### SaNSDE Parameters

```ini
sansde.populationSize=40         # Population size
sansde.maxGenerations=200        # Max generations
sansde.learningPeriod=20         # Strategy adaptation period
sansde.F_min=0.1                 # Min mutation factor
sansde.F_max=1.0                 # Max mutation factor
sansde.CR_min=0.0                # Min crossover rate
sansde.CR_max=1.0                # Max crossover rate
sansde.useLocalSearch=true       # Enable L-BFGS
```

### Run Control

```ini
numRuns=10                       # Independent runs
runAllCompositions=false         # All compositions?
verbose=true                     # Detailed output
randomSeed=-1                    # -1 = time-based
```

## 💡 Performance Tips

### Potential Selection
- ✅ Start with Gupta for FCC metals
- ✅ Use Finnis-Sinclair for BCC metals
- ✅ Compare potentials if unsure

### Population Size
- Small (<20 atoms): 20-30
- Medium (20-55 atoms): 40-60
- Large (>55 atoms): 60-100

### Generations
- Quick test: 100-200
- Production: 500-1000
- Difficult: 2000+

### Local Search
- ✅ Always enable for best results
- Frequency=1 for small clusters
- Frequency=5-10 for large clusters

### Multiple Runs
- Minimum: 10 runs
- Good: 50 runs
- Publication: 100+ runs

## 🐛 Troubleshooting

### Build Issues

**Cannot find lbfgs.h**
```
Solution: Add lib/ to include directories
Project Properties → C/C++ → General → Additional Include Directories
Add: $(ProjectDir)lib
```

**Linking errors**
```
Solution: Add lbfgs.lib to linker
Project Properties → Linker → Input → Additional Dependencies
Add: $(ProjectDir)lib\lbfgs.lib
```

### Runtime Issues

**Cannot load potential file**
- Check file exists: `dir data\gupta_PtCo.txt`
- Verify path in config.txt
- Ensure format matches potential type

**NaN or Inf energies**
- Check potential parameters are reasonable
- Verify element symbols match file
- Try different potential type

**Program crashes**
- Reduce population size
- Check available memory
- Disable runAllCompositions for testing

## 📚 References

### Potentials

1. **Gupta**: Cleri & Rosato (1993), *Phys. Rev. B*, 48, 22-33
2. **Finnis-Sinclair**: Finnis & Sinclair (1984), *Phil. Mag. A*, 50, 45-55
3. **Sutton-Chen**: Sutton & Chen (1990), *Phil. Mag. Lett.*, 61, 139-146

## 📝 Version History

### v2.0 (Current)
- ✨ Added Finnis-Sinclair potential
- ✨ Added Sutton-Chen potential
- ✨ Runtime potential selection
- ✨ Command-line potential specification
- 🔧 Fixed Chinese comment issues
- 🔧 Improved error handling

### v1.0
- Initial release with Gupta potential
- CDE and SaNSDE algorithms
- Basic configuration system

## 📄 License

MIT License - See LICENSE file for details

---

**Happy Optimizing! 🎉**

For questions: Check config examples, try different potentials, review this README.
