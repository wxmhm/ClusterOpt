# Binary Alloy Cluster Global Optimization

A C++ implementation of evolutionary algorithms for optimizing binary alloy nanoparticle structures using the Gupta potential.

## Features

- **Two Optimization Algorithms:**
  - **IDE**: Improved Differential Evolution with multiple mutation strategies
  - **SaNSDE**: Self-adaptive Differential Evolution with Neighborhood Search

- **Local Optimization**: L-BFGS method for structure refinement
- **Flexible Configuration**: Support for various element combinations and cluster sizes
- **Multi-run Capability**: Statistical analysis across multiple independent runs
- **Historical Best Tracking**: Automatically saves and compares with previous results

## Requirements

- Windows 10 or later
- Visual Studio 2019 or later
- liblbfgs library (included in project)

## Quick Start

### 1. Build the Project

1. Open `BinaryAlloyOptimizer.sln` in Visual Studio 2019
2. Select **Release** configuration and **x64** platform
3. Build → Build Solution (or press `Ctrl+Shift+B`)
4. Executable will be in `x64/Release/BinaryAlloyOptimizer.exe`

### 2. Prepare Configuration File

Edit `data/config.txt`:
totalAtoms=38
numElementA=19
numElementB=19
elementA=Pt
elementB=Co
algorithm=IDE
ide.populationSize=40
ide.maxGenerations=200
numRuns=10

### 3. Prepare Gupta Potential Parameters
Create/edit data/gupta_PtCo.txt:
Pt Co
2.6975, 0.2975, 1.64249, 10.612, 4.004
2.4985, 0.1765, 1.18997, 10.960, 3.814
2.5980, 0.2370, 1.41623, 10.786, 3.909
Format: Three lines with r0, A, xi, p, q for Pt-Pt, Co-Co, and Pt-Co interactions

### 4. Run the Program
Method 1: Double-click (Recommended)
1. Simply double-click BinaryAlloyOptimizer.exe
2. It will automatically load data/config.txt
3. Results will be saved in results/ folder

Method 2: Command Prompt
1. cd path\to\your\project
2. x64\Release\BinaryAlloyOptimizer.exe

Method 3: With command line parameters
1. x64\Release\BinaryAlloyOptimizer.exe -atoms 38 -composition 19 -elements Pt Co -runs 10

## Quick Start Summary:
1. Build in Visual Studio 2019 (Release x64)
2. Prepare data/config.txt and data/gupta_PtCo.txt
3. Double-click BinaryAlloyOptimizer.exe
4. Find results in results/ folder
5. Visualize with Ovito or other tools

## Study Cases
### Case 1: Pt19Co19 (38 atoms) - Core-Shell Structure
Configuration (data/config.txt):
totalAtoms=38
numElementA=19
numElementB=19
elementA=Pt
elementB=Co
algorithm=IDE
ide.populationSize=40
ide.maxGenerations=200
ide.useLocalSearch=true
numRuns=50
runAllCompositions=false

### Case 2: Pt7Co6 (13 atoms) - Small Cluster
Configuration (data/config.txt):
totalAtoms=38
numElementA=19
numElementB=19
elementA=Pt
elementB=Co
algorithm=IDE
ide.populationSize=40
ide.maxGenerations=200
ide.useLocalSearch=true
numRuns=10
runAllCompositions=false

### Case 3: Pure Pt38 Cluster
Configuration (data/config.txt):
totalAtoms=38
numElementA=38
numElementB=0
elementA=Pt
elementB=Pt
algorithm=IDE
numRuns=5
runAllCompositions=false

### Case 4: Composition Screening Pt-Co (38 atoms)
totalAtoms=38
numElementA=38
numElementB=0
numRuns=5
algorithm=IDE
ide.maxGenerations=200
elementA=Pt
elementB=Co
runAllCompositions=true

### Case 5: Au-Cu System (55 atoms)
Step 1: Create data/gupta_AuCu.txt:
Au Cu
2.8840, 0.2061, 1.7900, 10.229, 4.036
2.5562, 0.0894, 1.2240, 10.960, 4.000
2.7201, 0.1355, 1.5070, 10.595, 4.018

Step 2: Configure(data/config.txt):
totalAtoms=55
numElementA=28
numElementB=27
elementA=Au
elementB=Cu
potentialFile=data/gupta_AuCu.txt
algorithm=IDE
ide.maxGenerations=300
numRuns=10

## Project Structure
BinaryAlloyOptimizer/
├── data/
│   ├── config.txt              # Main configuration file
│   ├── gupta_PtCo.txt         # Gupta potential parameters
│   └── initial_structures/     # (Optional) Initial structures
├── include/                    # Header files
│   ├── Common.h
│   ├── BinaryAlloyCluster.h
│   ├── GuptaPotential.h
│   ├── IDE.h
│   ├── SaNSDE.h
│   └── ...
├── src/                        # Source files
│   ├── main.cpp
│   ├── BinaryAlloyCluster.cpp
│   ├── IDE.cpp
│   └── ...
├── x64/
│   └── Release/
│       └── BinaryAlloyOptimizer.exe
└── results/                    # Output directory (auto-created)
    └── Pt19Co19/
        ├── best.xyz
        ├── best.txt
        ├── historical_best.txt
        └── ...

## Output Files
For each composition (e.g., results/Pt19Co19/):
1. best.xyz：Best structure in XYZ format
2. best.txt：Best structure in Diamond format
3. historical_best.txt：All-time best structure
4. N_best_energy_per_generation.txt：Energy convergence for run N
5. energy.txt：Population energies during optimization

## Algorithm Details
### IDE (Improved Differential Evolution)
1. Uses multiple mutation strategies (RAND/1, BEST/1, RAND/2, etc.)
2. Multiple sub-populations with different strategies
3. Population exchange mechanism
4. Sphere-cut-splice crossover operator
5. Atom swap operator for binary alloys

### SaNSDE (Self-adaptive DE with Neighborhood Search)
1. Self-adaptive F and CR parameters
2. Strategy success-based probability adjustment
3. Cauchy and Gaussian distributions for parameter generation
4. Neighborhood search for local refinement
5. Success memory for parameter learning

## Adding New Metal Systems
1. Find Gupta Parameters
2. Create Parameter File
Format data/gupta_MetalA_MetalB.txt:
MetalA MetalB
r0_AA, A_AA, xi_AA, p_AA, q_AA
r0_BB, A_BB, xi_BB, p_BB, q_BB
r0_AB, A_AB, xi_AB, p_AB, q_AB

3. Update Configuration
elementA=YourMetalA
elementB=YourMetalB
potentialFile=data/gupta_MetalA_MetalB.txt

## Visualization
Recommended tools for visualizing .xyz files:
VESTA - Free, Windows/Mac/Linux
Ovito - Free basic version available
VMD - Free, powerful
ASE - Python library

## Configuration Parameters Reference
### General Parameters
totalAtoms=38              # Total number of atoms
numElementA=19             # Number of element A atoms
numElementB=19             # Number of element B atoms
elementA=Pt                # Element A symbol
elementB=Co                # Element B symbol
algorithm=IDE              # IDE or SaNSDE

### IDE Parameters
ide.populationSize=40              # Population size per sub-population
ide.maxGenerations=200             # Maximum generations
ide.exchangeInterval=20            # Interval for population exchange
ide.useLocalSearch=true            # Enable L-BFGS local search
ide.localSearchFrequency=1         # Apply local search every N evaluations
ide.useMultiPopulation=true        # Use multiple sub-populations

### SaNSDE Parameters
sansde.populationSize=40           # Population size
sansde.maxGenerations=200          # Maximum generations
sansde.learningPeriod=20           # Period for strategy adaptation
sansde.F_min=0.1                   # Minimum mutation factor
sansde.F_max=1.0                   # Maximum mutation factor
sansde.CR_min=0.0                  # Minimum crossover rate
sansde.CR_max=1.0                  # Maximum crossover rate
sansde.p_min=0.1                   # Minimum strategy probability
sansde.neighborhoodSizeMin=2       # Min neighborhood search size
sansde.neighborhoodSizeMax=15      # Max neighborhood search size
sansde.memorySize=100              # Success memory size

### Run Control
numRuns=10                 # Number of independent runs
runAllCompositions=false   # Optimize all compositions
verbose=true               # Show detailed output
randomSeed=-1              # -1 for time-based, or fixed integer

## License
MIT License - See LICENSE file for details
