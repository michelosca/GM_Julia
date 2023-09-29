# GM Julia
Portable chemical-kinetics global model for low temperature plasmas written in Julia

# Content
- [Run GM Julia](#run-gm-julia)
- [Input deck](#input-deck)
  * [Constants](#constants)
  * [System](#system)
  * [Species](#species)
  * [Chemical reactions](#chemical-reactions)
  * [Outputs](#outputs)

# Run GM Julia
1. Download GM Julia from
```
https://github.com/michelosca/GM_Julia/archive/refs/heads/master.zip
```
or clone the git repository
```bash
git clone https://github.com/michelosca/GM_Julia.git
```

2. Install the following libraries
```
pkg> add CSV, DataFrames, Printf, Dates, TimerOutputs, Roots, DifferentialEquations
```

3. Open the Julia REPL and type the following commands
```Julia
push!(LOAD_PATH, "path/to/GM-Julia/src")
using GM_main: run_GM
species_list, reaction_list, system, speciesIDs, output_list = run_GM(input)
```
where ``input`` is a variable of type ``String`` with the path to the input deck file

# Input deck
## Constants
## System
## Species
## Chemical reactions
## Outputs
