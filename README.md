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

2. Open the Julia REPL and type the following commands
```Julia
push!(LOAD_PATH, "path/to/GM-Julia/src")
using GM_main
species_list, reaction_list, system, speciesIDs, output_list = run_GM(input)
```

# Input deck
## Constants
## System
## Species
## Chemical reactions
## Outputs
