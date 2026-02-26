# nbody

An interactive 2D gravitational N-body simulator written in Julia.

## Features

- **RK4 integration** — 4th-order Runge-Kutta for accurate orbit calculation
- **Barnes-Hut tree** — optional O(n log n) force approximation via a quadtree, can be toggled at runtime
- **Softening** — prevents forces from diverging at close range
- **Interactive GUI** — built with GLMakie; control speed, load/save simulations and generate random particle clouds at runtime

## Running

```bash
julia --project=. src/main.jl
```

First run will take a few minutes to precompile GLMakie. Subsequent runs start quickly.

## GUI Controls

| Control | Description |
|---|---|
| Vertical slider (left) | Simulation speed — number of concurrent timer threads (0 = paused, max 20) |
| **Use Barnes-Hut?** toggle | Switch between O(n log n) BH and O(n²) direct-sum force calculation |
| **MAC** textbox | Multipole Acceptance Criterion for Barnes-Hut (default 1; lower = more accurate, slower) |
| **Particle max speed** textbox | Clamp particle velocities to this value |
| **Save simulation to:** | Type a filename (no extension) and click the button to save to `data/<name>.txt` |
| **Load simulation from:** | Type a filename (no extension) and click the button to load from `data/<name>.txt` |
| **Generate particles** | Add a random cloud of particles using the n / mass mean / mass stddev / velocity stddev fields |

## Project Structure

```
src/
  main.jl                 — GUI, simulation loop, controls
  iteration-algorithms.jl — Barnes-Hut tree, RK4 integration (BH and direct-sum)
  particle.jl             — Particle struct, direct-sum acceleration, RK4 integrator
  saving.jl               — Binary save/load
  utils.jl                — Physical constants (G, Δt, softening, etc.)
  utils/maybe.jl          — Maybe/Just/Nothing type
data/                     — Saved simulations (.txt)
```

## Dependencies

- [GLMakie.jl](https://github.com/MakieOrg/Makie.jl) — rendering and GUI
- [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) — fixed-size vectors for positions/velocities
- [DataStructures.jl](https://github.com/JuliaCollections/DataStructures.jl) — queue for BH traversal
- [Combinatorics.jl](https://github.com/JuliaMath/Combinatorics.jl) — pairwise combinations for direct-sum
- [PartialFunctions.jl](https://github.com/archermarx/PartialFunctions.jl) — `$` partial application operator
