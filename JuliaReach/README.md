# JuliaReach experiments

This manual describes the steps to run the JuliaReach experiments.



## Installation

- Install Julia (download [here](https://julialang.org/downloads/)).

- Start the REPL.

- Add the package implementing our method plus some additional dependencies for simulation and plotting.

    - The recommended way is to start Julia from the `JuliaReach` folder with the parameter

    ```bash
    julia --project=.
    ```
    Here `.` is the path to the files `Manifest.toml` and `Project.toml`. This way Julia will use the package versions defined in `Manifest.toml`, which were used at the time when this package was prepared.
    - To install the packages in a fresh version, enter the package mode (via typing `]`) and then enter the following commands:

    ```julia
    pkg> dev https://github.com/sisl/NeuralVerification.jl

    pkg> add https://github.com/JuliaReach/NeuralNetworkAnalysis.jl

    pkg> add DifferentialEquations Plots LaTeXStrings
    ```


## Running the experiments

Go to the `JuliaReach` folder, open a REPL, and load the experiments like this:

```julia
julia> include("Sherlock-Benchmark-10-Unicycle/Sherlock-Benchmark-10-Unicycle.jl")

julia> include("Sherlock-Benchmark-9-TORA/Sherlock-Benchmark-9-TORA.jl")

julia> include("ACC/ACC.jl")

julia> include("Single-Pendulum/Single-Pendulum.jl")

julia> include("Double-Pendulum/Double-Pendulum.jl")

julia> include("Airplane/Airplane.jl")
```

Note that Julia is just-in-time compiled and hence the first run is slow. That is why we run a short warm-up version in the benchmark scripts.



## Plots

The plots are automatically stored in the respective subfolders.
