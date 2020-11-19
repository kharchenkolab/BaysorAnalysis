# BaysorAnalysis

Analysis for the Baysor paper. See the [Baysor](https://github.com/hms-dbmi/Baysor) package and the [pre-print](https://doi.org/10.1101/2020.10.05.326777).

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/).

To (locally) reproduce this project, do the following:

1. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
2. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts.


## Visualization of the results

- [Allen sm-FISH](http://vitessce.io/?url=https%3A%2F%2Fsealver.in%2Fvitessce%2Fallen_sm_fish.json&theme=dark)
