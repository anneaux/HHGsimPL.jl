# HHGsimPL

[![Build Status](https://github.com/anneaux/HHGsimPL.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/anneaux/HHGsimPL.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/anneaux/HHGsimPL.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/anneaux/HHGsimPL.jl)

Find a small installation help [here](docs/Installation.md), as well as a jupyter notebook that showcases the main functionalities for getting HHG dipoles in [the same folder](docs)

## Common pitfalls
- when using check_contribution(), make sure that the action is defined with the correct sign. I.e., the integral is expected to be of form exp(-im f) *I THINK* and then the first check is whether real(f(ts)) < 0. 
- check the flowstepfactor, subdividethreshold etc.
