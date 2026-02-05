# HHGsimPL -- Simulation of High-order harmonic generation using Picard-Lefschetz integration methods

[![Build Status](https://github.com/anneaux/HHGsimPL.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/anneaux/HHGsimPL.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/anneaux/HHGsimPL.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/anneaux/HHGsimPL.jl)

Find a small installation help [here](docs/Installation.md), as well as a jupyter notebook that showcases the main functionalities for getting HHG dipoles in [the same folder](docs).

This is the (very rudimentary) code that I've used to solve the SFA-based HHG integral with the methods based on Picard-Lefschetz theory. 
At the moment, this repository contains all the necessary methods. 
In the very near future, the actual integration techniques will be separated out into their own julia package, to be found here: [https://github.com/anneaux/PicardLefschetz.jl].

This code here has been used for the publication [https://arxiv.org/abs/2510.12545], and there are explanations of the underlying theory, the implementation and application given there.
At some point there will be a proper documentation here. For now, please have a look at this paper, or just message me!


## Common pitfalls
- when using check_contribution(), make sure that the action is defined with the correct sign. I.e., the integral is expected to be of form exp(-im f) *I THINK* and then the first check is whether real(f(ts)) < 0. 
- check the flowstepfactor, subdividethreshold etc.
