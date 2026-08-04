#!/bin/bash
julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate(); Pkg.add("PalmerPenguins")'