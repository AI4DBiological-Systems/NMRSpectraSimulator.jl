module NMRSpectraSimulator

using Distributed, LinearAlgebra, FFTW, SharedArrays, OffsetArrays

import DelimitedFiles, Unicode, FileIO, Random, Printf, PyCall, FiniteDiff, Optim, NLopt, Interpolations, BSON, JLD, Kronecker, BenchmarkTools, NearestNeighbors, HTTP, JSON3, Dates

#Pkg.add("https://github.com/AI4DBiological-Systems/GISSMOReader")
import GISSMOReader # https://github.com/AI4DBiological-Systems/GISSMOReader

include("../src/types.jl")
include("../src/utils.jl")

include("../src/SH/SH_front_end.jl")
#include("../src/SH/molecule.jl")
include("../src/SH/SH.jl")
include("../src/SH/Hamiltonian.jl")
include("../src/SH/operators.jl")
include("../src/SH/resonance_partitioning.jl")

#include("../src/IO/load_configs.jl")
#include("../src/library/GISSMO_entries.jl")
include("../src/conversions/cs_parse.jl")

include("../src/model_evals/lorentzian.jl")


#include("../src/library/database_helpers.jl")
#include("../src/model/parse.jl")
#include("../src/misc/utilities.jl")
#include("../src/misc/front_end.jl")
end
