module NMRSpectraSimulator

using Distributed, LinearAlgebra, FFTW, SharedArrays, OffsetArrays, Statistics

import Kronecker, Interpolations, NearestNeighbors, JSON, Graphs
# import FiniteDiff, Optim, NLopt, Random, Unicode, FileIO, DelimitedFiles
# import Printf, PyCall, BSON, JLD, , BenchmarkTools, , HTTP, JSON3, Dates

# Julia quirk: must manually install unregistered packages (like GISSMOReader) on local machine before installing this NMRSpectraSimulator package.
#Pkg.add("https://github.com/AI4DBiological-Systems/GISSMOReader")
#import GISSMOReader # https://github.com/AI4DBiological-Systems/GISSMOReader.jl

include("../src/types.jl")
include("../src/types2.jl")

include("../src/utils.jl")

#include("../src/SH/SH_front_end.jl")
#include("../src/SH/molecule.jl")
include("../src/SH/SH.jl")
include("../src/SH/Hamiltonian.jl")
include("../src/SH/operators.jl")
include("../src/SH/resonance_partitioning.jl")

#include("../src/IO/load_configs.jl")
#include("../src/library/GISSMO_entries.jl")
include("../src/conversions/cs_parse.jl")

include("../src/model_evals/lorentzian.jl")
#include("../src/model_evals/lorentzian_z.jl")

include("../src/model_evals/setup_itp.jl")
include("../src/model_evals/itp_evals.jl")

include("../src/conversions/IO_csJ.jl")
include("../src/conversions/IO_mag_eq.jl")


#include("../src/library/database_helpers.jl")
#include("../src/model/parse.jl")
#include("../src/misc/utilities.jl")
#include("../src/misc/front_end.jl")

export combinevectors,
    CompoundFIDType,

    fitproxies!,
    setupmixtureproxies,
    evalmixture,
    evalitpproxymixture

end
