module Gallus

# Master include list 
using Unitful, UnitfulAstro
function convert_units(data_array::Vector, new_unit::Unitful.FreeUnits)
    ustrip.(new_unit, data_array) * new_unit
end

include("io/ReaderMaster.jl")
include("Simulation.jl")
include("Halos.jl")


end