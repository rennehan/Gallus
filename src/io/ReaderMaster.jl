module Reader

available_fields_by_type = Dict{String, Tuple}([
    "mass" => ("mass", "subgrid_mass"),
    "length" => ("coordinates", "size"),
    "time" => ("stellar_age",),
    "velocity" => ("velocities",),
    "density" => ("density",),
    "metal_mass_fraction" => ("metal_mass_fraction", "hydrogen_mass_fraction", "helium_mass_fraction"),
    "mass_rate" => ("star_formation_rate", "accretion_rate"),
    "specific_energy" => ("specific_internal_energy",),
    "temperature" => ("temperature",),
    "unitless" => ("id", "electron_abundance")
])

data_types_by_field_name = Dict{String, Type}([
    "id" => UInt64,
    "mass" => Float64,
    "subgrid_mass" => Float64,
    "stellar_age" => Float64,
    "velocities" => Float64,
    "coordinates" => Float64,
    "density" => Float64,
    "metal_mass_fraction" => Float64,
    "hydrogen_mass_fraction" => Float64,
    "helium_mass_fraction" => Float64,
    "star_formation_rate" => Float64,
    "accretion_rate" => Float64,
    "specific_internal_energy" => Float64,
    "temperature" => Float64,
    "electron_abundance" => Float64
])

include("ReaderGallus.jl")
include("ReaderGIZMO.jl")
include("ReaderRockstar.jl")

end