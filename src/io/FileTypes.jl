module FileTypes

# Abstract types allow us to not have to know
# which file type is being used. We simply use
# the type to determine how we process the data.
#
# Each implemented file format needs its own
# struct in order to work properly in the
# ReaderMaster and WriterMaster modules.
export File, Rockstar, AHF
using Parameters

abstract type File end

mutable struct Gallus <: File
    file_name::String
    parameter_file_name::String
end

mutable struct Rockstar <: File
    file_name::String
    parameter_file_name::String
end

mutable struct AHF <: File
    file_name::String
    parameter_file_name::String
end

mutable struct GIZMO <: File
    file_name::String
    parameter_file_name::String
    all_key_transform::Dict
    gas_key_transform::Dict
    stars_key_transform::Dict
    dark_key_transform::Dict
    bhs_key_transform::Dict
    internal_key_to_file_key::Dict

    function GIZMO(file_name::String, parameter_file_name::String)
        gizmo = new()
        gizmo.file_name = file_name
        gizmo.parameter_file_name = parameter_file_name

        gizmo.gas_key_transform = Dict([
            ("id", "ParticleIDs"),
            ("mass", "Masses"),
            ("velocities", "Velocities"),
            ("coordinates", "Coordinates"),
            ("density", "Density"),
            ("temperature", "Temperature"),
            ("specific_internal_energy", "InternalEnergy"),
            ("star_formation_rate", "StarFormationRate"),
            ("size", "SmoothingLength"),
            ("electron_abundance", "ElectronAbundance"),
            ("metal_mass_fraction", "Metallicity"),
            ("hydrogen_mass_fraction", "Metallicity"),
            ("helium_mass_fraction", "Metallicity")
        ])

        gizmo.stars_key_transform = Dict([
            ("id", "ParticleIDs"),
            ("mass", "Masses"),
            ("velocities", "Velocities"),
            ("coordinates", "Coordinates"),
            ("stellar_age", "StellarFormationTime"),
            ("metal_mass_fraction", "Metallicity"),
            ("hydrogen_mass_fraction", "Metallicity"),
            ("helium_mass_fraction", "Metallicity")
        ])

        gizmo.dark_key_transform = Dict([
            ("id", "ParticleIDs"),
            ("mass", "Masses"),
            ("velocities", "Velocities"),
            ("coordinates", "Coordinates")
        ])

        gizmo.bhs_key_transform = Dict([
            ("id", "ParticleIDs"),
            ("mass", "Masses"),
            ("velocities", "Velocities"),
            ("coordinates", "Coordinates"),
            ("subgrid_mass", "BH_Mass")
        ])

        gizmo.internal_key_to_file_key = Dict([
            ("gas", "PartType0"),
            ("stars", "PartType4"),
            ("dark", "PartType1"),
            ("bhs", "PartType5")
        ])

        gizmo
    end
end

end