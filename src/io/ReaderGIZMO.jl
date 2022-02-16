# ReaderGIZMO.jl
#
# Requires the following files included beforehand:
# FileTypes
# Cells
using HDF5, Printf, StructArrays, ..FileTypes

# These don't really need to be super generic 
# and available across the simulation, so we
# modularize them here.
module ReaderGIZMO
global header
# These are the only keys we will keep in the internal particle_dict
# so that we don't use too much memory. We won't use the other keys, 
# unless the code is modified in some way.
global readable_keys = ["Masses",
                        "Coordinates",
                        "Velocities",
                        "Density",
                        "SmoothingLength",
                        "Metallicity",
                        "ElectronAbundance",
                        "InternalEnergy",
                        "StarFormationRate",
                        "StellarFormationTime",
                        "BH_Mass",
                        "ParticleIDs"]

set_header(data) = (global header = data)
import ...FileTypes
import PhysicalConstants.CODATA2018: k_B, m_p
using HDF5, StructArrays, Unitful, UnitfulAstro, UnitfulEquivalences, Printf
using ProgressMeter

"""
    particle_key_to_idx(key::String)

Converts PartTypeN to N+1

...
# Arguments
- `key::String`: PartType0, PartType1, etc.
...
"""
function particle_key_to_idx(key::String) parse(UInt8, last(split(key, ""))) + 1 end

"""
    read_header(file_container::FileTypes.GIZMO)

Read the HDF5 Header key from the simulation snapshot into a dictionary.

...
# Arguments
- `file_container::FileTypes.GIZMO`: The prepared GIZMO container (file name/param file)
...
"""
function read_header(file_container::FileTypes.GIZMO)
    file_handle = h5open(file_container.file_name, "r")

    # Open the header and determine if we are reading
    # more than one file.
    header_dict = Dict()
    for key in keys(attributes(file_handle["Header"]))
        header_dict[key] = read_attribute(file_handle["Header"], key)
    end
    close(file_handle)
    set_header(header_dict)
    header_dict
end

"""
    preallocate_particle_dict(particle_dict::Dict, file_container::FileTypes.GIZMO)

When a simulation has multiple snapshot files, it is better to preallocate the entire
internal dictionary that stores the particle data. Loop through all of the particle
types and grab the sizes of the Arrays, and then allocate them internally. 

...
# Arguments
- `particle_dict::Dict`: The dictionary containing all of the snapshot data.
- `file_container::FileTypes.GIZMO`: The prepared GIZMO container (file name/param file)
...
"""
function preallocate_particle_dict(particle_dict::Dict, file_container::FileTypes.GIZMO)
    print("ReaderGIZMO::preallocate_particle_dict\n")
    global readable_keys
    # Loop through all of the available keys in the first file to determine the size, dim, and 
    # type of the datasets.
    message = "Preallocating internal particle dictionary..."
    @showprogress 0.01 message for (i, number_of_particles) in enumerate(header["NumPart_Total"])
        if number_of_particles == 0 continue end
        pt = i - 1
        file_handle = h5open(file_container.file_name)
        key = "PartType$pt"
        particle_dict[key] = Dict()
        for property_key in keys(file_handle[key])
            if !(property_key in readable_keys) continue end
            full_key = "$key/$property_key"
            data_type = eltype(file_handle[full_key])
            # Something like Metallicity has shape (11, NumPart_ThisFile)
            # but we want to save (11, NumPart_Total) for everything.
            if ndims(file_handle[full_key]) > 1
                particle_dict[key][property_key] = Array{data_type}(undef,
                    (
                        size(file_handle[full_key])[1], 
                        header["NumPart_Total"][i]
                    )
                )
            else
                particle_dict[key][property_key] = Array{data_type}(undef,
                    header["NumPart_Total"][i]
                )
            end
        end
        close(file_handle)
    end
end

"""
    build_particle_dict_multi(particle_dict::Dict, file_container::FileTypes.GIZMO)

Copies the data from multiple HDF5 files into a dictionary reflecting the exact structure.

...
# Arguments
- `particle_dict::Dict`: The dictionary containing all of the snapshot data.
- `file_container::FileTypes.GIZMO`: The prepared GIZMO container (file name/param file)
...
"""
function build_particle_dict_multi(particle_dict::Dict, file_container::FileTypes.GIZMO)
    print("ReaderGIZMO::build_particle_dict_multi\n")
    global readable_keys
    file_base, dummy, file_extension = split(file_container.file_name, ".")

    # This will contain the starting offset for the preallocated arrays.
    # Update every time we open a new file by the NumPart_ThisFile
    # value in the header to get the correct offset.
    idx_slice = [1 for i in 1:length(header["NumPart_Total"])]

    # First, we will allocate all of the arrays.
    preallocate_particle_dict(particle_dict, file_container)

    for file_number in 0:header["NumFilesPerSnapshot"]-1
        file_handle = h5open("$file_base.$file_number.$file_extension", "r")
        particle_count_this_file = read_attribute(file_handle["Header"], 
                                             "NumPart_ThisFile")
        for (i, number_of_particles) in enumerate(particle_count_this_file)
            # Particle type index, since Julia starts at 1
            if number_of_particles == 0 continue end
            pt = i - 1
            key = "PartType$pt"
            lower_bound = idx_slice[i]
            upper_bound = idx_slice[i] + number_of_particles - 1
            for property_key in keys(file_handle[key])
                if !(property_key in readable_keys) continue end
                if ndims(particle_dict[key][property_key]) > 1
                    particle_dict[key][property_key][:, lower_bound:upper_bound] = read(
                        file_handle["$key/$property_key"]
                    )
                else
                    particle_dict[key][property_key][lower_bound:upper_bound] = read(
                        file_handle["$key/$property_key"]
                    )
                end
            end

            idx_slice[i] += number_of_particles
        end

        close(file_handle)
    end
end

"""
    build_particle_dict_single(particle_dict::Dict, file_container::FileTypes.GIZMO)

Copies the data from the HDF5 file into a dictionary reflecting the exact structure.

...
# Arguments
- `particle_dict::Dict`: The dictionary containing all of the snapshot data.
- `file_container::FileTypes.GIZMO`: The prepared GIZMO container (file name/param file)
...
"""
function build_particle_dict_single(particle_dict::Dict, file_container::FileTypes.GIZMO)
    print("ReaderGIZMO::build_particle_dict_single\n")
    global readable_keys
    file_handle = h5open(file_container.file_name, "r")
    for key in keys(file_handle)
        if key == "Header" continue end
        if header["NumPart_Total"][particle_key_to_idx(key)] == 0 continue end
        particle_dict[key] = Dict()
        for property_key in keys(file_handle[key])
            if !(property_key in readable_keys) continue end
            particle_dict[key][property_key] = read(file_handle["$key/$property_key"])
        end
    end
    close(file_handle)
    particle_dict
end

"""
    get_units_from_file(parameter_file_name::String)

Reads the GIZMO parameter file to look for the 4 unit
definitions. If they are not defined, they are assumed
to be unity in the g, cm, (km/s), and Gauss system.

...
# Arguments
- `parameter_file_name::String`: The parameter file name (can include directory)
...
"""
function get_units_from_file(parameter_file_name::String)
    print("ReaderGIZMO::get_units_from_file\n")
    file_handle = open(parameter_file_name)
    for i in eachline(file_handle)
        if occursin(r"UnitLength_in_cm", i) global unit_length = split(i)[2] end
        if occursin(r"UnitMass_in_g", i) global unit_mass = split(i)[2] end
        if occursin(r"UnitVelocity_in_cm_per_s", i) global unit_velocity = split(i)[2] end
        if occursin(r"UnitMagneticField_in_gauss", i) global unit_magnetic = split(i)[2] end
    end
    close(file_handle)

    # Let's just get rid of little h now so we never have to think about
    # it again.
    Dict([("unit_length", parse(Float64, unit_length) * 1.0u"cm" / header["HubbleParam"]), 
         ("unit_mass", parse(Float64, unit_mass) * 1.0u"g" / header["HubbleParam"]), 
         ("unit_velocity", parse(Float64, unit_velocity) * 1.0u"cm/s"), 
         ("unit_magnetic", parse(Float64, unit_magnetic) * 1.0u"Gauss")])
end

"""
    convert_raw_data_to_units(particle_dict::Dict, file_container::FileTypes.GIZMO)

Convert the raw data from the simulation snapshot into its
own units, given by the parameter file that was provided.
We will convert to internal units in build_cell_list(). 

Everything will be in physical units with no little h.

readable_keys inside of the ReaderGIZMO module must reflect the keys
present in this function.  These are the only keys that translate into
the internal structure.

...
# Arguments
- `particle_dict::Dict`: The dictionary containing all of the snapshot data.
- `file_container::FileTypes.GIZMO`: The prepared GIZMO container (file name/param file)
...
"""
function convert_raw_data_to_units(particle_dict::Dict, file_container::FileTypes.GIZMO)
    print("ReaderGIZMO::convert_raw_data_to_units\n")
    # We have to parse the parameter file to find the units used in the simulation
    units_dict = get_units_from_file(file_container.parameter_file_name)

    scale_factor = 1.0
    try
        if header["OmegaLambda"] != 0
            scale_factor = header["Time"]
        end
    catch KeyError
        # Newer versions of GIZMO have this!
        header["OmegaLambda"] = copy(header["OmegaLambda"])
        if header["Omega_Lambda"] != 0
            scale_factor = header["Time"]
        end
    end

    mass_conv = units_dict["unit_mass"]
    length_conv = units_dict["unit_length"] * scale_factor
    velocity_conv = units_dict["unit_velocity"] / sqrt(scale_factor)
    time_conv = length_conv / (velocity_conv * sqrt(scale_factor))
    metal_conv = 1.0 / 0.0134  # Asplund 2009
    specific_energy_conv = (velocity_conv * sqrt(scale_factor))^2

    # Assume gamma = 5.0 / 3.0
    m_p_gamma_over_k_B = m_p * (5.0 / 3.0 - 1.0) / k_B
    y_helium(x) = x / (1.0 - x)
    mu(x, y) = (1.0 + 4.0 * x) / (1.0 + x + y)
    # x: helium mass fraction, y: electron abundance, z: specific internal energy
    temperature(x, y, z) = m_p_gamma_over_k_B * mu(y_helium(x), y) * z
    output_dict = Dict() 
    function convert_and_delete(output_dict::Dict, input_dict::Dict, 
                                key::String, conv::Number)
        output_dict[key] = broadcast(.*, input_dict[key], conv)
        delete!(input_dict, key)
    end

    message = "Converting data to Unitful quantities..."
    @showprogress 0.01 message for key in keys(particle_dict)
        if key == "Header" continue end
        output_dict[key] = Dict()
        convert_and_delete(output_dict[key], particle_dict[key], "ParticleIDs", 1)
        convert_and_delete(output_dict[key], particle_dict[key], "Masses", mass_conv)
        convert_and_delete(output_dict[key], particle_dict[key], "Coordinates", length_conv)
        convert_and_delete(output_dict[key], particle_dict[key], "Velocities", velocity_conv)

        if key == "PartType0"
            convert_and_delete(output_dict[key], particle_dict[key], "InternalEnergy", specific_energy_conv)
            output_dict[key]["Temperature"] = temperature.(
                particle_dict[key]["Metallicity"][2],
                particle_dict[key]["ElectronAbundance"],
                output_dict[key]["InternalEnergy"]
            )
            convert_and_delete(output_dict[key], particle_dict[key], "ElectronAbundance", 1.0)
            # This provides a way out if they don't provide the proper units somehow
            output_dict[key]["Temperature"] = broadcast(
                .*, 
                ustrip.(u"K", output_dict[key]["Temperature"]), 
                u"K"
            )

            convert_and_delete(output_dict[key], particle_dict[key], "Density", mass_conv / length_conv^3)
            convert_and_delete(output_dict[key], particle_dict[key], "SmoothingLength", length_conv)
            convert_and_delete(output_dict[key], particle_dict[key], "StarFormationRate", mass_conv / time_conv)
        end

        if key == "PartType0" || key == "PartType4"
            convert_and_delete(output_dict[key], particle_dict[key], "Metallicity", metal_conv)
        end

        if key == "PartType4"
            convert_and_delete(output_dict[key], particle_dict[key], "StellarFormationTime", 1.0)
        end

        if key == "PartType5"
            convert_and_delete(output_dict[key], particle_dict[key], "BH_Mass", mass_conv)
        end
    end

    output_dict
end

"""
    restructure_particle_dict(particle_dict::Dict)

Convert the particle dictionary from GIZMO format to the internal 
dictionary format. The keys are gas, stars, dark, and bhs.

...
# Arguments
- `particle_dict::Dict`: The dictionary containing all of the snapshot data.
...
"""
function restructure_particle_dict(particle_dict::Dict)
    output_dict = Dict([
        ("gas", Dict()),
        ("stars", Dict()),
        ("dark", Dict()),
        ("bhs", Dict())
    ])

    all_key_transform = Dict([
        ("ParticleIDs", "ids"),
        ("Masses", "masses"),
        ("Velocities", "velocities"),
        ("Coordinates", "coordinates")
    ])

    gas_key_transform = Dict([
        ("Density", "densities"),
        ("Temperature", "temperatures"),
        ("StarFormationRate", "starformationrates"),
        ("SmoothingLength", "smoothinglengths"),
        ("Metallicity", "metallicities")
    ])

    stars_key_transform = Dict([
        ("Metallicity", "metallicities"),
        ("StellarFormationTime", "stellarages")
    ])

    bhs_key_transform = Dict([
        ("BH_Mass", "subgridmasses")
    ])

    dark_key_transform = Dict([])

    function internal_key_to_part_type(internal::String)
        if internal == "gas" return "PartType0"
        elseif internal == "dark" return "PartType1"
        elseif internal == "stars" return "PartType4"
        elseif internal == "bhs" return "PartType5"
        else return "" end
    end


    message = "Copying data to Simulation structure..."
    @showprogress 0.01 message for key in keys(output_dict)
        # First copy all of the common data from the common keys list
        part_type_key = internal_key_to_part_type(key)
        for property_key in keys(all_key_transform)
            output_dict[key][all_key_transform[property_key]] = particle_dict[part_type_key][property_key]
        end

        # Now copy the unique data to each data set
        if key == "gas"
            key_transform = gas_key_transform
        elseif key == "stars"
            key_transform = stars_key_transform
        elseif key == "bhs"
            key_transform = bhs_key_transform
        else
            key_transform = dark_key_transform
        end

        for property_key in keys(key_transform)
            if property_key == "Metallicity"
                # TODO: Accept all metallicities!
                output_dict[key][key_transform[property_key]] = particle_dict[part_type_key][property_key][1, :]
            else
                output_dict[key][key_transform[property_key]] = particle_dict[part_type_key][property_key]
            end
        end
    end

    output_dict
end

end

"""
Implementation of the read_file(FilesTypes.XYZ) method for the GIZMO
hydrodynamical simulation code. 

...
# Arguments
- `file_container::FileTypes.GIZMO`: The prepared GIZMO container (file name/param file)
...
"""
function read_file(file_container::FileTypes.GIZMO)
    # Read the header. We need this to know how many files
    # to open.
    header = ReaderGIZMO.read_header(file_container)

    # particle_dict contains a full representation of the hdf5
    # file but inside of our module. We will keep a representation 
    # of the simulation data according to our definitions.
    particle_dict = Dict()

    if header["NumFilesPerSnapshot"] > 1
        ReaderGIZMO.build_particle_dict_multi(particle_dict, file_container)
    else
        ReaderGIZMO.build_particle_dict_single(particle_dict, file_container)
    end

    particle_dict = ReaderGIZMO.convert_raw_data_to_units(copy(particle_dict), file_container)
    ReaderGIZMO.restructure_particle_dict(particle_dict)
end