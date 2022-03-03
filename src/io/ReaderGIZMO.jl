using ..FileTypes, ..Tools
using Unitful, UnitfulAstro
import PhysicalConstants.CODATA2018: k_B, m_p

# These don't really need to be super generic 
# and available across the simulation, so we
# modularize them here.
module ReaderGIZMO
global header

set_header(data) = (global header = data)
set_units(data) = (global units = data)
import ...FileTypes
using HDF5, Unitful, UnitfulAstro
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
    unit_dict = Dict([
        ("unit_length", parse(Float64, unit_length) * 1.0u"cm" / header["HubbleParam"]), 
        ("unit_mass", parse(Float64, unit_mass) * 1.0u"g" / header["HubbleParam"]), 
        ("unit_velocity", parse(Float64, unit_velocity) * 1.0u"cm/s"), 
        ("unit_magnetic", parse(Float64, unit_magnetic) * 1.0u"Gauss")
    ])

    set_units(unit_dict)
    unit_dict
end

"""
    build_field_array_single(data_types_by_field_name::Dict,
                             cell_type::String, field_name::String,
                             file_container::FileTypes.GIZMO)

Returns an array of the field data from a single GIZMO file. If the field_name
is not in the list of available keys, then it loads the data directory from the
key given by the user. 

...
# Arguments
- `data_types_by_field_name::Dict`: Defined in FileTypes.jl.
- `cell_type::String`: For example, "gas" for gas particles.
- `field_name::String`: For example, "mass" for masses.
- `file_container::FileTypes.GIZMO`: The prepared GIZMO file container.
...
"""
function build_field_array_single(data_types_by_field_name::Dict,
                                  cell_type::String, field_name::String, 
                                  file_container::FileTypes.GIZMO)
    key = file_container.internal_key_to_file_key[cell_type]
    i = ReaderGIZMO.particle_key_to_idx(key)
    number_of_particles = header["NumPart_Total"][i]

    if number_of_particles == 0 return end

    field_transform = getfield(file_container, Symbol("$cell_type" * "_key_transform"))
    field_key = field_transform[field_name]

    if field_name == "coordinates" || field_name == "velocities"
        field_array = zeros(
            data_types_by_field_name[field_name], 
            3, 
            number_of_particles
        )
    else
        field_array = zeros(data_types_by_field_name[field_name], number_of_particles)
    end

    file_handle = h5open(file_container.file_name, "r") 
    data_from_file = read(file_handle["$key/$field_key"])
    if ndims(data_from_file) > 1 
        if field_name == "metal_mass_fraction"
            field_array = data_from_file[1, :]
        elseif field_name == "hydrogen_mass_fraction"
            field_array = 1.0 .- data_from_file[1, :] .- data_from_file[2, :]
        elseif field_name == "helium_mass_fraction"
           field_array = data_from_file[2, :]
        else
            # Coordinates/velocities can just be sent out!
            field_array = data_from_file
        end
    else
        field_array = data_from_file
    end
    close(file_handle)

    field_array
end

"""
    build_field_array_multi(data_types_by_field_name::Dict,
                            cell_type::String, field_name::String,
                            file_container::FileTypes.GIZMO)

Returns an array of the field data from multiple GIZMO files. If the field_name
is not in the list of available keys, then it loads the data directory from the
key given by the user. 

...
# Arguments
- `data_types_by_field_name::Dict`: Defined in FileTypes.jl.
- `cell_type::String`: For example, "gas" for gas particles.
- `field_name::String`: For example, "mass" for masses.
- `file_container::FileTypes.GIZMO`: The prepared GIZMO file container.
...
"""
function build_field_array_multi(data_types_by_field_name::Dict, 
                                 cell_type::String, field_name::String, 
                                 file_container::FileTypes.GIZMO)
    file_base, dummy, file_extension = split(file_container.file_name, ".")

    # Converts "gas" to "PartType0", for example
    key = file_container.internal_key_to_file_key[cell_type]
    i = ReaderGIZMO.particle_key_to_idx(key)

    # Converts "mass" to "Masses", for example
    field_transform = getfield(file_container, Symbol("$cell_type" * "_key_transform"))
    # field_key is the HDF5 file key, e.g. Masses
    field_key = field_transform[field_name]

    # Append to these instead of filling them up.
    # For now, only coordinates and velocities have different
    # shapes. 
    if field_name == "coordinates" || field_name == "velocities"
        field_array = zeros(
            data_types_by_field_name[field_name], 
            3, 
            0
        )
    else
        field_array = zeros(data_types_by_field_name[field_name], 0)
    end

    for file_number in 0:header["NumFilesPerSnapshot"]-1
        file_handle = h5open("$file_base.$file_number.$file_extension", "r")

        number_of_particles = read_attribute(file_handle["Header"], "NumPart_ThisFile")[i]

        if number_of_particles == 0 continue end

        data_from_file = read(file_handle["$key/$field_key"])
        if ndims(data_from_file) > 1
            if field_name == "coordinates" || field_name == "velocities"
                field_array = hcat(field_array, data_from_file)
            elseif field_name == "metal_mass_fraction"
                append!(field_array, data_from_file[1, :])
            elseif field_name == "hydrogen_mass_fraction"
                append!(field_array, 1.0 .- data_from_file[1, :] .- data_from_file[2, :])
            elseif field_name == "helium_mass_fraction"
                append!(field_array, data_from_file[2, :])
            end
        else
            append!(field_array, data_from_file)
        end

        close(file_handle)
    end

    field_array
end

"""
    convert_raw_field_to_units(field_array::AbstractArray, field_name::String,
                               available_fields_by_type::Dict{String, Tuple})

Converts a loaded array from the HDF5 file into the unit system from Unitful.
The units are converted based on the values read in from the provided
parameter file in the file_container::FileTypes.GIZMO. 

...
# Arguments
- `field_array::AbstractArray`: The data to convert.
- `field_name::String`: For example, "mass", "star_formation_rate", etc.
- `available_fields_by_type::Dict{String, Tuple}`: The fields that are known in Gallus.
...
"""
function convert_raw_field_to_units(field_array::AbstractArray, field_name::String,
                                    available_fields_by_type::Dict{String, Tuple})
    scale_factor = 1.0
    try
        if header["OmegaLambda"] != 0
            scale_factor = header["Time"]
        end
    catch KeyError
        # Newer versions of GIZMO have this!
        header["OmegaLambda"] = copy(header["Omega_Lambda"])
        if header["Omega_Lambda"] != 0
            scale_factor = header["Time"]
        end
    end

    conv = 1.0
    if field_name in available_fields_by_type["mass"]
        conv = units["unit_mass"]
    elseif field_name in available_fields_by_type["length"]
        conv = units["unit_length"] * scale_factor
    elseif field_name in available_fields_by_type["density"]
        conv = units["unit_mass"] / (units["unit_length"] * scale_factor)^3.0
    elseif field_name in available_fields_by_type["velocity"]
        conv = units["unit_velocity"] * sqrt(scale_factor)
    elseif field_name in available_fields_by_type["time"]
        length_conv = units["unit_length"]
        velocity_conv = units["unit_velocity"]
        conv = length_conv / velocity_conv
    elseif field_name in available_fields_by_type["specific_energy"]
        velocity_conv = units["unit_velocity"]
        conv = velocity_conv^2.0
    elseif field_name in available_fields_by_type["temperature"]
        conv = u"K"
    elseif field_name in available_fields_by_type["mass_rate"]
        mass_conv = units["unit_mass"]
        length_conv = units["unit_length"]
        velocity_conv = units["unit_velocity"]
        time_conv = length_conv / velocity_conv
        conv = mass_conv / time_conv
    elseif field_name in available_fields_by_type["unitless"]
        return field_array
    end

    field_array .* conv
end

end

"""
    read_field_from_file(file_container::FileTypes.GIZMO, cell_type::String,
                         field_name::String)

The GIZMO implementation of the function used for reading all simulation types.
Returns an array of data with Unitful units of the correct type, unless the 
requested field is not in the master list.  In that case, the returned data
is the raw HDF5 data from the file(s).

...
# Arguments
- `file_container::FileTypes.GIZMO`: The prepared GIZMO file container.
- `cell_type::String`: For example, "gas" for gas particles.
- `field_name::String`: For example, "mass" for masses.
...
"""
function read_field_from_file(file_container::FileTypes.GIZMO, cell_type::String, 
                              field_name::String)
    global available_fields_by_type, data_types_by_field_name

    # Read the header. We need this to know how many files
    # to open.
    header = ReaderGIZMO.read_header(file_container)
    units = ReaderGIZMO.get_units_from_file(file_container.parameter_file_name)

    # We will have to compute temperature ourselves!
    returning_gas_temperature = false
    if cell_type == "gas" && field_name == "temperature"
        field_name = "specific_internal_energy"
        returning_gas_temperature = true
    end

    if header["NumFilesPerSnapshot"] > 1
        field_array = ReaderGIZMO.build_field_array_multi(
            data_types_by_field_name,
            cell_type, 
            field_name, 
            file_container
        )
    else
        field_array = ReaderGIZMO.build_field_array_single(
            data_types_by_field_name,
            cell_type, 
            field_name, 
            file_container
        )
    end

    # GIZMO doesn't store temperature, so we need to compute it ourselves
    if returning_gas_temperature
        # Assume gamma = 5.0 / 3.0
        m_p_gamma_over_k_B = m_p * (5.0 / 3.0 - 1.0) / k_B
        y_helium(x) = x / (1.0 - x)
        mu(x, y) = (1.0 + 4.0 * x) / (1.0 + x + y)
        # x: helium mass fraction, y: electron abundance, z: specific internal energy
        temperature(x, y, z) = m_p_gamma_over_k_B * mu(y_helium(x), y) * z

        # TODO: Figure out a way to not have to read the file again
        # if we have already read in this information. Really, this is
        # only a problem if we need additional information to translate
        # simulation data such as translating internal energy to 
        # temperature in GIZMO.
        # 
        # To make that general would be difficult given that the ReaderMaster module
        # does not know about the Simulation module, since we need the reader first!
        helium_mass_fraction = read_field_from_file(
            file_container, 
            "gas", 
            "helium_mass_fraction"
        )
        electron_abundance = read_field_from_file(
            file_container, 
            "gas", 
            "electron_abundance"
        )
        Tools.convert_units(
            temperature.(
                helium_mass_fraction,
                electron_abundance,
                ReaderGIZMO.convert_raw_field_to_units(
                    field_array, 
                    field_name, 
                    available_fields_by_type
                )
            ),
            u"K"
        )
    else
        # Everything other than temperature is handled here.
        if field_name in available_fields_by_type["unitless"]
            field_array
        else
            ReaderGIZMO.convert_raw_field_to_units(
                field_array, 
                field_name, 
                available_fields_by_type
            )
        end
    end
end