using DataFrames, CSV, StructArrays

module ReaderRockstar
import ...FileTypes
using DataFrames, Unitful, UnitfulAstro

function read_parameter_file(file_container::FileTypes.Rockstar)
    print("ReaderRockstar::read_parameter_file\n")
    parameters = Dict()
    # Need to read in the scale factor!
    open(file_container.file_name, "r") do file_handle
        for line in eachline(file_handle)
            if match(r"^#a", line) !== nothing
                parameters["a"] = split(line, "=")[2]
                break
            end
        end
    end

    open(file_container.parameter_file_name, "r") do file_handle
        for line in eachline(file_handle)
            if match(r"^#", line) !== nothing
                continue
            end
            params = split(line, "=")
            
            try
                # Remove any trailing comments after value.
                # values[1] will contain the number to parse.
                values = split(params[2])

                parameters[strip(params[1])] = strip(values[1])
            catch
                continue
            end
        end
    end

    parameters
end

function get_halo_dict_from_raw_data(raw_data::DataFrames.DataFrame)
    print("ReaderRockstar::get_halo_dict_from_raw_data\n")
    number_of_halos = nrow(raw_data)
    particle_count = Dict{String, Vector{Int64}}(
        "gas" => zeros(number_of_halos),
        "dark" => zeros(number_of_halos),
        "stars" => zeros(number_of_halos),
        "bhs" => zeros(number_of_halos)
    )
    masses = Dict{String, Vector{Float64}}(
        "mvir" => zeros(number_of_halos),
        "m200" => zeros(number_of_halos),
        "m500" => zeros(number_of_halos),
        "m2500" => zeros(number_of_halos)
    )
    radii = Dict{String, Vector{Float64}}(
        "rvir" => zeros(number_of_halos),
        "r200" => zeros(number_of_halos),
        "r500" => zeros(number_of_halos),
        "r2500" => zeros(number_of_halos)
    )

    ids = zeros(number_of_halos)
    coordinates = zeros((3, number_of_halos))
    velocities = zeros((3, number_of_halos))
    angular_momenta = zeros((3, number_of_halos))
    
    ids = raw_data[!, "#id"]

    coordinates[1, :] = raw_data[!, "x"]
    coordinates[2, :] = raw_data[!, "y"]
    coordinates[3, :] = raw_data[!, "z"]

    velocities[1, :] = raw_data[!, "vx"]
    velocities[2, :] = raw_data[!, "vy"]
    velocities[3, :] = raw_data[!, "vz"]

    angular_momenta[1, :] = raw_data[!, "Jx"]
    angular_momenta[2, :] = raw_data[!, "Jy"]
    angular_momenta[3, :] = raw_data[!, "Jz"]

    particle_count["dark"] = raw_data[!, "num_p"]

    masses["mvir"] = raw_data[!, "mvir"]
    masses["m200"] = raw_data[!, "m200c"]
    masses["m500"] = raw_data[!, "m500c"]
    masses["m2500"] = raw_data[!, "m2500c"]

    radii["rvir"] = raw_data[!, "rvir"]

    halo_dict = Dict([
        ("ids", ids),
        ("coordinates", coordinates),
        ("velocities", velocities),
        ("angular_momenta", angular_momenta),
        ("particle_count", particle_count),
        ("masses", masses),
        ("radii", radii)
    ])

    halo_dict
end

function convert_raw_data_to_units(halo_dict::Dict, params::Dict)    
    print("ReaderRockstar::convert_raw_data_to_units\n")
    output_dict = Dict()
    # Let us get rid of the Hubble parameter now and forever
    h0 = parse(Float64, params["h0"])
    a = parse(Float64, params["a"])
    mass_conv = 1.0u"Msun" / h0
    length_conv = (1.0u"Mpc" / h0) * a
    radius_conv = (1.0u"kpc" / h0) * a
    # Velocity is PHYSICAL already
    velocity_conv = 1.0u"km" / 1.0u"s"
    # Angular momentum is PHYSICAL already
    angular_momentum_conv = ((1.0u"Msun" * 1.0u"Mpc") / h0^2) * velocity_conv

    output_dict["ids"] = halo_dict["ids"]
    output_dict["coordinates"] = broadcast(.*, halo_dict["coordinates"], length_conv)
    output_dict["velocities"] = broadcast(.*, halo_dict["velocities"], velocity_conv)
    output_dict["angular_momenta"] = broadcast(.*, halo_dict["angular_momenta"], angular_momentum_conv)

    output_dict["masses"] = Dict()
    mass_keys = ["mvir", "m200", "m500", "m2500"]
    for mass_key in mass_keys
        output_dict["masses"][mass_key] = broadcast(.*, halo_dict["masses"][mass_key], mass_conv)
    end
    output_dict["radii"] = Dict()
    radius_keys = ["rvir", "r200", "r500", "r2500"]
    for radius_key in radius_keys
        output_dict["radii"][radius_key] = broadcast(.*, halo_dict["radii"][radius_key], radius_conv)
    end
    output_dict
end

end

# Rockstar ASCII is pretty simple to read since we can just
# use a CSV reader package directly into a data frame.
# All that is left is to convert the raw data into the 
# default Gallus format for halos.
#
# We will return the number of halos found, and 0 galaxies.
# There are zero galaxies because this function is only
# called when the data hasn't been post-processed yet.
# There should be a FileTypes.Gallus file saved somewhere,
# or else we must recompute it!
#
# The conversion from ASCII column to Gallus.Halos.Halo format is:
# 1  : halo_id
# 2  : particle_count[1]  (dark matter particle count)
# 3  : masses['mvir']
# 5  : radii['rvir']
# 9  : position[1]
# 10 : position[2]
# 11 : position[3]
# 12 : velocity[1]
# 13 : velocity[2]
# 14 : velocity[3]
# 15 : angular_momentum[1]
# 16 : angular_momentum[2]
# 17 : angular_momentum[3]
# 29 : masses['m200']
# 30 : masses['m500']
# 31 : masses['m2500']

# All other fields must be computed, or ignored until implemented.
function read_file(file_container::FileTypes.Rockstar)
    raw_data = DataFrame(CSV.File(file_container.file_name, header=1, skipto=20))
    params = ReaderRockstar.read_parameter_file(file_container)
    halo_dict = ReaderRockstar.get_halo_dict_from_raw_data(raw_data)
    halo_dict = ReaderRockstar.convert_raw_data_to_units(copy(halo_dict), params)
    halo_dict
end