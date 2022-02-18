import Base.show, Base.print
using DataFrames
mutable struct Galaxies
    obj::Galaxies
    galaxies::DataFrames.DataFrame

    """
        Galaxies(file_container::FileTypes.Gallus)

    Load in a previously saved galaxy list.

    ...
    # Arguments
    - `file_container::FileTypes.Gallus`: The file name struct for the galaxy data.
    ...
    """
    function Galaxies(file_container::FileTypes.Gallus)

    end

    """
        Galaxies(halos::Halos)

    Create a new galaxy list from a given set of Gallus Halos.

    ...
    # Arguments
    - `halos::Halos`: The master Halos object loaded by Gallus.
    ...
    """
    function Galaxies(halos::Halos)
        galaxies = new()
        galaxies.data = DataFrame()
        galaxies.gas = DataFrame()
        galaxies.stars = DataFrame()
        galaxies.dark = DataFrame()
        galaxies.bhs = DataFrame()

        num_halos = nrow(halos.halos)
        coordinates = halos[!, "coordinates"]
        velocities = halos[!, "velocities"]
        sfrs = zeros(num_halos)
        masses = Dict{String, Array{Float64}}([
            "gas" => zeros(num_halos),
            "stars" => zeros(num_halos),
            "dark" => zeros(num_halos),
            "bhs" => zeros(num_halos),
            "dust" => zeros(num_halos),
            "HI" => zeros(num_halos),
            "H2" => zeros(num_halos),
            "bulge" => zeros(num_halos)
        ])
        half_mass_radii = Dict{String, Array{Float64}}([
            "gas" => zeros(num_halos),
            "stars" => zeros(num_halos),
            "dark" => zeros(num_halos),
            "bhs" => zeros(num_halos),
            "baryon" => zeros(num_halos),
            "total" => zeros(num_halos)
        ])
        twenty_percent_radii = Dict{String, Array{Float64}}([
            "gas" => zeros(num_halos),
            "stars" => zeros(num_halos),
            "dark" => zeros(num_halos),
            "bhs" => zeros(num_halos),
            "baryon" => zeros(num_halos),
            "total" => zeros(num_halos)
        ])
        eighty_percent_radii = Dict{String, Array{Float64}}([
            "gas" => zeros(num_halos),
            "stars" => zeros(num_halos),
            "dark" => zeros(num_halos),
            "bhs" => zeros(num_halos),
            "baryon" => zeros(num_halos),
            "total" => zeros(num_halos)
        ])
        metallicities = Dict{String, Array{Float64}}([
            "gas_mass_weighted" => zeros(num_halos),
            "gas_sfr_weighted" => zeros(num_halos),
            "stars" => zeros(num_halos)
        ])
        velocity_dispersions = Dict{String, Array{Float64}}([
            "gas" => zeros(num_halos),
            "stars" => zeros(num_halos),
            "dark" => zeros(num_halos),
            "bhs" => zeros(num_halos),
            "baryon" => zeros(num_halos),
            "total" => zeros(num_halos)
        ])
        angular_momenta = Dict{String, Array{Float64}}([
            "gas" => zeros(num_halos),
            "stars" => zeros(num_halos),
            "dark" => zeros(num_halos),
            "bhs" => zeros(num_halos),
            "baryon" => zeros(num_halos),
            "total" => zeros(num_halos)
        ])
        Threads.@threads for i in 1:num_halos
            
        end

        galaxies
    end
end
show(io::IO, h::Galaxies) = show(io, "Created Galaxies structure!")