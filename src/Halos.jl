module Halos
import Base.show, Base.print
using DataFrames, Unitful, UnitfulAstro, ProgressMeter, NearestNeighbors
using ..FileTypes, ..Reader, ..Writer, ..Tools, ..Simulation

mutable struct HaloData{T <: FileTypes.File}
    obj::HaloData
    halo_files::T
    halos::DataFrames.DataFrame
    gas_list::Dict{UInt32, Vector{Int64}}
    stars_list::Dict{UInt32, Vector{Int64}}
    dark_list::Dict{UInt32, Vector{Int64}}
    bhs_list::Dict{UInt32, Vector{Int64}}
    galaxies::Ref
    subhalos::Ref
    function HaloData(halo_files::T) where T
        halo_data = new{T}()
        halo_dict = Reader.read_file(halo_files)
        halo_data.halos = DataFrame()
        halo_data.halos[!, "id"] = halo_dict["ids"]
        halo_data.gas_list = Dict([])
        halo_data.stars_list = Dict([])
        halo_data.dark_list = Dict([])
        halo_data.bhs_list = Dict([]) 
        num_halos = nrow(halo_data.halos)
        for i in 1:num_halos
            halo_data.gas_list[i] = Vector{Int64}(undef, 0)
            halo_data.stars_list[i] = Vector{Int64}(undef, 0)
            halo_data.dark_list[i] = Vector{Int64}(undef, 0)
            halo_data.bhs_list[i] = Vector{Int64}(undef, 0)
        end

        mass_keys = ["mvir", "m200", "m500", "m2500"]
        for mass_key in mass_keys
            halo_data.halos[!, mass_key] = halo_dict["masses"][mass_key]
        end
        radius_keys = ["rvir", "r200", "r500", "r2500"]
        for radius_key in radius_keys
            halo_data.halos[!, radius_key] = halo_dict["radii"][radius_key]
        end

        halo_data.halos[!, "x"] = halo_dict["coordinates"][1, :]
        halo_data.halos[!, "y"] = halo_dict["coordinates"][2, :]
        halo_data.halos[!, "z"] = halo_dict["coordinates"][3, :]

        halo_data.halos[!, "vx"] = halo_dict["velocities"][1, :]
        halo_data.halos[!, "vy"] = halo_dict["velocities"][2, :]
        halo_data.halos[!, "vz"] = halo_dict["velocities"][3, :]

        halo_data.halos[!, "jx"] = halo_dict["angular_momenta"][1, :]
        halo_data.halos[!, "jy"] = halo_dict["angular_momenta"][2, :]
        halo_data.halos[!, "jz"] = halo_dict["angular_momenta"][3, :]

        halo_data
    end
end
show(io::IO, h::HaloData) = show(io, "Created HaloData structure!")

function save_particles_in_halo(halo_data::HaloData, file_container::FileTypes.Gallus)
    keys = ["gas", "stars", "dark", "bhs"]
    data = Dict{String, Dict{UInt32, Vector{Int64}}}([])
    for key in keys
        data[key] = getfield(halo_data, Symbol("$key" * "_list"))
    end

    Writer.save_file(file_container, data)
end

"""
    Halos.find_particles_in_halo(halo_data::HaloData, sim_data::Simulation.SimulationData)

Takes all of the particle coordinates to determine which particles are in which 
halos. The main loop is threaded so that it should be reasonably fast.

...
# Arguments
- `halo_data::HaloData`: The previously constructed HaloData struct.
- `sim_data::Simulation.SimulationData`: The main simulation structure containing particle data.
- `file_container::FileTypes.Gallus`: File container for Gallus files to read/save halos.
...
"""
function find_particles_in_halo(halo_data::HaloData, sim_data::Simulation.SimulationData, 
                                file_container::FileTypes.Gallus)
    if isfile(file_container.file_name)
        data = Reader.read_file(file_container)
        keys = ["gas", "stars", "dark", "bhs"]
        for key in keys
            setfield!(halo_data, Symbol("$key" * "_list"), data[key]) 
        end

        return
    end

    keys = ["gas", "stars", "dark", "bhs"]
    for key in keys
        particle_tree = Simulation.get_tree(sim_data, key)

        # TODO: Load balance. Probably randomly shuffling indices is good enough.
        num_halos = nrow(halo_data.halos)
        p = ProgressMeter.Progress(num_halos, 0.5, "Finding $key particles in halos...")
        @views @inbounds Threads.@threads for i in 1:num_halos
            halo_position = Array{Float64}(undef, 3)
            halo_position[1] = ustrip(Float64, u"kpc", halo_data.halos[i, "x"])
            halo_position[2] = ustrip(Float64, u"kpc", halo_data.halos[i, "y"])
            halo_position[3] = ustrip(Float64, u"kpc", halo_data.halos[i, "z"])

            rvir = ustrip(Float64, u"kpc", halo_data.halos[i, "rvir"])

            idxs = inrange(particle_tree, halo_position, rvir, true)

            # TODO: Can probably get rid of this ifelse tree if I am smart
            if key == "gas"
                halo_data.gas_list[i] = idxs
            elseif key == "stars"
                halo_data.stars_list[i] = idxs
            elseif key == "dark"
                halo_data.dark_list[i] = idxs
            else
                halo_data.bhs_list[i] = idxs
            end

            ProgressMeter.next!(p)
        end
    end

    save_particles_in_halo(halo_data, file_container)
end

end
