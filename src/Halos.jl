module Halos
using DataFrames, Unitful, UnitfulAstro, ProgressMeter
using ..FileTypes, ..Reader
import ..Simulation

mutable struct HaloData{T <: FileTypes.File}
    obj::HaloData
    halo_files::T
    halos::DataFrames.DataFrame
    gas_list::Dict
    stars_list::Dict
    dark_list::Dict
    bhs_list::Dict
    galaxies::Ref
    subhalos::Ref
    function HaloData(halo_files::T) where T
        halo_data = new{T}()
        halo_dict = Reader.read_file(halo_files)
        halo_data.halos = DataFrame()
        halo_data.halos[!, "ids"] = halo_dict["ids"]
        halo_data.gas_list = Dict{UInt32, Array{UInt32}}([])
        halo_data.stars_list = Dict{UInt32, Array{UInt32}}([])
        halo_data.dark_list = Dict{UInt32, Array{UInt32}}([])
        halo_data.bhs_list = Dict{UInt32, Array{UInt32}}([]) 
        for i in nrow(halo_data.halos)
            halo_data.gas_list[i] = Array{UInt32}(undef, 0)
            halo_data.stars_list[i] = Array{UInt32}(undef, 0)
            halo_data.dark_list[i] = Array{UInt32}(undef, 0)
            halo_data.bhs_list[i] = Array{UInt32}(undef, 0)
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

"""
    Halos.find_particles_in_halo(halo_data::HaloData, simulation::Simulation)

Takes all of the particle coordinates to determine which particles are in which 
halos. The main loop is threaded so that it should be reasonably fast.

...
# Arguments
- `halo_data::HaloData`: The previously constructed HaloData struct.
- `simulation::Simulation`: The main simulation structure containing particle data.
...
"""
function find_particles_in_halo(halo_data::HaloData, simulation::Simulation)
    function find_particles_in_cube(coords::Array{Float64}, 
                                    radius::Float64)
        # Save some time by only selecting those within the cube
        x = ifelse.(abs.(@view(coords[1, :])) .< radius, true, false)
        y = ifelse.(abs.(@view(coords[2, :])) .< radius, true, false)
        z = ifelse.(abs.(@view(coords[3, :])) .< radius, true, false)
        x .&& y .&& z
    end

    function generate_particle_lists(all_particles_mask::BitVector, sphere_mask::BitVector)
        idx_list = Array{UInt32}(undef, 0)
        # running_idx corresponds to the index in XYZ_sphere_idx
        # since everything is aligned in memory
        running_idx = 1
        @inbounds for j in 1:length(all_particles_mask)
            if all_particles_mask[j] === true
                if sphere_mask[running_idx] === true
                    push!(idx_list, j)
                end
                running_idx += 1
            end
        end
        idx_list
    end

    gas_coords = Array{Unitful.Quantity{Float64}}(undef, 3, length(simulation.gas[!, "x"]))
    star_coords = Array{Unitful.Quantity{Float64}}(undef, 3, length(simulation.stars[!, "x"]))
    dark_coords = Array{Unitful.Quantity{Float64}}(undef, 3, length(simulation.dark[!, "x"]))
    bh_coords = Array{Unitful.Quantity{Float64}}(undef, 3, length(simulation.bhs[!, "x"]))

    function conv(var)
        if var == 1
            "x"
        elseif var == 2
            "y"
        else
            "z"
        end
    end

    for i in 1:3
        gas_coords[i, :] = simulation.gas[!, conv(i)]
        star_coords[i, :] = simulation.stars[!, conv(i)]
        dark_coords[i, :] = simulation.dark[!, conv(i)]
        bh_coords[i, :] = simulation.bhs[!, conv(i)]
    end
    
    gas_coords = ustrip.(Float64, u"kpc", gas_coords)
    star_coords = ustrip.(Float64, u"kpc", star_coords)
    dark_coords = ustrip.(Float64, u"kpc", dark_coords)
    bh_coords = ustrip.(Float64, u"kpc", bh_coords)

    r2(coords) = @views @inbounds [sum(abs2, coords[:, j]) for j=1:size(coords)[2]]
    within_r2(radii2, rvir2) = ifelse.(radii2 .< rvir2, true, false)

    # TODO: Load balance. Probably randomly shuffling indices is good enough.
    num_halos = nrow(halo_data.halos)
    @views @inbounds Threads.@threads for i in 1:num_halos
        halo_position = Array{Float64}(undef, 3)
        halo_position[1] = ustrip(Float64, u"kpc", halo_data.halos[i, "x"])
        halo_position[2] = ustrip(Float64, u"kpc", halo_data.halos[i, "y"])
        halo_position[3] = ustrip(Float64, u"kpc", halo_data.halos[i, "z"])

        rvir = ustrip(Float64, u"kpc", halo_data.halos[i, "rvir"])
        rvir_squared = rvir^2

        offset_gas_coords = gas_coords .- halo_position
        if Threads.threadid() == 1
            println(halo_position)
            println("i=$i")
        end
        gas_idx = find_particles_in_cube(offset_gas_coords, rvir)
        if true in gas_idx
            gas_sphere_idx = within_r2(r2(offset_gas_coords[:, gas_idx]), rvir_squared)
            if true in gas_sphere_idx
                halo_data.gas_list[i] = generate_particle_lists(gas_idx, gas_sphere_idx)
            end
        end
        offset_star_coords = star_coords .- halo_position
        stars_idx = find_particles_in_cube(offset_star_coords, rvir)
        if true in stars_idx
            stars_sphere_idx = within_r2(r2(offset_star_coords[:, stars_idx]), rvir_squared)
            if true in stars_sphere_idx
                halo_data.stars_list[i] = generate_particle_lists(stars_idx, stars_sphere_idx)
            end
        end
        offset_dark_coords = dark_coords .- halo_position
        dark_idx = find_particles_in_cube(offset_dark_coords, rvir)
        if true in dark_idx
            dark_sphere_idx = within_r2(r2(offset_dark_coords[:, dark_idx]), rvir_squared)
            if true in dark_sphere_idx
                halo_data.dark_list[i] = generate_particle_lists(dark_idx, dark_sphere_idx)
            end
        end
        offset_bh_coords = bh_coords .- halo_position
        bhs_idx = find_particles_in_cube(offset_bh_coords, rvir)
        if true in bhs_idx
            bhs_sphere_idx = within_r2(r2(offset_bh_coords[:, bhs_idx]), rvir_squared)
            if true in bhs_sphere_idx
                halo_data.bhs_list[i] = generate_particle_lists(bhs_idx, bhs_sphere_idx)
            end
        end
    end
end

end
