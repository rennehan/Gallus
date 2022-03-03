module Galaxies
import Base.show, Base.print
using DataFrames, Unitful, UnitfulAstro, ..FileTypes, ..Halos, ..Tools
using Plots, ProgressMeter, Statistics, LinearAlgebra
using ..Simulation

mutable struct GalaxyData
    obj::GalaxyData
    galaxies::DataFrames.DataFrame
    total_galaxies::UInt32
    galaxy_bitmask::BitVector

    gas_total::DataFrames.DataFrame
    stars_total::DataFrames.DataFrame
    dark_total::DataFrames.DataFrame
    bhs_total::DataFrames.DataFrame

    gas_half::DataFrames.DataFrame
    stars_half::DataFrames.DataFrame
    dark_half::DataFrames.DataFrame
    bhs_half::DataFrames.DataFrame
    """
        GalaxyData(file_container::FileTypes.Gallus)

    Load in a previously saved galaxy list.

    ...
    # Arguments
    - `file_container::FileTypes.Gallus`: The file name struct for the galaxy data.
    ...
    """
    function GalaxyData(file_container::FileTypes.Gallus)

    end

    """
        GalaxyData(halo_data::Halos.HaloData, sim_data::Simulation.SimulationData)

    Create a new galaxy list from a given set of Gallus Halos.

    ...
    # Arguments
    - `halo_data::Halos.HaloData`: The master HaloData object loaded by Gallus.
    - `sim_data::Simulation.SimulationData`: The master Simulation structure.
    ...
    """
    function GalaxyData(halo_data::Halos.HaloData, sim_data::Simulation.SimulationData)
        galaxy_data = new()
        galaxy_data.galaxies = DataFrame()

        num_halos = nrow(halo_data.halos)

        # "stars" /must/ be the first key because all of the quantities
        # are calculated within Rhalf,* and 3Rhalf,*
        keys = ["gas", "stars", "dark", "bhs"]

        @inbounds for i=1:3
            ax = Tools.idx_to_axis(i)
            galaxy_data.galaxies[!, ax] = Tools.convert_units(halo_data.halos[!, ax], u"kpc")
        end
        @inbounds for i=1:3
            vax = Tools.idx_to_axis(i, "v")
            galaxy_data.galaxies[!, vax] = Tools.convert_units(halo_data.halos[!, vax], u"km/s")
        end

        galaxy_data.galaxies[!, "half_mass_radius"] = zeros(Float64, num_halos) .* u"kpc"
        galaxy_data.galaxies[!, "radius"] = zeros(Float64, num_halos) .* u"kpc"
        galaxy_data.galaxies[!, "bulge_mass"] = zeros(Float64, num_halos) .* u"Msun"
        galaxy_data.galaxies[!, "black_hole_mass"] = zeros(Float64, num_halos) .* u"Msun"

        bulge_masses = zeros(Float64, num_halos)
        black_hole_masses = zeros(Float64, num_halos)
        half_mass_radius = zeros(Float64, num_halos)
        total_mass_radius = zeros(Float64, num_halos)
        masses_half = Dict{String, Array{Float64}}([])
        masses_total = Dict{String, Array{Float64}}([])
        veldisp_half = Dict{String, Array{Float64}}([])
        veldisp_total = Dict{String, Array{Float64}}([])
        angular_momenta_half = Dict{String, Matrix{Float64}}([])
        angular_momenta_total = Dict{String, Matrix{Float64}}([])
        sfr_half = Dict{String, Array{Float64}}([])
        sfr_total = Dict{String, Array{Float64}}([])
        save_keys = ["gas", "stars", "dark", "bhs", "baryon", "total"]
        for key in save_keys
            masses_half[key] = zeros(Float64, num_halos)
            masses_total[key] = zeros(Float64, num_halos)
            veldisp_half[key] = zeros(Float64, num_halos)
            veldisp_total[key] = zeros(Float64, num_halos)
            angular_momenta_half[key] = zeros(Float64, 3, num_halos)
            angular_momenta_total[key] = zeros(Float64, 3, num_halos)
            
            if key == "gas"
                sfr_half[key] = zeros(Float64, num_halos)
                sfr_total[key] = zeros(Float64, num_halos)
            end
        end
        
        z_half = Dict{String, Array{Float64}}([])
        z_total = Dict{String, Array{Float64}}([])
        save_keys = ["z_mass", "z_sfr", "stars"]
        for key in save_keys
            z_half[key] = zeros(Float64, num_halos)
            z_total[key] = zeros(Float64, num_halos)
        end

        # These contain the ustrip versions of the data frame columns
        # for all particles
        all_metallicities = Dict{String, Vector{Float64}}([])
        all_coords = Dict{String, Matrix{Float64}}([])
        all_masses = Dict{String, Vector{Float64}}([])
        all_velocities = Dict{String, Matrix{Float64}}([])
        all_sfrs = Dict{String, Vector{Float64}}([])
        all_bh_masses = Dict{String, Vector{Float64}}([])
        for key in keys
            all_coords[key] = ustrip.(
                u"kpc",
                Simulation.get_field(sim_data, key, "coordinates")
            )

            all_velocities[key] = ustrip.(
                u"km/s",
                Simulation.get_field(sim_data, key, "velocities")
            )

            all_masses[key] = ustrip.(
                u"Msun",
                Simulation.get_field(sim_data, key, "mass")
            )

            if key == "gas" || key == "stars"
                all_metallicities[key] = Simulation.get_field(
                    sim_data,
                    key,
                    "metal_mass_fraction"
                )

                if key == "gas"
                    all_sfrs[key] = ustrip.(
                        u"Msun/yr", 
                        Simulation.get_field(sim_data, "gas", "star_formation_rate")
                    )
                end
            end
            if key == "bhs"
                all_bh_masses[key] = ustrip.(
                    u"Msun", 
                    Simulation.get_field(sim_data, "bhs", "subgrid_mass")
                )
            end
        end

        function check_done(a::Float64, b::Float64, tol::Float64)
            ifelse(a != 0.0, (abs(a - b) / a) < tol, true)
        end

        p = ProgressMeter.Progress(num_halos, desc = "Finding galaxy properties...")
        # By this point we should have all of the particles in each halo.
        @views @inbounds Threads.@threads for i in 1:num_halos
            # Some halos might not have particles of a certain type.
            # In that case, we only do the try..catch statement
            # once and store a boolean for whether we should ignore it
            # in all future loops for this halo.
            halo_particle_keys = copy(keys)
            to_delete = []
            for key in halo_particle_keys
                if length(getfield(halo_data, Symbol("$key" * "_list"))[i]) < 1
                    # If there are not enough particles, just remove it from the list
                    push!(to_delete, key)
                end
            end

            if length(to_delete) == length(keys) continue end
            for key in to_delete
                deleteat!(halo_particle_keys, findall(x->x==key,halo_particle_keys))
            end

            # A "galaxy" must have stars
            if !("stars" in halo_particle_keys) continue end

            halo_rvir = ustrip(Float64, u"kpc", halo_data.halos[i, "rvir"])
            halo_position = Array{Float64}(undef, 3)
            halo_position[1] = ustrip(Float64, u"kpc", halo_data.halos[i, "x"])
            halo_position[2] = ustrip(Float64, u"kpc", halo_data.halos[i, "y"])
            halo_position[3] = ustrip(Float64, u"kpc", halo_data.halos[i, "z"])

            offset_coords = Dict{String, Matrix{Float64}}([])
            masses_in_halo = Dict{String, Vector{Float64}}([])
            vels_in_halo = Dict{String, Matrix{Float64}}([])
            particle_radii2 = Dict{String, Vector{Float64}}([])
            vel_norms = Dict{String, Vector{Float64}}([])
            # TODO: Needs to be a matrix?
            z_in_halo = Dict{String, Vector{Float64}}([])
            sfrs_in_halo = Dict{String, Vector{Float64}}([])
            bh_masses_in_halo = Dict{String, Vector{Float64}}([])
            for key in halo_particle_keys
                # E.g., slice_idx = halo_data.gas_list[i]
                slice_idx = getfield(halo_data, Symbol("$key" * "_list"))[i]
                offset_coords[key] = all_coords[key][:, slice_idx] .- halo_position
                masses_in_halo[key] = all_masses[key][slice_idx]
                vels_in_halo[key] = all_velocities[key][:, slice_idx]
                particle_radii2[key] = Tools.vector_norm_squared(offset_coords[key])
                vel_norms[key] = sqrt.(Tools.vector_norm_squared(vels_in_halo[key]))

                if key == "gas" || key == "stars"
                    z_in_halo[key] = all_metallicities[key][slice_idx]

                    if key == "gas"
                        sfrs_in_halo[key] = all_sfrs[key][slice_idx]
                    end
                end
                if key == "bhs"
                    bh_masses_in_halo[key] = all_bh_masses[key][slice_idx]
                end
            end          

            tolerance = Float64(0.01)  # % tolerance for a change in mass

            # In kpc
            # TODO: Make starting radius dependent on softening scale.
            delta_r = Float64(0.5)
            radius = Float64(1.0)
            max_radius = halo_rvir

            # Store more information than necessary right now
            # so that we don't have to push!() to all of the radial
            # vectors. The final value is the last element that
            # is nonzero.
            number_of_bins = floor(Int, max_radius / delta_r)
            radii = zeros(Float64, number_of_bins)
            cumulative_mass = zeros(Float64, number_of_bins)
            done_flag = false
            cumulative_mass[1] = Tools.get_cumulative_property(
                masses_in_halo["stars"],
                particle_radii2["stars"],
                radius
            )
            radii[1] = radius
            radius += delta_r

            # Keep an index to access the current mass values at radius
            j = 2
            while !done_flag
                shell_mass = Tools.get_shell_property(
                    masses_in_halo["stars"], 
                    particle_radii2["stars"],
                    radius, 
                    radius + delta_r
                )
                cumulative_mass[j] = cumulative_mass[j - 1] + shell_mass

                # We just integrate all the way out to max_radius
                # for the SMBHs since they are quite scattered
                # throughout the halo.
                done_flag = check_done(
                    cumulative_mass[j], 
                    cumulative_mass[j - 1], 
                    tolerance
                )

                if radius > max_radius
                    done_flag = true
                end
                
                radii[j] = radius
                radius += delta_r
                j += 1
            end

            half_total_mass = 0.0
            total_mass_idx = 1
            half_mass_idx = 0
            # Find the stellar total and half radii indices
            @inbounds for j=1:number_of_bins
                k = (number_of_bins - j) + 1
                if cumulative_mass[k] != 0.0
                    half_total_mass = cumulative_mass[k] / 2.0
                    total_mass_idx = k
                    break
                end
            end

            if half_total_mass == 0.0 continue end

            @inbounds for j=1:number_of_bins
                if cumulative_mass[j] > half_total_mass
                    half_mass_idx = j - 1
                    break
                end
            end

            # First bin could have all of the mass. Therefore,
            # won't be able to access the 0th index (j-1). 
            # When that happens, just store 0.5 * full_radius
            # as the half mass radius, and don't keep track of 
            # any of the Rhalf properties.
            if half_mass_idx != 0
                half_mass_radius[i] = mean([
                    radii[half_mass_idx], 
                    radii[half_mass_idx] + 1
                ])
            else
                half_mass_radius[i] = radii[1] / 2.0
            end

            # Now compute all of the properties within the half mass radius
            # and our definition of "total galaxy", 3*Rhalf.
            #
            total_mass_radius[i] = mean([
                radii[total_mass_idx],
                radii[total_mass_idx] - 1
            ])
            
            for key in halo_particle_keys
                # TODO: Remove hard-coded number 3.0
                full_idx = Tools.less_than_bitmask(
                    particle_radii2[key],
                    (3.0 * half_mass_radius[i])^2.0           
                )
                half_idx = Tools.less_than_bitmask(
                    particle_radii2[key],
                    half_mass_radius[i]^2.0               
                )

                # Sometimes there could be no particles in Rhalf
                if true in half_idx
                    masses_half[key][i] = sum(masses_in_halo[key][half_idx])

                    true_count = 0
                    key_map = Array{Int}(undef, 0)
                    @inbounds for j=1:length(half_idx) 
                        if half_idx[j] 
                            true_count += 1 
                            push!(key_map, j)
                        end 
                    end

                    # std() gives NaN by default for 1-element instead
                    # of crashing
                    if true_count > 1
                        veldisp_half[key][i] = std(vel_norms[key][half_idx])
                    end
                    center_of_mass = zeros(Float64, 3)
                    @views @inbounds for j=1:3
                        center_of_mass[j] = Tools.weighted_sum(
                            masses_in_halo[key][half_idx],
                            offset_coords[key][j, half_idx]
                        )
                    end

                    # The angular momentum for each component is with
                    # respect to the center of mass of each component.
                    Tools.shift_coordinates(
                        offset_coords[key][:, half_idx],
                        center_of_mass
                    )

                    cross_half = zeros(Float64, 3, true_count)
                    @views @inbounds for j=1:3
                        cross_half[j, :] = Tools.cross(
                            offset_coords[key][:, half_idx], 
                            vels_in_halo[key][:, half_idx], 
                            j
                        )
                        angular_momenta_half[key][j, i] = sum(cross_half[j])
                    end

                    Tools.shift_coordinates(
                        offset_coords[key][:, half_idx],
                        -1.0 .* center_of_mass
                    )

                    if key == "gas"
                        sfr_half[key][i] = sum(sfrs_in_halo[key][half_idx])

                        z_half["z_sfr"][i] = Tools.weighted_sum(
                            sfrs_in_halo[key][half_idx], 
                            z_in_halo[key][half_idx]
                        )

                        z_half["z_mass"][i] = Tools.weighted_sum(
                            masses_in_halo[key][half_idx], 
                            z_in_halo[key][half_idx]
                        )
                    end
                    if key == "stars"
                        z_half[key][i] = median(z_in_halo[key][half_idx])

                        # Determine how much stellar mass is counter-rotating
                        # to determine the bulge mass.
                        # M_bulge ~ 2 * M_counter (<Rhalf)
                        counter_rotating_mass = 0.0
                        @views @inbounds for j=1:size(cross_half)[2]
                            if sum(angular_momenta_half[key][:, i] .* cross_half[:, j]) < 0
                                counter_rotating_mass += masses_in_halo[key][key_map[j]]
                            end
                        end

                        bulge_masses[i] = 2.0 * counter_rotating_mass
                    end
                    if key == "bhs"
                        core_idx = Tools.less_than_bitmask(
                            particle_radii2[key],
                            1.0
                        )
                        if true in core_idx
                            black_hole_masses[i] = maximum(bh_masses_in_halo[key][core_idx])
                        end
                    end
                end

                if true in full_idx
                    masses_total[key][i] = sum(masses_in_halo[key][full_idx])

                    # std() gives NaN by default for 1 element instead of 
                    # crashing.
                    if length(vel_norms[key][full_idx]) > 1
                        veldisp_total[key][i] = std(vel_norms[key][full_idx])
                    end

                    center_of_mass = zeros(Float64, 3)
                    @views @inbounds for j=1:3
                        center_of_mass[j] = Tools.weighted_sum(
                            masses_in_halo[key][full_idx],
                            offset_coords[key][j, full_idx]
                        )
                    end
                    # The angular momentum for each component is with
                    # respect to the center of mass of each component.
                    Tools.shift_coordinates(
                        offset_coords[key][:, full_idx],
                        center_of_mass
                    )

                    @views @inbounds for j=1:3
                        cross_total_j = Tools.cross(
                            offset_coords[key][:, full_idx], 
                            vels_in_halo[key][:, full_idx], 
                            j
                        )
                        angular_momenta_total[key][j, i] = sum(cross_total_j)
                    end

                    Tools.shift_coordinates(
                        offset_coords[key][:, full_idx],
                        -1.0 .* center_of_mass
                    )

                    if key == "gas"
                        sfr_total[key][i] = sum(sfrs_in_halo[key][full_idx])

                        z_total["z_sfr"][i] = Tools.weighted_sum(
                            sfrs_in_halo[key][full_idx], 
                            z_in_halo[key][full_idx]
                        )

                        z_total["z_mass"][i] = Tools.weighted_sum(
                            masses_in_halo[key][full_idx], 
                            z_in_halo[key][full_idx]
                        )
                    end
                    if key == "stars"
                        z_total[key][i] = median(z_in_halo[key][full_idx])
                    end
                end
            end

            ProgressMeter.next!(p)
        end

        galaxy_data.galaxies[!, "half_mass_radius"] = half_mass_radius .* u"kpc"
        galaxy_data.galaxies[!, "radius"] = 3.0 .* half_mass_radius .* u"kpc"
        galaxy_data.galaxies[!, "bulge_mass"] = bulge_masses .* u"Msun"
        galaxy_data.galaxies[!, "black_hole_mass"] = black_hole_masses .* u"Msun"

        # I do not see a way to set DataFrame columns using the setfield()
        # method, so I have to do it manually at this point.
        for key in keys
            half_df = DataFrame()
            total_df = DataFrame()

            half_df[!, "mass"] = masses_half[key] .* u"Msun"
            half_df[!, "radius"] = half_mass_radius .* u"kpc"
            half_df[!, "velocity_dispersion"] = veldisp_half[key] .* u"km/s"
            half_df[!, "jx"] = angular_momenta_half[key][1, :] .* u"kpc * km/s"
            half_df[!, "jy"] = angular_momenta_half[key][2, :] .* u"kpc * km/s"
            half_df[!, "jz"] = angular_momenta_half[key][3, :] .* u"kpc * km/s"

            total_df[!, "mass"] = masses_total[key] .* u"Msun"
            total_df[!, "radius"] = total_mass_radius .* u"kpc"
            total_df[!, "velocity_dispersion"] = veldisp_total[key] .* u"km/s"
            total_df[!, "jx"] = angular_momenta_total[key][1, :] .* u"kpc * km/s"
            total_df[!, "jy"] = angular_momenta_total[key][2, :] .* u"kpc * km/s"
            total_df[!, "jz"] = angular_momenta_total[key][3, :] .* u"kpc * km/s"

            if key == "gas"
                half_df[!, "sfr"] = sfr_half[key]
                total_df[!, "sfr"] = sfr_total[key]

                half_df[!, "metallicity_sfr"] = z_half["z_sfr"]
                half_df[!, "metallicity_mass"] = z_half["z_mass"]

                total_df[!, "metallicity_sfr"] = z_total["z_sfr"]
                total_df[!, "metallicity_mass"] = z_total["z_mass"]
            end
            if key == "stars"
                half_df[!, "metallicity"] = z_half["stars"]
                total_df[!, "metallicity"] = z_total["stars"]
            end

            setfield!(galaxy_data, Symbol("$key" * "_half"), half_df)
            setfield!(galaxy_data, Symbol("$key" * "_total"), total_df)
        end

        galaxy_data.galaxy_bitmask = Tools.greater_than_bitmask(
            galaxy_data.stars_total[!, "mass"],
            0.0 * u"Msun"
        )
        galaxy_data.total_galaxies = length(
            galaxy_data.stars_total[
                galaxy_data.galaxy_bitmask,
                "mass"
            ]
        )

        galaxy_data
    end
end
show(io::IO, g::GalaxyData) = show(io, "Created GalaxyData structure!")

end