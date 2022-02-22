module Galaxies
import Base.show, Base.print
using DataFrames, Unitful, UnitfulAstro, ..FileTypes, ..Halos, ..Tools
using Plots, ProgressMeter, Statistics, LinearAlgebra
import ..Simulation

mutable struct GalaxyData
    obj::GalaxyData
    galaxies::DataFrames.DataFrame
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
        GalaxyData(halo_data::Halos.HaloData, simulation::Simulation)

    Create a new galaxy list from a given set of Gallus Halos.

    ...
    # Arguments
    - `halo_data::Halos.HaloData`: The master HaloData object loaded by Gallus.
    - `simulation::Simulation`: The master Simulation structure.
    ...
    """
    function GalaxyData(halo_data::Halos.HaloData, simulation::Simulation)
        galaxy_data = new()
        galaxy_data.galaxies = DataFrame()

        num_halos = nrow(halo_data.halos)
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

        half_mass_radii = Dict{String, Array{Float64}}([])
        total_mass_radii = Dict{String, Array{Float64}}([])
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
            half_mass_radii[key] = zeros(Float64, num_halos)
            total_mass_radii[key] = zeros(Float64, num_halos)
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
        for key in keys
            all_coords[key] = ustrip.(
                u"kpc",
                Tools.build_vector_from_columns(
                    Tools.convert_units(
                        getfield(simulation, Symbol(key))[!, "x"], 
                        u"kpc"
                    ),
                    Tools.convert_units(
                        getfield(simulation, Symbol(key))[!, "y"], 
                        u"kpc"
                    ),
                    Tools.convert_units(
                        getfield(simulation, Symbol(key))[!, "z"], 
                        u"kpc"
                    ),
                    u"kpc"
                )
            )

            all_velocities[key] = ustrip.(
                u"km/s",
                Tools.build_vector_from_columns(
                    Tools.convert_units(
                        getfield(simulation, Symbol(key))[!, "vx"], 
                        u"km/s"
                    ),
                    Tools.convert_units(
                        getfield(simulation, Symbol(key))[!, "vy"], 
                        u"km/s"
                    ),
                    Tools.convert_units(
                        getfield(simulation, Symbol(key))[!, "vz"], 
                        u"km/s"
                    ),
                    u"km/s"
                )
            )

            all_masses[key] = ustrip.(
                u"Msun",
                getfield(simulation, Symbol(key))[!, "mass"]
            )

            if key == "gas" || key == "stars"
                all_metallicities[key] = getfield(simulation, Symbol(key))[!, "metallicity"]

                if key == "gas"
                    all_sfrs[key] = ustrip.(u"Msun/yr", simulation.gas[!, "star_formation_rate"])
                end
            end
        end

        function get_cumulative_masses(mass::Vector{Float64}, 
                                       coords::Matrix{Float64}, 
                                       radius::Float64)
            sum(mass[
                Tools.less_than_bitmask(
                    Tools.vector_norm_squared(coords),
                    radius^2
                )
            ])
        end
        function get_cumulative_masses(mass::Vector{Float64},
                                       radii_squared::Vector{Float64},
                                       radius::Float64)
            sum(mass[
                Tools.less_than_bitmask(
                    radii_squared,
                    radius^2
                )
            ])
        end
        function get_shell_mass(mass::Vector{Float64},
                                coords::Matrix{Float64},
                                radius_start::Float64,
                                radius_end::Float64)
            sum(mass[
                Tools.between_bitmask(
                    Tools.vector_norm_squared(coords),
                    radius_start^2,
                    radius_end^2
                )
            ])
        end
        function get_shell_mass(mass::Vector{Float64},
                                radii_squared::Vector{Float64},
                                radius_start::Float64,
                                radius_end::Float64)
            sum(mass[
                Tools.between_bitmask(
                    radii_squared,
                    radius_start^2,
                    radius_end^2
                )
            ])
        end
        function check_done(a::Float64, b::Float64, tol::Float64)
            if a != 0.0
                (abs(a - b) / a) < tol
            else
                true
            end
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
                try
                    if length(getfield(halo_data, Symbol("$key" * "_list"))[i]) < 64
                        # If there are not enough particles, just remove it from the list
                        push!(to_delete, key)
                    end
                catch
                    # If there is no list, just ignore?
                    push!(to_delete, key)
                end
            end

            if length(to_delete) == length(keys)
                continue
            end
            for key in to_delete
                deleteat!(halo_particle_keys, findall(x->x==key,halo_particle_keys))
            end

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
            cumulative_mass = Dict{String, Vector{Float64}}([])
            done_flags = Dict{String, Bool}([])
            halo_half_mass_radii = Dict{String, Float64}([])
            half_total_mass = Dict{String, Float64}([])
            half_mass_idx = Dict{String, Int}([])
            total_mass_idx = Dict{String, Int}([])
            for key in halo_particle_keys
                cumulative_mass[key] = zeros(Float64, number_of_bins) # Msun
                done_flags[key] = false
                halo_half_mass_radii[key] = 0.0  # kpc
                half_total_mass[key] = 0.0  # Msun
                half_mass_idx[key] = 1
                total_mass_idx[key] = 1

                cumulative_mass[key][1] = get_cumulative_masses(
                    masses_in_halo[key],
                    particle_radii2[key],
                    radius
                )
                if cumulative_mass[key][1] == 0.0
                    cumulative_mass[key][1] = tolerance
                end
            end

            # Keep an index to access the current mass values at radius
            j = 2
            while false in values(done_flags)
                for key in halo_particle_keys
                    if done_flags[key]
                        continue
                    end

                    shell_mass = get_shell_mass(masses_in_halo[key], particle_radii2[key],
                                                radius, radius + delta_r)
                    cumulative_mass[key][j] = cumulative_mass[key][j - 1] + shell_mass

                    if check_done(cumulative_mass[key][j], cumulative_mass[key][j - 1], tolerance)
                        done_flags[key] = true
                    end
                    if radius > max_radius
                        done_flags[key] = true
                    end
                end
                
                radii[j] = radius
                radius += delta_r
                j += 1
            end

            for key in halo_particle_keys
                @inbounds for j=1:number_of_bins
                    k = (number_of_bins - j) + 1
                    if cumulative_mass[key][k] != 0.0
                        half_total_mass[key] = cumulative_mass[key][k] / 2.0
                        total_mass_idx[key] = k
                        break
                    end
                end

                @inbounds for j=1:number_of_bins
                    if cumulative_mass[key][j] > half_total_mass[key]
                        half_mass_idx[key] = j - 1
                        break
                    end
                end

                # No mass?
                if half_mass_idx[key] == 0 continue end

                half_mass_radii[key][i] = mean([
                    radii[half_mass_idx[key]], 
                    radii[half_mass_idx[key]] + 1
                ])
                total_mass_radii[key][i] = mean([
                    radii[total_mass_idx[key]],
                    radii[total_mass_idx[key]] - 1
                ])

                # Now compute all of the properties within the half mass radius
                # and our definition of "total galaxy", 3*Rhalf.
                #
                # TODO: Remove hard-coded number 3.0
                half_idx = Tools.less_than_bitmask(
                    particle_radii2[key],
                    half_mass_radii[key][i]                   
                )

                # No mass?
                if !(true in half_idx) continue end

                full_idx = Tools.less_than_bitmask(
                    particle_radii2[key],
                    3.0 * half_mass_radii[key][i]                 
                )
                masses_half[key][i] = half_total_mass[key]
                masses_total[key][i] = cumulative_mass[key][total_mass_idx[key]]
                veldisp_half[key][i] = std(vel_norms[key][half_idx])
                veldisp_total[key][i] = std(vel_norms[key][full_idx])

                @inbounds for j=1:3
                    cross_half_j = Tools.cross(
                        offset_coords[key][:, half_idx], 
                        vels_in_halo[key][:, half_idx], 
                        j
                    )
                    cross_total_j = Tools.cross(
                        offset_coords[key][:, full_idx], 
                        vels_in_halo[key][:, full_idx], 
                        j
                    )
                    angular_momenta_half[key][j, i] = sum(cross_half_j)
                    angular_momenta_total[key][j, i] = sum(cross_total_j)
                end

                if key == "gas"
                    sfr_half[key][i] = sum(sfrs_in_halo[key][half_idx])
                    sfr_total[key][i] = sum(sfrs_in_halo[key][full_idx])

                    z_half["z_sfr"][i] = Tools.weighted_sum(
                        sfrs_in_halo[key][half_idx], 
                        z_in_halo[key][half_idx]
                    )
                    z_total["z_sfr"][i] = Tools.weighted_sum(
                        sfrs_in_halo[key][full_idx], 
                        z_in_halo[key][full_idx]
                    )

                    z_half["z_mass"][i] = Tools.weighted_sum(
                        masses_in_halo[key][half_idx], 
                        z_in_halo[key][half_idx]
                    )
                    z_total["z_mass"][i] = Tools.weighted_sum(
                        masses_in_halo[key][full_idx], 
                        z_in_halo[key][full_idx]
                    )
                end
                if key == "stars"
                    z_half[key][i] = median(z_in_halo[key][half_idx])
                    z_total[key][i] = median(z_in_halo[key][full_idx])
                end
            end

            ProgressMeter.next!(p)
        end

        galaxy_data.galaxies[!, "half_mass_radius"] = half_mass_radii["stars"] .* u"kpc"
        galaxy_data.galaxies[!, "radius"] = 3.0 .* half_mass_radii["stars"] .* u"kpc"

        # I do not see a way to set DataFrame columns using the setfield()
        # method, so I have to do it manually at this point.
        for key in keys
            half_df = DataFrame()
            total_df = DataFrame()

            half_df[!, "mass"] = masses_half[key] .* u"Msun"
            half_df[!, "radius"] = half_mass_radii[key] .* u"kpc"
            half_df[!, "velocity_dispersion"] = veldisp_half[key] .* u"km/s"
            half_df[!, "jx"] = angular_momenta_half[key][1, :] .* u"kpc * km/s"
            half_df[!, "jy"] = angular_momenta_half[key][2, :] .* u"kpc * km/s"
            half_df[!, "jz"] = angular_momenta_half[key][3, :] .* u"kpc * km/s"

            total_df[!, "mass"] = masses_total[key] .* u"Msun"
            total_df[!, "radius"] = total_mass_radii[key] .* u"kpc"
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

        galaxy_data
    end
end
show(io::IO, g::GalaxyData) = show(io, "Created GalaxyData structure!")

end