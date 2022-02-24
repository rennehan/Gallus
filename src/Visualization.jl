module Visualization
using Distributions, Plots, PyPlot, Unitful, UnitfulAstro
using ..Halos, ..Galaxies, ..Tools
import ..Simulation

function quick_view(halo_data::Halos.HaloData, 
                    simulation::Simulation,
                    idx::Int, 
                    key::String)
    halo_position = Array{Unitful.Quantity{Float64}}(undef, 3)
    halo_position[1] = ustrip(u"kpc", halo_data.halos[idx, "x"])
    halo_position[2] = ustrip(u"kpc", halo_data.halos[idx, "y"])
    halo_position[3] = ustrip(u"kpc", halo_data.halos[idx, "z"])

    particle_idx = getfield(halo_data, Symbol("$key" * "_list"))[idx]
    if length(particle_idx) == 0
        println("No particles to show!")
    else
        key_sym = Symbol(key)
        particle_coords = ustrip.(
            u"kpc",
            Tools.build_vector_from_columns(
                Tools.convert_units(getfield(simulation, key_sym)[particle_idx, "x"], u"kpc"),
                Tools.convert_units(getfield(simulation, key_sym)[particle_idx, "y"], u"kpc"),
                Tools.convert_units(getfield(simulation, key_sym)[particle_idx, "z"], u"kpc"),
                u"kpc"
            )
        )

        view_coords = particle_coords .- halo_position
        plt.scatter(view_coords[1, :], view_coords[2, :], s = 0.5, facecolor = "k", 
                    edgecolor = "k")
    end
end

# TODO: Volume must be saved in the Simulation structure
function galaxy_stellar_mass_function(galaxy_data::Galaxies.GalaxyData, volume::Float64)
    step_size = 0.5
    sample_space = 8.0:step_size:13
    stellar_masses = log10.(ustrip.(u"Msun", galaxy_data.stars_total[!, "mass"]))
    bins = [ count(x->(lb <= x < lb+step_size), stellar_masses) for lb in sample_space]

    Plots.plot(
        sample_space .+ (step_size / 2.0), 
        bins ./ (step_size * volume), 
        yaxis=:log,
        ylim=(1e-5, 1e0)
    )
end

function black_hole_stellar_mass(galaxy_data::Galaxies.GalaxyData)
    all_bh_masses = ustrip.(u"Msun", galaxy_data.galaxies[!, "black_hole_mass"])
    all_stellar_masses = ustrip.(u"Msun", galaxy_data.galaxies[!, "bulge_mass"])

    galaxy_idx = Tools.greater_than_bitmask(all_stellar_masses, 0.0)
    if !(true in galaxy_idx)
        println("No valid galaxies!")
        return
    end

    function schutte(log_masses::AbstractArray)
        alpha = 8.80
        beta = 1.24
        scaled_stellar_masses = 10.0.^log_masses ./ 1.0e11
        
        return alpha .+ beta .* log10.(scaled_stellar_masses)
    end

    log_bh_masses = log10.(all_bh_masses[galaxy_idx])
    log_stellar_masses = log10.(all_stellar_masses[galaxy_idx])

    log_obs_stellar_masses = range(minimum(log_stellar_masses), maximum(log_stellar_masses), 100)
    log_obs_bh_masses = schutte(log_obs_stellar_masses)

    plt.plot(log_obs_stellar_masses, log_obs_bh_masses, ls = ":")
    plt.scatter(log_stellar_masses, log_bh_masses, s = 10, facecolor = "k", edgecolor = "k")
end

end