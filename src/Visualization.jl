module Visualization
using Distributions, Plots, PyPlot, Unitful, UnitfulAstro
using ..Halos, ..Galaxies, ..Tools, ..Simulation, ..Observations

function quick_view(halo_data::Halos.HaloData, 
                    sim_data::Simulation.SimulationData,
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
        particle_coords = ustrip.(
            u"kpc",
            Simulation.get_field(sim_data, key, "coordinates"),
        )

        view_coords = particle_coords .- halo_position
        plt.scatter(view_coords[1, :], view_coords[2, :], s = 0.5, facecolor = "k", 
                    edgecolor = "k")
    end
end

# TODO: Volume must be saved in the Simulation structure.
# TODO: Make this better!
function galaxy_stellar_mass_function(galaxy_data::Galaxies.GalaxyData, redshift::Float64, 
                                      author::String, volume::Float64)
    step_size = 0.5
    sample_space = 6:step_size:13
    stellar_masses = log10.(ustrip.(u"Msun", galaxy_data.stars_total[!, "mass"]))
    bins = [count(x->(lb <= x < lb+step_size), stellar_masses) for lb in sample_space]
    sm_points = sample_space .+ (step_size / 2.0)

    plt.figure(figsize = (7, 6), dpi = 120)
    plt.yscale("log")
    plt.ylim([1.0e-5, 1.0e0])
    plt.plot(sm_points, Observations.gsmf(sm_points, redshift, author), label = "Observations")
    plt.plot(
        sm_points,
        bins ./ (step_size * volume), 
        label = "Simulation")
    plt.legend(loc = "upper right")
    plt.show()
end

function black_hole_stellar_mass(galaxy_data::Galaxies.GalaxyData)
    all_bh_masses = ustrip.(u"Msun", galaxy_data.galaxies[!, "black_hole_mass"])
    all_stellar_masses = ustrip.(u"Msun", galaxy_data.galaxies[!, "bulge_mass"])

    galaxy_idx = Tools.greater_than_bitmask(all_stellar_masses, 0.0)
    if !(true in galaxy_idx)
        println("No valid galaxies!")
        return
    end

    log_bh_masses = log10.(all_bh_masses[galaxy_idx])
    log_stellar_masses = log10.(all_stellar_masses[galaxy_idx])

    log_obs_stellar_masses = range(
        minimum(log_stellar_masses), 
        maximum(log_stellar_masses), 
        100
    )

    plt.figure(figsize = (7, 6), dpi = 120)
    plt.plot(
        log_obs_stellar_masses, 
        Observations.schutte18(log_obs_stellar_masses), 
        ls = ":", 
        lw = 3, 
        label = "Schutte+'18"
    )
    plt.scatter(
        log_stellar_masses, 
        log_bh_masses, 
        s = 10, 
        lw = 3,
        facecolor = "k", 
        edgecolor = "k",
        label = "Simulation"
    )
    plt.legend(loc = "upper left")
    plt.show()
end

end