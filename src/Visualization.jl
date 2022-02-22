module Visualization
using Distributions, Plots, Unitful, UnitfulAstro
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
        Plots.plot(view_coords[1, :], view_coords[2, :], seriestype = :scatter)
    end
end

function galaxy_stellar_mass_function(galaxy_data::Galaxies.GalaxyData, volume::Float64)
    step_size = 0.5
    sample_space = 8.0:step_size:13
    stellar_masses = log10.(ustrip.(u"Msun", galaxy_data.stars_total[!, "mass"]))
    bins = [ count(x->(lb <= x < lb+step_size), stellar_masses) for lb in sample_space]

    plot(
        sample_space .+ (step_size / 2.0), 
        bins ./ (step_size * volume), 
        yaxis=:log,
        ylim=(1e-5, 1e0)
    )
end

end