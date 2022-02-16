import Base.show, Base.print
using DataFrames
mutable struct Simulation{T <: FileTypes.File}
    obj::Simulation
    simulation_files::T
    gas::DataFrames.DataFrame
    stars::DataFrames.DataFrame
    dark::DataFrames.DataFrame
    bhs::DataFrames.DataFrame
    function Simulation(simulation_files::T) where T
        simulation = new{T}()
        particle_dict = Reader.read_file(simulation_files)
        simulation.gas = DataFrame()
        simulation.stars = DataFrame()
        simulation.dark = DataFrame()
        simulation.bhs = DataFrame()
        for key in keys(particle_dict["gas"])
            if key == "coordinates"
                simulation.gas[!, "x"] = particle_dict["gas"][key][1, :]
                simulation.gas[!, "y"] = particle_dict["gas"][key][2, :]
                simulation.gas[!, "z"] = particle_dict["gas"][key][3, :]
            elseif key == "velocities"
                simulation.gas[!, "vx"] = particle_dict["gas"][key][1, :]
                simulation.gas[!, "vy"] = particle_dict["gas"][key][2, :]
                simulation.gas[!, "vz"] = particle_dict["gas"][key][3, :]
            else
                simulation.gas[!, key] = particle_dict["gas"][key]
            end
        end
        for key in keys(particle_dict["stars"])
            if key == "coordinates"
                simulation.stars[!, "x"] = particle_dict["stars"][key][1, :]
                simulation.stars[!, "y"] = particle_dict["stars"][key][2, :]
                simulation.stars[!, "z"] = particle_dict["stars"][key][3, :]
            elseif key == "velocities"
                simulation.stars[!, "vx"] = particle_dict["stars"][key][1, :]
                simulation.stars[!, "vy"] = particle_dict["stars"][key][2, :]
                simulation.stars[!, "vz"] = particle_dict["stars"][key][3, :]
            else
                simulation.stars[!, key] = particle_dict["stars"][key]
            end
        end
        for key in keys(particle_dict["dark"])
            if key == "coordinates"
                simulation.dark[!, "x"] = particle_dict["dark"][key][1, :]
                simulation.dark[!, "y"] = particle_dict["dark"][key][2, :]
                simulation.dark[!, "z"] = particle_dict["dark"][key][3, :]
            elseif key == "velocities"
                simulation.dark[!, "vx"] = particle_dict["dark"][key][1, :]
                simulation.dark[!, "vy"] = particle_dict["dark"][key][2, :]
                simulation.dark[!, "vz"] = particle_dict["dark"][key][3, :]
            else
                simulation.dark[!, key] = particle_dict["dark"][key]
            end
        end
        for key in keys(particle_dict["bhs"])
            if key == "coordinates"
                simulation.bhs[!, "x"] = particle_dict["bhs"][key][1, :]
                simulation.bhs[!, "y"] = particle_dict["bhs"][key][2, :]
                simulation.bhs[!, "z"] = particle_dict["bhs"][key][3, :]
            elseif key == "velocities"
                simulation.bhs[!, "vx"] = particle_dict["bhs"][key][1, :]
                simulation.bhs[!, "vy"] = particle_dict["bhs"][key][2, :]
                simulation.bhs[!, "vz"] = particle_dict["bhs"][key][3, :]
            else
                simulation.bhs[!, key] = particle_dict["bhs"][key]
            end
        end
        simulation
    end
end
show(io::IO, s::Simulation) = show(io, "Created Simulation structure!")
