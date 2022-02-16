import Base.show, Base.print
mutable struct Simulation{T <: FileTypes.File}
    obj::Simulation
    simulation_files::T
    gas::Dict
    stars::Dict
    dark::Dict
    bhs::Dict
    function Simulation(simulation_files::T) where T
        simulation = new{T}()
        particle_dict = Reader.read_file(simulation_files)
        simulation.gas = Dict()
        simulation.stars = Dict()
        simulation.dark = Dict()
        simulation.bhs = Dict()
        
        simulation.gas = particle_dict["gas"]
        simulation.stars = particle_dict["stars"]
        simulation.dark = particle_dict["dark"]
        simulation.bhs = particle_dict["bhs"]
        simulation
    end
end
show(io::IO, s::Simulation) = show(io, "Created Simulation structure!")
