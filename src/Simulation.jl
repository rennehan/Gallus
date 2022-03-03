module Simulation
import Base.show, Base.print
using DataFrames, Unitful, UnitfulAstro
using ..Tools, ..FileTypes, ..Reader

mutable struct SimulationData{T <: FileTypes.File}
    obj::SimulationData
    simulation_files::T
    cell_data::Dict
    loaded_fields::Dict
    function SimulationData(simulation_files::T) where T
        sim_data = new{T}()
        sim_data.simulation_files = simulation_files
        sim_data.cell_data = Dict()
        keys = ["gas", "stars", "dark", "bhs"]
        for key in keys
            sim_data.cell_data[key] = DataFrame()
        end

        sim_data.loaded_fields = Dict{String, Bool}([])
        sim_data
    end
end
show(io::IO, s::SimulationData) = show(io, "Created SimulationData structure!")

function field_present(sim_data::SimulationData, cell_type::String, field_name::String)
    ifelse("$cell_type-$field_name" in keys(sim_data.loaded_fields), true, false)
end

function get_field(sim_data::SimulationData, cell_type::String, field_name::String)
    if field_present(sim_data, cell_type, field_name)
        if field_name == "coordinates"
            Tools.build_vector_from_columns(
                sim_data.cell_data[cell_type][!, "x"],
                sim_data.cell_data[cell_type][!, "y"],
                sim_data.cell_data[cell_type][!, "z"],
                u"kpc"
            )
        elseif field_name == "velocities"
            Tools.build_vector_from_columns(
                sim_data.cell_data[cell_type][!, "vx"],
                sim_data.cell_data[cell_type][!, "vy"],
                sim_data.cell_data[cell_type][!, "vz"],
                u"km/s"
            )
        else
            # cell_data is a dictionary with keys "gas", "stars", "dark", "bhs"
            # with each element containing a DataFrame
            sim_data.cell_data[cell_type][!, field_name] 
        end
    else
        data_from_file = Reader.read_field_from_file(
            sim_data.simulation_files, 
            cell_type, 
            field_name
        )
        sim_data.loaded_fields["$cell_type-$field_name"] = true
        if field_name == "coordinates"
            for i=1:3
                sim_data.cell_data[cell_type][!, Tools.idx_to_axis(i)] = data_from_file[i, :]
            end

            Tools.build_vector_from_columns(
                sim_data.cell_data[cell_type][!, "x"],
                sim_data.cell_data[cell_type][!, "y"],
                sim_data.cell_data[cell_type][!, "z"],
                u"kpc"
            )
        elseif field_name == "velocities"
            for i=1:3
                sim_data.cell_data[cell_type][!, Tools.idx_to_axis(i, "v")] = data_from_file[i, :]
            end

            Tools.build_vector_from_columns(
                sim_data.cell_data[cell_type][!, "vx"],
                sim_data.cell_data[cell_type][!, "vy"],
                sim_data.cell_data[cell_type][!, "vz"],
                u"km/s"
            )
        else
            println(cell_type)
            println(field_name)
            sim_data.cell_data[cell_type][!, field_name] = data_from_file
            sim_data.cell_data[cell_type][!, field_name]
        end
    end
end

end