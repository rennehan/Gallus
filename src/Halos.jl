import Base.show, Base.print
using DataFrames
mutable struct Halos{T <: FileTypes.File}
    obj::Halos
    halo_files::T
    halos::DataFrames.DataFrame
    subhalos::Ref
    ids::Vector{UInt64}
    masses::Dict
    radii::Dict
    coordinates::Matrix
    velocities::Matrix
    angular_momenta::Matrix
    velocity_dispersions::Dict
    # Contains all of the possible calculations one might want to do eventually
    derived::Dict
    function Halos(halo_files::T) where T
        halos = new{T}()
        halo_dict = Reader.read_file(halo_files)
        halos.halos = DataFrame()
        halos.halos[!, "ids"] = halo_dict["ids"]

        mass_keys = ["mvir", "m200", "m500", "m2500"]
        for mass_key in mass_keys
            halos.halos[!, mass_key] = halo_dict["masses"][mass_key]
        end
        radius_keys = ["rvir", "r200", "r500", "r2500"]
        for radius_key in radius_keys
            halos.halos[!, radius_key] = halo_dict["radii"][radius_key]
        end

        halos.halos[!, "x"] = halo_dict["coordinates"][1, :]
        halos.halos[!, "y"] = halo_dict["coordinates"][2, :]
        halos.halos[!, "z"] = halo_dict["coordinates"][3, :]

        halos.halos[!, "vx"] = halo_dict["velocities"][1, :]
        halos.halos[!, "vy"] = halo_dict["velocities"][2, :]
        halos.halos[!, "vz"] = halo_dict["velocities"][3, :]

        halos.halos[!, "jx"] = halo_dict["angular_momenta"][1, :]
        halos.halos[!, "jy"] = halo_dict["angular_momenta"][2, :]
        halos.halos[!, "jz"] = halo_dict["angular_momenta"][3, :]

        halos
    end
end
show(io::IO, h::Halos) = show(io, "Created Halos structure!")
