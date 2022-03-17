using FileIO, HDF5
using ..FileTypes

module ReaderGallus
using HDF5
using ...FileTypes

function read_data(next_handle::Union{HDF5.File, HDF5.Group}, data::Dict)
    for key in keys(next_handle)
        if isa(next_handle[key], HDF5.Group)
            data[key] = Dict([])
            read_data(next_handle[key], data[key])
        else
            data[key] = read(next_handle[key])
        end
    end
end

end

function read_file(file_container::FileTypes.Gallus)
    data = Dict([])
    h5open(file_container.file_name, "r") do file_handle
        ReaderGallus.read_data(file_handle, data)
    end
    data
end
