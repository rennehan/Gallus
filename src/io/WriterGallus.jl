using FileIO, HDF5
using ...FileTypes

module WriterGallus
using FileIO, HDF5
using ...FileTypes

function save_data(next_handle::Union{HDF5.File,HDF5.Group}, next_data::Dict)
    for key in keys(next_data)
        key_str = string(key)
        if !isa(next_data[key], Dict)
            next_handle[key_str] = next_data[key]
        else
            create_group(next_handle, key_str)
            save_data(next_handle[key_str], next_data[key])
        end
    end
end

end

function save_file(file_container::FileTypes.Gallus, data::Dict)
    h5open(file_container.file_name, "w") do file_handle
        WriterGallus.save_data(file_handle, data)
    end
end
