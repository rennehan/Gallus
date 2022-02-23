using FileIO
using ...FileTypes

module WriterGallus
using FileIO
using ...FileTypes

end

function save_file(file_container::FileTypes.Gallus, data::Dict)
    save(file_container.file_name, data)
end