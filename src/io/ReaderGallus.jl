using FileIO
using ..FileTypes

module ReaderGallus
using ...FileTypes

end

function read_file(file_container::FileTypes.Gallus)
    load(file_container.file_name)
end