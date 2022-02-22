module FileTypes

# Abstract types allow us to not have to know
# which file type is being used. We simply use
# the type to determine how we process the data.
#
# Each implemented file format needs its own
# struct in order to work properly in the
# ReaderMaster and WriterMaster modules.
export File, Rockstar, AHF
using Parameters

abstract type File end

mutable struct Gallus <: File
    file_name::String
    parameter_file_name::String
end

mutable struct Rockstar <: File
    file_name::String
    parameter_file_name::String
end

mutable struct AHF <: File
    file_name::String
    parameter_file_name::String
end

mutable struct GIZMO <: File
    file_name::String
    parameter_file_name::String
end

end