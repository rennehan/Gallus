module Tools
using Unitful, UnitfulAstro

function idx_to_axis(idx::Int)
    if idx == 1
        "x"
    elseif idx == 2
        "y"
    else
        "z"
    end
end

function idx_to_axis(idx::Int, type::String)
    if idx == 1
        "$type" * "x"
    elseif idx == 2
        "$type" * "y"
    else
        "$type" * "z"
    end
end

function build_vector_from_columns(x::AbstractArray, 
                                   y::AbstractArray, 
                                   z::AbstractArray, 
                                   unit::Unitful.FreeUnits)
    new_array = zeros(3, length(x)) .* unit

    new_array[1, :] = x
    new_array[2, :] = y
    new_array[3, :] = z

    new_array
end

function centered_cube_bitmask(coords::Array{Float64}, radius::Float64)
    # Save some time by only selecting those within the cube
    x = ifelse.(abs.(@view(coords[1, :])) .< radius, true, false)
    y = ifelse.(abs.(@view(coords[2, :])) .< radius, true, false)
    z = ifelse.(abs.(@view(coords[3, :])) .< radius, true, false)
    x .&& y .&& z
end

convert_units(data::AbstractArray, unit::Unitful.FreeUnits) = ustrip.(unit, data) .* unit
convert_units(data::Number, unit::Unitful.FreeUnits) = ustrip(unit, data) * unit
vector_norm_squared(vector) = @views @inbounds [sum(abs2, vector[:, j]) for j=1:size(vector)[2]]
less_than_bitmask(a::AbstractArray, b::Number) = ifelse.(a .< b, true, false)
between_bitmask(a::AbstractArray, b::Number, c::Number) = ifelse.(b .< a .< c, true, false)
greater_than_bitmask(a::AbstractArray, b::Number) = ifelse.(a .> b, true, false)

function shift_coordinates(a::AbstractArray, b::AbstractArray)
    a .= a .- b
end

function weighted_sum(a::AbstractArray, b::AbstractArray)
    ifelse(sum(a) != 0.0, sum(a .* b) / sum(a), 0.0)
end

@views function cross(a::AbstractArray, b::AbstractArray, comp::Int)
    if comp == 1
        a[2, :].*b[3, :] .- a[3, :].*b[2, :]
    elseif comp == 2
        a[3, :].*b[1, :] .- a[1, :].*b[3, :]
    else
        a[1, :].*b[2, :] .- a[2, :].*b[1, :]
    end
end

function get_cumulative_property(property::Vector{Float64}, 
                                 coords::Matrix{Float64}, 
                                 radius::Float64)
    sum(property[
        Tools.less_than_bitmask(
            Tools.vector_norm_squared(coords),
            radius^2
        )
    ])
end
function get_cumulative_property(property::Vector{Float64},
                                 radii_squared::Vector{Float64},
                                 radius::Float64)
    sum(property[
        Tools.less_than_bitmask(
            radii_squared,
            radius^2
        )
    ])
end
function get_shell_property(property::Vector{Float64},
                            coords::Matrix{Float64},
                            radius_start::Float64,
                            radius_end::Float64)
    sum(property[
        Tools.between_bitmask(
            Tools.vector_norm_squared(coords),
            radius_start^2,
            radius_end^2
        )
    ])
end
function get_shell_property(property::Vector{Float64},
                            radii_squared::Vector{Float64},
                            radius_start::Float64,
                            radius_end::Float64)
    sum(property[
        Tools.between_bitmask(
            radii_squared,
            radius_start^2,
            radius_end^2
        )
    ])
end

end