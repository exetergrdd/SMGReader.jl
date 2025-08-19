"""
    mutable struct VectorBuffer{T} <: AbstractVector{T}

A resizable buffer that wraps a `Vector{T}` and tracks a logical length (`len`).
Efficient for repeated resizing and reuse, e.g., when processing records in a loop.

    Recommendation is to allocate a large enough vector buffer to work with, e.g.
    following constants are set in moddata.jl

    const MOD_6MA_BUFFER_LENGTH = 70_000
    const MOD_5MC_BUFFER_LENGTH = 70_000
    const READ_BUFFER_LENGTH = 200_000
    const ALIGN_BLOCK_LENGTH = 5_000
"""
mutable struct VectorBuffer{T} <: AbstractVector{T}
    data::Vector{T}
    len::Int
end

"""
    VectorBuffer{T}(backcapacity::Int, initcapacity=0)

Construct a `VectorBuffer` with backing storage of size `backcapacity` and logical length `initcapacity`.
"""
function VectorBuffer{T}(backcapacity::Int, initcapacity=0) where T
    VectorBuffer(Vector{T}(undef, backcapacity), initcapacity)
end

### abstract vector interface
@inline Base.size(buff::VectorBuffer) = (buff.len,)
@inline Base.axes(buff::VectorBuffer) = (Base.OneTo(buff.len),)
@inline Base.eachindex(buf::VectorBuffer) = Base.OneTo(buf.len)
@inline Base.IndexStyle(::Type{<:VectorBuffer}) = IndexLinear()

"""
    setlength!(buff::VectorBuffer, newlen, backfactor=1.2)

Resize the logical length of the buffer to `newlen`. If the backing storage is too small,
resize it by a factor of `backfactor` or double the current length, whichever is larger.
Returns the buffer.
"""
@inline function setlength!(buff::VectorBuffer, newlen, backfactor=1.2)
    if newlen > length(buff.data)
        newbufferlen = max(2*buff.len, Int(round(newlen*backfactor)))
        resize!(buff.data, newbufferlen)
    end
    buff.len = newlen
    buff
end

"""
    checkbounds(buff::VectorBuffer, i)

Check if index `i` is within the logical bounds of the buffer.
Throws `BoundsError` if out of bounds.
"""
@inline function checkbounds(buff::VectorBuffer, i)
    if 1 <= i <= buff.len
        return true
    else
        throw(BoundsError(buff, i))
    end
end

"""
    Base.getindex(buff::VectorBuffer{T}, i::Int)

Return the element at logical index `i` in the buffer.
"""
@inline function Base.getindex(buff::VectorBuffer{T}, i::Int) where {T}
    @boundscheck checkbounds(buff, i)
    @inbounds buff.data[i]
end

"""
    Base.setindex!(buff::VectorBuffer{T}, val, i::Int)

Set the element at logical index `i` to `val`.
"""
@inline function Base.setindex!(buf::VectorBuffer{T}, val, i::Int) where {T}
    @boundscheck checkbounds(buf, i)
    @inbounds buf.data[i] = val
end

"""
    Base.getindex(buff::VectorBuffer{T}, inds::AbstractVector{Bool})

Return a regular array of elements where `inds` is `true`.
Logical indexing must match the logical length of the buffer.
"""
@inline function Base.getindex(buff::VectorBuffer{T}, inds::AbstractVector{Bool}) where {T}
    @boundscheck @assert length(inds) == buff.len "Logical index must match buffer length"
    return buff.data[1:buff.len][inds]  # Returns regular array
end

"""
    Base.getindex(buff::VectorBuffer{T}, r::UnitRange{Int})

Return a slice of the buffer for the given range `r`.
"""
@inline function Base.getindex(buff::VectorBuffer{T}, r::UnitRange{Int}) where {T}
    @boundscheck checkbounds(buff, first(r)) && checkbounds(buff, last(r)) 
    return buff.data[r]
end
