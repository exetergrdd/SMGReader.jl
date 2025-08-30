##### Constants for initial buffer sizes
const MOD_6MA_BUFFER_LENGTH = 70_000
const MOD_5MC_BUFFER_LENGTH = 70_000
const READ_BUFFER_LENGTH    = 200_000
const ALIGN_BLOCK_LENGTH    = 5_000

### types and structures for preallocating and handling readdata

abstract type HTSReadData end

"""
    StencillingData(auxmap::AuxMap) 

Construct preallocated HTS read data building a read to genome map for chromatin stencilling setting an BAM auxillary data map `AuxMap`.
e.g.

  - For generic tracking of modifications `StencillingData(AuxMapMod())`
  - For tracking FIRE fields `StencillingData(AuxMapModFire())`
"""
struct StencillingData{T} <: HTSReadData
    alignmap::VectorBuffer{Int}
    auxmap::T
end
StencillingData(am::T) where {T} = StencillingData{T}(VectorBuffer{Int}(READ_BUFFER_LENGTH), am)

"""
    DirectRNA(auxmap::AuxMap)

Construct preallocated HTS read data building a read to genome map for DirectRNA setting an BAM auxillary data map `AuxMap`.
e.g.

  - For generic tracking of modifications `DirectRNA(AuxMapMod())`
"""
struct DirectRNA{T} <: HTSReadData
    alignmap::VectorBuffer{Int}
    auxmap::T
end
DirectRNA(am::T) where {T} = DirectRNA{T}(VectorBuffer{Int}(READ_BUFFER_LENGTH), am)


"""
    DirectRNAAlignBlocks(auxmap::AuxMap)

Construct preallocated HTS read data building both a read to genome map and an align blocks to genome map for spliced alignments. 
e.g.

  - For generic tracking of modifications `DirectRNAAlignBlocks(AuxMapMod())`
"""
struct DirectRNAAlignBlocks{T} <: HTSReadData
    alignmap::VectorBuffer{Int}
    alignblocks::VectorBuffer{AlignBlock}
    auxmap::T
end
DirectRNAAlignBlocks(am::T) where {T} = DirectRNAAlignBlocks{T}(VectorBuffer{Int}(READ_BUFFER_LENGTH), VectorBuffer{AlignBlock}(ALIGN_BLOCK_LENGTH), am)



"""
    processread!(record::BamRecord, recorddata::HTSReadData)

Process read mapping auxillary fields and building read to genome map using preallocated `recorddata`.

"""
@inline function processread!(record::BamRecord, recorddata::HTSReadData)
    readtogenome!(record, recorddata.alignmap)
    auxmap!(record, recorddata.auxmap)
end

@inline function processread!(r, rdata::DirectRNAAlignBlocks{T}) where {T}
    readtogenome!(r, rdata.alignmap)
    alignmentblocks!(r, rdata.alignblocks)
    auxmap!(r, rdata.auxmap)
end


### accessor functions
"""
    haplotype(record::BamRecord, recorddata::HTSReadData)

Return haplotype of record and nothing if haplotype missing
"""
@inline haplotype(record::BamRecord, recorddata::HTSReadData, default=nothing) = isnothing(rdata.auxmap.hp) ? default : r.data[rdata.auxmap.hp.start]
@inline hashaplotype(r, rdata) = isnothing(rdata.auxmap.hp)