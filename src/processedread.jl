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
@inline haplotype(record::BamRecord, recorddata::HTSReadData, default=nothing) = isnothing(recorddata.auxmap.hp) ? default : record.data[recorddata.auxmap.hp.start]

"""
    hashaplotype(record::BamRecord, recorddata::HTSReadData)

Return true if `record` has haplotype field
"""
@inline hashaplotype(record::BamRecord, recorddata::HTSReadData) = isnothing(recorddata.auxmap.hp)



### functions to map to genome
@inline function genomecoords(pos, len, r::BamRecord, rdata; onebased=true)
    if ispositive(r)
        return rdata.alignmap[pos] + onebased, rdata.alignmap[pos + len - 1] + onebased
    else
        return rdata.alignmap[r.core.l_qseq - (pos + len - 1) + 1] + onebased, rdata.alignmap[r.core.l_qseq - pos + 1] + onebased
    end
end

### accessor functions
@inline firenucpos(r, am::AuxMapModFire) = reinterpret(Int32, @view(r.data[am.ns.start:am.ns.stop]))
@inline firenuclen(r, am::AuxMapModFire) = reinterpret(Int32, @view(r.data[am.nl.start:am.nl.stop]))
@inline firemsppos(r, am::AuxMapModFire) = reinterpret(Int32, @view(r.data[am.as.start:am.as.stop]))
@inline firemsplen(r, am::AuxMapModFire) = reinterpret(Int32, @view(r.data[am.al.start:am.al.stop]))
@inline firemspqual(r, am::AuxMapModFire) = @view r.data[am.aq.start:am.aq.stop]


@inline firemsps(r, rdata) = zip(firemsppos(r, rdata.auxmap), firemsplen(r, rdata.auxmap), firemspqual(r, rdata.auxmap))
@inline firenucs(r, rdata) = zip(firenucpos(r, rdata.auxmap), firenuclen(r, rdata.auxmap))


@inline haplotype(r, rdata, default=nothing) = isnothing(rdata.auxmap.hp) ? default : r.data[rdata.auxmap.hp.start]
@inline hashaplotype(r, rdata) = isnothing(rdata.auxmap.hp)



@inline firegenomeelement(r, rdata) = x -> (genomecoords(x[1],x[2], r, rdata)..., x[3:end]...)
@inline firefilt(x, onebased=true) = (x[1] > onebased) && (x[2] > onebased)
@inline firegenome(r, rdata) = it -> Iterators.filter(firefilt, Iterators.map(firegenomeelement(r, rdata), it))