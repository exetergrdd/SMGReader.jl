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

@inline function processread!(record::BamRecord, recorddata::DirectRNAAlignBlocks)
    readtogenome!(record, recorddata.alignmap)
    alignmentblocks!(record, recorddata.alignblocks)
    auxmap!(record, recorddata.auxmap)
end


### accessor functions
"""
    rightposition(recorddata::HTSReadData)

Get the right position of alignment on genome.
"""
@inline rightposition(recorddata::HTSReadData) = recorddata.alignmap[findlast(!iszero, recorddata.alignmap)]


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
"""
    genomecoords(pos::Int, record::BamRecord, recorddata::HTSReadData; onebased=true)

Map read coordinates to genome coordates that uses `alignmap` in `recorddata`, `onebased=true` for 1-based coordinates
"""
@inline function genomecoords(pos::Int, record::BamRecord, recorddata::HTSReadData; onebased=true)
    if ispositive(r)
        return recorddata.alignmap[pos] + onebased
    else
        return recorddata.alignmap[record.core.l_qseq - (pos + len - 1) + 1] + onebased
    end
end


"""
    genomecoords(pos::Int, len::Int, record::BamRecord, recorddata::HTSReadData; onebased=true)

Map read coordinates to genome coordates that uses `alignmap` in `recorddata`,

  - `pos` is one-based to match the alignment map
  - `one-based = true`  returns one based inclusive coords
  - `one-based = false` returns zero based exclusive coords

`onebased=true` for 1-based coordinates
"""
@inline function genomecoords(pos, len, record::BamRecord, recorddata::HTSReadData; onebased=true)
    if ispositive(record)
        start = recorddata.alignmap[pos]
        stop  = recorddata.alignmap[pos + len - 1]
    else
        start = recorddata.alignmap[record.core.l_qseq - (pos + len - 1) + 1]
        stop  = recorddata.alignmap[record.core.l_qseq - pos + 1]
    end
    start = ifelse(iszero(start), start, start + onebased)
    stop  = ifelse(iszero(stop), stop, stop + 1)

    start, stop
end

### accessor functions for fire

"""
    firenucpos(r::BamRecord, recorddata::StencillingData{AuxMapModFire})

Return fiber tools annotated nucleosome positions on read coordinates 0-based coordinates.
"""
@inline  firenucpos(r::BamRecord, recorddata::StencillingData{AuxMapModFire}) =  firenucpos(r, recorddata.auxmap)

"""
    firenuclen(r::BamRecord, recorddata::StencillingData{AuxMapModFire})

Return fiber tools annotated nucleosome lengths on read coordinates.
"""
@inline  firenuclen(r::BamRecord, recorddata::StencillingData{AuxMapModFire}) =  firenuclen(r, recorddata.auxmap)

"""
    firemsppos(r::BamRecord, recorddata::StencillingData{AuxMapModFire})

Return fiber tools annotated MSP positions on read coordinates. 0-based coordinates.
"""
@inline  firemsppos(r::BamRecord, recorddata::StencillingData{AuxMapModFire}) =  firemsppos(r, recorddata.auxmap)

"""
    firemsplen(r::BamRecord, recorddata::StencillingData{AuxMapModFire})

Return fiber tools annotated MSP lengths on read coordinates.
"""
@inline  firemsplen(r::BamRecord, recorddata::StencillingData{AuxMapModFire}) =  firemsplen(r, recorddata.auxmap)

"""
    firemspqual(r::BamRecord, recorddata::StencillingData{AuxMapModFire})

Return fiber tools annotated MSP quality scores on read coordinates.
"""
@inline firemspqual(r::BamRecord, recorddata::StencillingData{AuxMapModFire}) = firemspqual(r, recorddata.auxmap)


@inline  firenucpos(r::BamRecord, am::AuxMapModFire) = reinterpret(Int32, @view(r.data[am.ns.start:am.ns.stop]))
@inline  firenuclen(r::BamRecord, am::AuxMapModFire) = reinterpret(Int32, @view(r.data[am.nl.start:am.nl.stop]))
@inline  firemsppos(r::BamRecord, am::AuxMapModFire) = reinterpret(Int32, @view(r.data[am.as.start:am.as.stop]))
@inline  firemsplen(r::BamRecord, am::AuxMapModFire) = reinterpret(Int32, @view(r.data[am.al.start:am.al.stop]))
@inline firemspqual(r::BamRecord, am::AuxMapModFire) = @view r.data[am.aq.start:am.aq.stop]


"""
    firemsps(r::BamRecord, recorddata::StencillingData{AuxMapModFire})

Construct iterator of fire msp positions, lengths, and quality scores in read coordinates
Usage:

    for (pos, len, qual) in firemsps(record, recorddata)
        ### do something 
    end
"""
@inline firemsps(record::BamRecord, recorddata::StencillingData{AuxMapModFire}) = zip(firemsppos(record, recorddata.auxmap), firemsplen(record, recorddata.auxmap), firemspqual(record, recorddata.auxmap))

"""
    firenucs(r::BamRecord, recorddata::StencillingData{AuxMapModFire})

Construct iterator of fire nucleosome positions, lengths in read coordinates
Usage:

    for (pos, len) in firenucs(record, recorddata)
        ### do something 
    end
"""
@inline firenucs(record::BamRecord, recorddata::StencillingData{AuxMapModFire}) = zip(firenucpos(record, recorddata.auxmap), firenuclen(record, recorddata.auxmap))




@inline firegenomeelement(r, rdata) = x -> (genomecoords(x[1],x[2], r, rdata)..., x[3:end]...)
@inline firefilt(x, onebased=true) = (x[1] > onebased) && (x[2] > onebased)

"""
    firegenome(record::BamRecord, recorddata::StencillingData{AuxMapModFire}

Construct a function for mapping fire nucleosomes and MSP coordiates from read coordinates to genome filtering out entries that do not map

Usage:
    firegenomecoords = firegenome(record, recorddata)

    for (pos, len, qual) in firegenomcoords(firemsps(record, recorddata))
        ### do something with genome coordinates
    end
"""
@inline firegenome(record::BamRecord, recorddata::StencillingData{AuxMapModFire}) = it -> Iterators.filter(firefilt, Iterators.map(firegenomeelement(record, recorddata), it))