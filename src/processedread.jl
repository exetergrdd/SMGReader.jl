##### Constants for initial buffer sizes
const MOD_6MA_BUFFER_LENGTH = 70_000
const MOD_5MC_BUFFER_LENGTH = 70_000
const READ_BUFFER_LENGTH = 200_000
const ALIGN_BLOCK_LENGTH = 5_000
const NUC_MSP_TF_LENGTH = 2_000

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
    StencillingDataKmer(auxmap::AuxMap, prob_threshold::UInt8=0x80)

Construct preallocated HTS read data building a read to genome map for chromatin stencilling setting an BAM auxillary data map `AuxMap`.
This type also tracks the modified sequence states of the read for fast k-mer extraction.
"""
struct StencillingDataKmer{T} <: HTSReadData
    alignmap::VectorBuffer{Int}
    auxmap::T
    seq_mods::VectorBuffer{UInt8}
    prob_threshold::UInt8
end
StencillingDataKmer(am::T, prob_threshold::UInt8=0x80) where {T} = StencillingDataKmer{T}(VectorBuffer{Int}(READ_BUFFER_LENGTH), am, VectorBuffer{UInt8}(READ_BUFFER_LENGTH), prob_threshold)

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
    FiberHMMData(auxmap::AuxMap)

Construct preallocated HTS read data building a read to genome map for FiberHMM setting an BAM auxillary data map `AuxMap`.
Contains preallocated buffers for FiberHMM-specific Molecular Annotations (nuc, msp, tf).
"""
struct FiberHMMData{T} <: HTSReadData
    alignmap::VectorBuffer{Int}
    auxmap::T

    nuc_starts::VectorBuffer{Int}
    nuc_lengths::VectorBuffer{Int}
    nuc_quals::VectorBuffer{UInt8}

    msp_starts::VectorBuffer{Int}
    msp_lengths::VectorBuffer{Int}

    tf_starts::VectorBuffer{Int}
    tf_lengths::VectorBuffer{Int}
    tf_tq::VectorBuffer{UInt8}
    tf_el::VectorBuffer{UInt8}
    tf_er::VectorBuffer{UInt8}
end
function FiberHMMData(am::T) where {T}
    FiberHMMData{T}(
        VectorBuffer{Int}(READ_BUFFER_LENGTH), am,
        VectorBuffer{Int}(NUC_MSP_TF_LENGTH), VectorBuffer{Int}(NUC_MSP_TF_LENGTH), VectorBuffer{UInt8}(NUC_MSP_TF_LENGTH),
        VectorBuffer{Int}(NUC_MSP_TF_LENGTH), VectorBuffer{Int}(NUC_MSP_TF_LENGTH),
        VectorBuffer{Int}(NUC_MSP_TF_LENGTH), VectorBuffer{Int}(NUC_MSP_TF_LENGTH),
        VectorBuffer{UInt8}(NUC_MSP_TF_LENGTH), VectorBuffer{UInt8}(NUC_MSP_TF_LENGTH), VectorBuffer{UInt8}(NUC_MSP_TF_LENGTH)
    )
end



function autodetecthtsdata(file::String)
    reader = open(HTSFileReader, file)
    dt = autodetecthtsdata(reader)
    close(reader)
    dt
end

"""
    autodetecthtsdata(reader::HTSFileReader)

Auto detect the HTSReadData type, with `StencillingData` or `DirectRNA`. `StencillingData` works with ONT WGS, with 5mC and 5hmC called. In the future this may be refactored
renamed to discriminate. 

Currently relies on the bam header `RG` tag `DS` entry specifying:
  - `basecall_model=dna_r10...`
  - `basecall_model=rna004...`

"""
function autodetecthtsdata(reader::HTSFileReader)
    n = @ccall libhts.sam_hdr_count_lines(reader.hdr::sam_hdr_t_p, "RG"::Cstring)::Cint

    buf = Vector{UInt8}(undef, 4096)
    ks = Ref(kstring_t(0, length(buf), pointer(buf)))

    for i in 0:(n-1)

        ret = @ccall libhts.sam_hdr_find_tag_pos(reader.hdr::sam_hdr_t_p, "RG"::Cstring, i::Cint, "DS"::Cstring, ks::Ref{kstring_t})::Cint
        ret != 0 && continue   # -1: no DS tag, -2: error

        # DS value starts here, e.g.:
        # "basecall_model=dna_r10..."
        # "basecall_model=rna004..."

        ### C potentially reallocates the ks buffer, through error if it does, this needs to be fixed
        ks[].m == length(buf) || @warn "kstring resized by htslib - memory has leaked"

        ds = unsafe_string(ks[].s, ks[].l)


        if occursin("basecall_model=dna", ds)
            return StencillingData
        elseif occursin("basecall_model=rna", ds)
            return DirectRNA
        end
    end

    error("HTS Data type could not be inferred from $(reader.file) with $n RG records.")
end

"""
    processread!(record::BamRecord, recorddata::HTSReadData)

Process read mapping auxillary fields and building read to genome map using preallocated `recorddata`.

"""
@inline function processread!(record::BamRecord, recorddata::HTSReadData)
    readtogenome!(record, recorddata.alignmap)
    auxmap!(record, recorddata.auxmap)
end

@inline function processread!(record::BamRecord, recorddata::StencillingDataKmer)
    readtogenome!(record, recorddata.alignmap)
    auxmap!(record, recorddata.auxmap)
    
    # 1. Resize the seq_mods buffer to query length
    len = querylength(record)
    setlength!(recorddata.seq_mods, len)
    
    # 2. Pass 1: Decode standard bases into the buffer
    populate_standard_bases!(record, recorddata.seq_mods)
    
    # 3. Pass 2: Overlay modifications using ModIterator
    for modinfo in ModIterator(record, recorddata)
        if modinfo.prob >= recorddata.prob_threshold
            if modinfo.mod == mod_5mC
                recorddata.seq_mods[modinfo.pos] = 5
            elseif modinfo.mod == mod_5hmC
                recorddata.seq_mods[modinfo.pos] = 6
            elseif modinfo.mod == mod_6mA
                recorddata.seq_mods[modinfo.pos] = 4
            end
        end
    end
    auxmap!(record, recorddata.auxmap)

end

@inline function processread!(record::BamRecord, recorddata::DirectRNAAlignBlocks)
    readtogenome!(record, recorddata.alignmap)
    alignmentblocks!(record, recorddata.alignblocks)
    auxmap!(record, recorddata.auxmap)
end

@inline function processread!(record::BamRecord, recorddata::FiberHMMData)
    setlength!(recorddata.nuc_starts, 0)
    setlength!(recorddata.nuc_lengths, 0)
    setlength!(recorddata.nuc_quals, 0)
    setlength!(recorddata.msp_starts, 0)
    setlength!(recorddata.msp_lengths, 0)
    setlength!(recorddata.tf_starts, 0)
    setlength!(recorddata.tf_lengths, 0)
    setlength!(recorddata.tf_tq, 0)
    setlength!(recorddata.tf_el, 0)
    setlength!(recorddata.tf_er, 0)

    readtogenome!(record, recorddata.alignmap)
    auxmap!(record, recorddata.auxmap) || return false
    parse_fiberhmm_ma!(record, recorddata)

    return true
end

## example MA format: MA:Z:1287;nuc+Q:106-149,315-131,506-193,710-110,836-54,906-382;msp+:1-105,255-60,446-60,699-11,820-86;tf+QQQ:18-59

@inline function advance_until(data, idx, length, stopat; skipppast=false)
    while idx <= length && (data[idx] != stopat)
        idx += 1
    end
    skipppast && (idx < length) && (idx += 1)
    return idx
end

@inline function parsenumber(data, index, length, termchar)
    val = 0
    while index <= length && (data[index] != termchar) && (data[index] != UInt8(';'))
        val = val * 10 + (data[index] - UInt8('0'))
        index += 1
    end

    return index, val
end

@inline function parse_fiberhmm_ma!(record::BamRecord, recorddata::FiberHMMData)
    ma_start = recorddata.auxmap.ma.start
    ma_stop = recorddata.auxmap.ma.stop - 1 # exclude null

    aq_vals = if !isnothing(recorddata.auxmap.AQ)
        @view record.data[recorddata.auxmap.AQ.start:recorddata.auxmap.AQ.stop]
    else
        view(UInt8[], 1:0)
    end
    aq_index = 1


    # skip read length todo: update to store this in read data
    index = advance_until(record.data, ma_start, ma_stop, UInt8(';'), skipppast=true)

    ### total counts for nucleosomes, msps, and tfs
    nuc_n = 0
    msp_n = 0
    tf_n = 0

    while index <= ma_stop
        # read type string
        typestart = index
        index = advance_until(record.data, index, ma_stop, UInt8(':'))
        typestop = index - 1
        index += 1 # skip ':'

        typelength = typestop - typestart + 1
        is_nuc = false
        is_msp = false
        is_tf = false
        if typelength >= 3 && record.data[typestart] == UInt8('n') && record.data[typestart+1] == UInt8('u') && record.data[typestart+2] == UInt8('c')
            is_nuc = true
        elseif typelength >= 3 && record.data[typestart] == UInt8('m') && record.data[typestart+1] == UInt8('s') && record.data[typestart+2] == UInt8('p')
            is_msp = true
        elseif typelength >= 2 && record.data[typestart] == UInt8('t') && record.data[typestart+1] == UInt8('f')
            is_tf = true
        end

        #### identify strand and number of quality values 
        ## type[strand][Q/P]*
        typestrand = UInt8('.')
        qualcount = 0
        for k in typestop:-1:typestart
            c = record.data[k]
            if c == UInt8('+') || c == UInt8('-') || c == UInt8('.')
                typestrand = c
                break
            end
            if c == UInt8('P') || c == UInt8('Q')
                qualcount += 1
            end
        end

        @assert (is_nuc && qualcount == 1) || (is_msp && qualcount == 0) || (is_tf && qualcount == 3)
        # parse coordinates: start-len,start-len
        while index <= ma_stop
            index, regstart = parsenumber(record.data, index, ma_stop, UInt8('-'))
            index += 1 # skip '-'

            index, reglen = parsenumber(record.data, index, ma_stop, UInt8(','))


            if is_nuc
                nuc_n += 1
                setlength!(recorddata.nuc_starts, nuc_n)
                setlength!(recorddata.nuc_lengths, nuc_n)
                setlength!(recorddata.nuc_quals, nuc_n)
                (aq_index > length(aq_vals)) && error("Error in MA field parsing. Expected more aq values in $(record.qname)")


                recorddata.nuc_starts[nuc_n] = regstart
                recorddata.nuc_lengths[nuc_n] = reglen
                recorddata.nuc_quals[nuc_n] = aq_vals[aq_index]
            elseif is_msp
                msp_n += 1
                setlength!(recorddata.msp_starts, msp_n)
                setlength!(recorddata.msp_lengths, msp_n)

                recorddata.msp_starts[msp_n] = regstart
                recorddata.msp_lengths[msp_n] = reglen

            elseif is_tf
                tf_n += 1
                setlength!(recorddata.tf_starts, tf_n)
                setlength!(recorddata.tf_lengths, tf_n)
                setlength!(recorddata.tf_tq, tf_n)
                setlength!(recorddata.tf_el, tf_n)
                setlength!(recorddata.tf_er, tf_n)
                (aq_index > length(aq_vals)) && error("Error in MA field parsing. Expected more aq values in $(record.qname)")

                recorddata.tf_starts[tf_n] = regstart
                recorddata.tf_lengths[tf_n] = reglen

                recorddata.tf_tq[tf_n] = aq_vals[aq_index]
                recorddata.tf_el[tf_n] = aq_vals[aq_index+1]
                recorddata.tf_er[tf_n] = aq_vals[aq_index+2]
            end

            aq_index += qualcount


            if index <= ma_stop && record.data[index] == UInt8(';')
                index += 1
                break
            elseif index <= ma_stop && record.data[index] == UInt8(',')
                index += 1
            else
                break
            end
        end
    end
end



### accessor functions
"""
    rightposition(recorddata::HTSReadData)

Get the right position of alignment on genome. 0-based exclusive coordinates
"""
@inline rightposition(recorddata::HTSReadData) = recorddata.alignmap[findlast(!iszero, recorddata.alignmap)] + 1


"""
    referencespan(record::BamRecord, recorddata::HTSReadData) 

Get the reference span on the genome of the read.
"""
@inline referencespan(record::BamRecord, recorddata::HTSReadData) = rightposition(recorddata) - leftposition(record)


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


"""
    polyatail(record::BamRecord, recorddata::DirectRNA{AuxMapModPolyA}, default=nothing)

    Return the polyA tail length and default=nothing if not present.
"""
@inline polyAtaillength(record::BamRecord, recorddata::DirectRNA{AuxMapModPolyA}, default=nothing) = isnothing(recorddata.auxmap.pt) ? default : reinterpret(Int32, (
    record.data[recorddata.auxmap.pt.start],
    record.data[recorddata.auxmap.pt.start+1],
    record.data[recorddata.auxmap.pt.start+2],
    record.data[recorddata.auxmap.pt.start+3]
))
"""
    haspolyatail(record::BamRecord, recorddata::DirectRNA{AuxMapModPolyA})

    Check if read has a polyA tail.
"""
@inline haspolyAtail(record::BamRecord, recorddata::DirectRNA{AuxMapModPolyA}) = isnothing(recorddata.auxmap.pt)



#### functions to access alignment and basecall (nanopore) qc fields 

"""
    editdistance(record::BamRecord, recorddata::StencillingData{AuxMapModFireQC})

    Return the edit distance
"""
@inline editdistance(record::BamRecord, recorddata::StencillingData{AuxMapModFireQC}) = getauxint_flexible(record, recorddata.auxmap.NM)

"""
    alignmentscore(record::BamRecord, recorddata::StencillingData{AuxMapModFireQC})

    Return the alignment score
"""
@inline alignmentscore(record::BamRecord, recorddata::StencillingData{AuxMapModFireQC}) = getauxint_flexible(record, recorddata.auxmap.AS)


"""
    gapcompresseddivergence(record::BamRecord, recorddata::StencillingData{AuxMapModFireQC})

    Return the minimap2 gapcompressed divergence
"""
@inline gapcompresseddivergence(record::BamRecord, recorddata::StencillingData{AuxMapModFireQC}) = getauxfloat(record, recorddata.auxmap.de)


"""
    basecallqscore(record::BamRecord, recorddata::StencillingData{AuxMapModFireQC}, default=nothing)

    Return the dorado basecall score
"""
@inline basecallqscore(record::BamRecord, recorddata::StencillingData{AuxMapModFireQC}, default=nothing) = isnothing(recorddata.auxmap.qs) ? default : getauxfloat(record, recorddata.auxmap.qs)




### functions to map to genome
"""
    genomecoords(pos::Int, record::BamRecord, recorddata::HTSReadData; onebased=true)

Map read coordinates to genome coordates that uses `alignmap` in `recorddata`, `onebased=true` for 1-based coordinates
"""
@inline function genomecoords(pos::Int, record::BamRecord, recorddata::HTSReadData; onebased=true)
    if ispositive(record)
        return recorddata.alignmap[pos] + onebased
    else
        return recorddata.alignmap[record.core.l_qseq-pos+1] + onebased
    end
end


"""
    genomecoords(pos::Int, len::Int, record::BamRecord, recorddata::HTSReadData; onebased=true)

Map read coordinates to genome coordates that uses `alignmap` in `recorddata`,

  - `pos` is one-based to match the alignment map
  - `one-based = true`  returns one based inclusive coords
  - `one-based = false` returns zero based exclusive coords

"""
@inline function genomecoords(pos, len, record::BamRecord, recorddata::HTSReadData; onebased=true, ftstop=false, verbose=false)
    if ispositive(record)
        starti = pos
        stopi = pos + len - 1
    else
        starti = record.core.l_qseq - (pos + len - 1) + 1
        stopi = record.core.l_qseq - pos + 1
    end

    ### find first matching start base of interval
    while iszero(recorddata.alignmap[starti]) && (starti < stopi)
        starti += 1
    end
    start = recorddata.alignmap[starti]

    ### if we found no matching bases in the interval then region is unaligned
    if (starti == stopi) && iszero(start)
        return 0, 0
    end
    start = ifelse(iszero(start), start, start + onebased)

    #### Currently fiber-tools has a bug https://github.com/fiberseq/fibertools-rs/issues/90
    #### ftstop = mirrors fiber-tools output for testing
    if !ftstop
        ### search for the matching stop coordinate
        ### work backwards, this is guaranteed to find a nonzero position >= starti
        while iszero(recorddata.alignmap[stopi]) && (stopi > starti)
            stopi -= 1
        end
        stop = recorddata.alignmap[stopi]
        stop = ifelse(iszero(stop), stop, stop + 1)
    else


        ### if current stop is zero work back until you find a mapped based
        ### this is guaranteed as if there are none it would have returned above
        ### then chose the next base for exclusive coords
        ### if current stop is not zero then walk forward to find the next mapped base for exclusive coords

        if iszero(recorddata.alignmap[stopi])
            while iszero(recorddata.alignmap[stopi]) && (stopi > starti)
                stopi -= 1
            end
            stop = recorddata.alignmap[stopi]
            stop = ifelse(iszero(stop), stop, stop + 1)
        else
            original_stopi = stopi
            stopi += 1
            verbose && println(stopi, "\t", recorddata.alignmap[stopi])
            while (stopi < length(recorddata.alignmap)) && iszero(recorddata.alignmap[stopi])
                stopi += 1
                verbose && println(stopi, "\t", recorddata.alignmap[stopi])
            end
            verbose && println("Complete:\t", stopi, "\t", recorddata.alignmap[stopi])
            stop = recorddata.alignmap[stopi]

            ### if we got all the way to the end without an alignment skip back to the first nonzero mapped
            if iszero(stop)
                stop = recorddata.alignmap[original_stopi] + 1
            end
        end
    end

    start, stop
end


### TODO: rework type hierachy for HTSReadData potentially using traits (THTT)
### currently you can call firegenomecoords() for any HTSReadData subtype 

"""
    firegenomecoords(pos::Int, len::Int, record::BamRecord, recorddata::HTSReadData; onebased=true)

Special version of `genomecoords` for fiber-tools and fire elements that are encoded 0-based.
    
  - `pos` is zero-based position on the read
  - `one-based = true`  returns one based inclusive coords
  - `one-based = false` returns zero based exclusive coords

"""
@inline firegenomecoords(pos, len, record, recorddata; onebased=true, ftstop=false, verbose=false) = genomecoords(pos + 1, len, record, recorddata, onebased=onebased, ftstop=ftstop, verbose=verbose)

### accessor functions for fire

"""
    firenucpos(r::BamRecord, recorddata::HTSReadData)

Return fiber tools annotated nucleosome positions on read coordinates 0-based coordinates.
"""
@inline firenucpos(r::BamRecord, recorddata::HTSReadData) = firenucpos(r, recorddata.auxmap)

"""
    firenuclen(r::BamRecord, recorddata::HTSReadData)

Return fiber tools annotated nucleosome lengths on read coordinates.
"""
@inline firenuclen(r::BamRecord, recorddata::HTSReadData) = firenuclen(r, recorddata.auxmap)

"""
    firemsppos(r::BamRecord, recorddata::HTSReadData)

Return fiber tools annotated MSP positions on read coordinates. 0-based coordinates.
"""
@inline firemsppos(r::BamRecord, recorddata::HTSReadData) = firemsppos(r, recorddata.auxmap)

"""
    firemsplen(r::BamRecord, recorddata::HTSReadData)

Return fiber tools annotated MSP lengths on read coordinates.
"""
@inline firemsplen(r::BamRecord, recorddata::HTSReadData) = firemsplen(r, recorddata.auxmap)

"""
    firemspqual(r::BamRecord, recorddata::HTSReadData)

Return fiber tools annotated MSP quality scores on read coordinates.
"""
@inline firemspqual(r::BamRecord, recorddata::HTSReadData) = firemspqual(r, recorddata.auxmap)


@inline firenucpos(r::BamRecord, am::AuxMap) = reinterpret(Int32, @view(r.data[am.ns.start:am.ns.stop]))
@inline firenuclen(r::BamRecord, am::AuxMap) = reinterpret(Int32, @view(r.data[am.nl.start:am.nl.stop]))
@inline firemsppos(r::BamRecord, am::AuxMap) = reinterpret(Int32, @view(r.data[am.as.start:am.as.stop]))
@inline firemsplen(r::BamRecord, am::AuxMap) = reinterpret(Int32, @view(r.data[am.al.start:am.al.stop]))
@inline firemspqual(r::BamRecord, am::AuxMap) = @view r.data[am.aq.start:am.aq.stop]


"""
    firemsps(r::BamRecord, recorddata::HTSReadData)

Construct iterator of fire msp positions, lengths, and quality scores in read coordinates
Usage:

    for (pos, len, qual) in firemsps(record, recorddata)
        ### do something 
    end
"""
@inline firemsps(record::BamRecord, recorddata::HTSReadData) = zip(firemsppos(record, recorddata.auxmap), firemsplen(record, recorddata.auxmap), firemspqual(record, recorddata.auxmap))

"""
    firenucs(r::BamRecord, recorddata::HTSReadData)

Construct iterator of fire nucleosome positions, lengths in read coordinates
Usage:

    for (pos, len) in firenucs(record, recorddata)
        ### do something 
    end
"""
@inline firenucs(record::BamRecord, recorddata::HTSReadData) = zip(firenucpos(record, recorddata.auxmap), firenuclen(record, recorddata.auxmap))

"""
    fiberhmm_nucs(recorddata::FiberHMMData)

Construct iterator of fiberhmm nucleosome positions, lengths, and posterior mean qualities in read coordinates
"""
@inline fiberhmm_nucs(recorddata::FiberHMMData) = zip(recorddata.nuc_starts, recorddata.nuc_lengths, recorddata.nuc_quals)

"""
    fiberhmm_msps(recorddata::FiberHMMData)

Construct iterator of fiberhmm MSP positions and lengths in read coordinates
"""
@inline fiberhmm_msps(recorddata::FiberHMMData) = zip(recorddata.msp_starts, recorddata.msp_lengths)

"""
    fiberhmm_tfs(recorddata::FiberHMMData)

Construct iterator of fiberhmm TF positions, lengths, and 3 quality scores in read coordinates
"""
@inline fiberhmm_tfs(recorddata::FiberHMMData) = zip(recorddata.tf_starts, recorddata.tf_lengths, recorddata.tf_quals1, recorddata.tf_quals2, recorddata.tf_quals3)





@inline firegenomeelement(r, rdata; onebased=false) = x -> (genomecoords(x[1] + !onebased, x[2], r, rdata)..., x[3:end]...)
@inline firefilt(x, onebased=true) = (x[1] > onebased) && (x[2] > onebased)

"""
    firegenome(record::BamRecord, recorddata::HTSReadData)

Construct a function for mapping fire nucleosomes and MSP coordiates from read coordinates to genome filtering out entries that do not map

Usage:
    firegenomecoords = firegenome(record, recorddata)

    for (pos, len, qual) in firegenomcoords(firemsps(record, recorddata))
        ### do something with genome coordinates
    end
"""
@inline firegenome(record::BamRecord, recorddata::HTSReadData; onebased=false) = it -> Iterators.filter(firefilt, Iterators.map(firegenomeelement(record, recorddata; onebased=onebased), it))