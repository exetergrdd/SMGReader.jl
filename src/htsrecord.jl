#### accessor functions for an hts record


"""
    qname(record::BamRecord)

Return query name of `record`. The read name.
"""
@inline function qname(record::BamRecord)
    return unsafe_string(pointer(record.data))  # exclude trailing null
end

"""
    refname(reader::HTSFileReader, record::BamRecord)

Return the reference name of `record` read by `reader`. `reader` is necessary to load the ref name from the header of the bam/cram file.
"""
@inline refname(reader::HTSFileReader, record::BamRecord) = refname(reader.hdr, record.core.tid)
@inline function refname(header::sam_hdr_t_p, tid)
    if tid == -1
        return "*"
    end
    chrom_pointer = @ccall libhts.sam_hdr_tid2name(header::sam_hdr_t_p, tid::Cint)::Ptr{Cchar}
    unsafe_string(chrom_pointer)
end

"""
    materefname(reader::HTSFileReader, record::BamRecord)

Return the reference name of the mate of `record` read by `reader`. `reader` is necessary to load the ref name from the header of the bam/cram file.
"""
@inline materefname(reader::HTSFileReader, record::BamRecord) = refname(reader.hdr, record.core.mtid)

"""
    cigarvec(record)

Get CIGAR description of alignment of `record` returned as lengths (Vector{Int}) and codes (Vector{Char}).

To build alignment map use `readtogenome` instead of working with cigar directly
"""
@inline function cigarvec(record, ops="MIDNSHP=XB")
    off = record.core.l_qname + 1
    n = record.core.n_cigar
    lengths = Vector{Int}(undef, n)
    codes = Vector{Char}(undef, n)
    
    cigview = @view record.data[off:(off + 4*n - 1)]
    cig = reinterpret(UInt32, cigview)

    for (k, c) in enumerate(cig)
        lengths[k] = c >> 4
        ind = (c & 0xF) + 1
        codes[k] = ind < length(ops) ? ops[(c & 0xF) + 1] : '?'
        
    end
    
    return lengths, codes
end

"""
    cigarstring(record)

Return CIGAR string of `record`.
"""
cigarstring(record) = join(string.(cigarvec(record)...))




# function cigar(record::BamRecord)
#     off = record.core.l_qname + 1
#     n = record.core.n_cigar
#     bytes = @view record.data[off:off + 4n - 1]
#     return reinterpret(UInt32, bytes)
# end

"""
    seq(record)

Return seq string of `record`.
"""
function seq(record)
    l = record.core.l_qseq
    (l == 0) && return ""
    off = record.core.l_qname + 4 * record.core.n_cigar + 1
    bam_bases = UInt8['=', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
                        'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N']
    buffer = Vector{UInt8}(undef, l)
    bytepos = 1
    pos = 1
    @inbounds while pos <= (l - 1)
        byte = record.data[off + bytepos - 1]
        base1 = byte >> 4
        base2 = byte & 0x0f
        
        buffer[pos] = bam_bases[base1 + 1]
        pos += 1
        buffer[pos] = bam_bases[base2 + 1]
        pos += 1
        bytepos += 1
    end
    @inbounds if pos == l
        byte = record.data[off + bytepos - 1]
        base1 = byte >> 4
        buffer[pos] = bam_bases[base1 + 1]
    end
    String(buffer)
end

"""
    qual(record::BamRecord)

Get quality scores as byte array
"""
@inline function qual(record::BamRecord)
    l = record.core.l_qseq
    off = record.core.l_qname + 4 * record.core.n_cigar + cld(l, 2) + 1
    return @view record.data[off:off + l - 1]
end

@inline function qualstring(record::BamRecord)
    q = qual(record)
    buff = Vector{UInt8}(undef, length(q))
    @inbounds for i in eachindex(q)
        buff[i] = q[i] + 33
    end
    return String(buff)
end


"""
    flag(record::BamRecord)

Get flag of record
"""
@inline flag(record::BamRecord) = record.core.flag

"""
    ispositive(record::BamRecord)

Determine if `record` is aligned to the positive strand
"""
@inline ispositive(record::BamRecord) = (record.core.flag & 0x10) == 0

"""
    validfrag(record::BamRecord)

Return true if flag is 0 or 16.
"""
@inline validflag(record::BamRecord) = (record.core.flag == 0x0000) || (record.core.flag == 0x0010)



"""
    leftposition(record::BamRecord)

Return the position of the record in 0-based coords
"""
@inline leftposition(record::BamRecord) = record.core.pos


"""
    mappingquality(record::BamRecord)

Return mapping quality of record
"""
@inline mappingquality(record::BamRecord) = record.core.qual

"""
    mateposition(record::BamRecord)
return the mate position in 0-based coords
"""
@inline matepos(record::BamRecord) = record.core.mpos


"""
    templatelength(record::BamRecord)

return the template length of the record
"""
@inline templatelength(record::BamRecord) = record.core.isize

"""
    querylength(record::BamRecord)

Return the query length of the record
"""
@inline querylength(record::BamRecord) = record.core.l_qseq


