@inline function _get_base_at(record::BamRecord, p::Int)
    l = record.core.l_qseq
    if 1 <= p <= l
        off = record.core.l_qname + (record.core.n_cigar << 2) + 1
        bytepos = 1 + ((p - 1) >> 1)
        byte = record.data[off+bytepos-1]
        return iseven(p) ? (byte & 0x0f) : (byte >> 4)
    end
    return 0x0f # N
end

@inline function _get_read_base(record::BamRecord, pos::Int, d::Int)
    l = record.core.l_qseq
    if ispositive(record)
        p = pos + d
        return _get_base_at(record, p)
    else
        p = l - (pos + d) + 1
        b = _get_base_at(record, p)
        # Bit-reversal complements BAM's 4-bit encoded bases (e.g., A=1 <-> T=8, C=2 <-> G=4)
        return ((b >> 3) & 0x01) | ((b >> 1) & 0x02) | ((b << 1) & 0x04) | ((b << 3) & 0x08)
    end
end

@inline function _get_kmer_bval(record::BamRecord, pos::Int, d::Int)
    base_idx = _get_read_base(record, pos, d)
    base_idx == 0x01 && return 0 # A
    base_idx == 0x02 && return 1 # C
    base_idx == 0x04 && return 2 # G
    base_idx == 0x08 && return 3 # T
    return 0                     # Default
end

"""
    kmer(record::BamRecord, pos::Int, k::Int)

Return the `k`-mer string of sequence centered at `pos` (where `k` is odd).
`pos` is interpreted as the 1-based position on the original 5'->3' read sequence.
If the read is mapped to the reverse strand, the stored BAM sequence is
automatically reverse-complemented to yield the original read sequence.
Positions outside the read sequence bounds are safely padded with 'N'.
"""
function kmer(record::BamRecord, pos::Int, k::Int)
    @assert isodd(k) "k must be odd"
    h = k >> 1
    bam_bases = UInt8['=', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
        'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N']
    buffer = Vector{UInt8}(undef, k)

    for (i, d) in enumerate(-h:h)
        base_idx = _get_read_base(record, pos, d)
        buffer[i] = bam_bases[base_idx+1]
    end
    return String(buffer)
end

"""
    kmer_index(record::BamRecord, pos::Int, k::Int)

Return an integer in `1:4^(k-1)` representing the `k-1` bases around the central position `pos`.
`pos` is interpreted as the 1-based position on the original 5'->3' read sequence.
If the read is mapped to the reverse strand, the sequence context is automatically
adjusted to match the original read strand.
The central base at `pos` is excluded from the encoding since it is assumed fixed (e.g., the modification).
Bases A, C, G, T are encoded as 0, 1, 2, 3. Any ambiguous bases or out-of-bounds positions default to A (0).
"""
function kmer_index(record::BamRecord, pos::Int, k::Int)
    k == 3 && return kmer_index(record, pos, Val(3))
    k == 5 && return kmer_index(record, pos, Val(5))
    k == 7 && return kmer_index(record, pos, Val(7))

    @assert isodd(k) "k must be odd"
    h = k >> 1
    idx = 0

    for d in -h:h
        d == 0 && continue
        bval = _get_kmer_bval(record, pos, d)
        idx = idx * 4 + bval
    end

    return idx + 1
end

"""
    kmer_index(record::BamRecord, pos::Int, ::Val{3})

Unrolled fast-path for `k=3` k-mer index calculation.
"""
@inline function kmer_index(record::BamRecord, pos::Int, ::Val{3})
    bval1 = _get_kmer_bval(record, pos, -1)
    bval2 = _get_kmer_bval(record, pos, 1)
    return (bval1 * 4 + bval2) + 1
end

"""
    kmer_index(record::BamRecord, pos::Int, ::Val{5})

Unrolled fast-path for `k=5` k-mer index calculation.
"""
@inline function kmer_index(record::BamRecord, pos::Int, ::Val{5})
    bval1 = _get_kmer_bval(record, pos, -2)
    bval2 = _get_kmer_bval(record, pos, -1)
    bval3 = _get_kmer_bval(record, pos, 1)
    bval4 = _get_kmer_bval(record, pos, 2)
    return (((bval1 * 4 + bval2) * 4 + bval3) * 4 + bval4) + 1
end

"""
    kmer_index(record::BamRecord, pos::Int, ::Val{7})

Unrolled fast-path for `k=7` k-mer index calculation.
"""
@inline function kmer_index(record::BamRecord, pos::Int, ::Val{7})
    bval1 = _get_kmer_bval(record, pos, -3)
    bval2 = _get_kmer_bval(record, pos, -2)
    bval3 = _get_kmer_bval(record, pos, -1)
    bval4 = _get_kmer_bval(record, pos, 1)
    bval5 = _get_kmer_bval(record, pos, 2)
    bval6 = _get_kmer_bval(record, pos, 3)
    return (((((bval1 * 4 + bval2) * 4 + bval3) * 4 + bval4) * 4 + bval5) * 4 + bval6) + 1
end

@inline kmer3_index(record::BamRecord, pos::Int) = kmer_index(record, pos, Val(3))
@inline kmer5_index(record::BamRecord, pos::Int) = kmer_index(record, pos, Val(5))
@inline kmer7_index(record::BamRecord, pos::Int) = kmer_index(record, pos, Val(7))