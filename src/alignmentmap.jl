
"""
    readtogenome!(record::BamRecord, alignmap::VectorBuffer{T}) where {T}

In place update of `alignmap` VectorBuffer with basepair level read to genome map.
Uses CIGAR information in `record`
"""
@inline function readtogenome!(record::BamRecord, alignmap::VectorBuffer{T}) where {T}
    setlength!(alignmap, record.core.l_qseq)
    alignmap .= 0
    _readtogenome!(record, alignmap)
    alignmap
end


### This is a allocating version of readto genome, commented as it should not be used
# @inline function readtogenome(record)
#     alignmap = zeros(Int, record.core.l_qseq)
#     _readtogenome!(record, alignmap)
#     alignmap
# end
# @inline function readtogenome!(record, alignmap=zeros(Int, record.core.l_qseq))
#     _readtogenome!(record, alignmap)
# end


@inline function _readtogenome!(record::BamRecord, alignmap::VectorBuffer{T}) where {T}
    # Extract CIGAR information
    offset = record.core.l_qname + 1  # Offset to CIGAR data
    n_cigar = record.core.n_cigar     # Number of CIGAR operations

    # View the CIGAR data as UInt32 (binary encoded CIGAR operations)
    # The view appears to produce no allocation overhead
    cigarview = @view record.data[offset:(offset + 4 * n_cigar - 1)]
    cigar = reinterpret(UInt32, cigarview)

    # Initialize positions
    readpos = 1
    genomepos = record.core.pos ## 0-based 

    # Iterate over CIGAR operations
    @inbounds @simd for c in cigar
        len = c >> 4                     # Extract length (upper 28 bits)
        op = c & 0xF                     # Extract operation (lower 4 bits)

        if op == 0  # M: Match or mismatch
            for j in 0:(len-1)
                alignmap[readpos + j] = genomepos + j
            end
            readpos += len
            genomepos += len
        elseif op == 1  # I: Insertion
            readpos += len  # Consumes read positions but not genome positions
        elseif op == 2  # D: Deletion
            genomepos += len  # Consumes genome positions but not read positions
        elseif op == 4  # S: Soft clipping
            readpos += len  # Soft clipping consumes read positions
        elseif op == 5  # H: Hard clipping
            # Hard clipping does not consume read or genome positions
            # continue
        else
            # Other operations (e.g., N, P) either skip or pad
            genomepos += ifelse(op == 3, len, 0)  # N: Skipped region
        end
    end
    alignmap
end



#### Read to genome map using alignment blocks

struct AlignBlock
    start::Int  # 0-based genomic start
    stop::Int   # 0-based exclusive genomic end
end

"""
    alignmentblocks!(record::BamRecord, alignblocks::VectorBuffer{AlignBlock})

In place update of alignment blocks for a `record`
"""
@inline function alignmentblocks!(record::BamRecord, alignblocks::VectorBuffer{AlignBlock})
    offset = record.core.l_qname + 1
    n_cigar = record.core.n_cigar
    cigarview = @view record.data[offset:(offset + 4 * n_cigar - 1)]
    cigar = reinterpret(UInt32, cigarview)
    
    setlength!(alignblocks, n_cigar) ### ensure alignment blocks has enough space

    genomepos = record.core.pos  # 0-based
    blockstart = -1
    blocki = 0

    @inbounds @simd for c in cigar
        len = c >> 4
        op = c & 0xF

        if op == 0  #|| op == 7 || op == 8  # M, =, X → aligned match block
            if blockstart < 0
                blockstart = genomepos
            end
            genomepos += len
        elseif op == 3 || op == 2  # N (skipped intron implied) or D (deletion) → block end
            if blockstart >= 0
                blocki += 1
                alignblocks[blocki] = AlignBlock(blockstart, genomepos)
                blockstart = -1
            end
            genomepos += len
        # elseif op == 1 || op == 4 || op == 5  # I/S/H → no genome consumption
        #     continue
        # else
        #     # skip unknown ops
        #     continue
        end
    end

    # Capture final block
    if blockstart >= 0
        blocki += 1
        alignblocks[blocki] = AlignBlock(blockstart, genomepos)
    end

    setlength!(alignblocks, blocki) ## restrict the alignment blocks to actual number present
    alignblocks
end