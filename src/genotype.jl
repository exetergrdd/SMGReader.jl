"""
    IndelInfo

Information about an insertion or deletion identified from an alignment map.

# Fields
  - `type::Symbol`: `:insertion` or `:deletion` (Note: spliced/skipped regions 'N' will also appear as `:deletion`)
  - `readpos::Int`: 1-based start of the insertion on the read, or `-1` for a deletion
  - `genomepos::Int`: 0-based genome position of the deletion, or 0-based genome position before the insertion
  - `length::Int`: Length of the insertion or deletion
"""
struct IndelInfo
    type::Symbol
    readpos::Int
    genomepos::Int
    length::Int
end

"""
    IndelIteratorState

Immutable internal state of `ReadIndels` iterator.

# Fields
  - `readpos::Int`: 1-based index of the current base in the read.
  - `prevgenomepos::Int`: 0-based genome position of the previous base.
"""
struct IndelIteratorState
    readpos::Int
    prevgenomepos::Int
end

"""
    ReadIndels{T<:HTSReadData}

An iterator over an `HTSReadData`'s `alignmap` that yields `IndelInfo` for each insertion or deletion.
"""
struct ReadIndels{T<:HTSReadData}
    recorddata::T
    len::Int32
end
ReadIndels(record::BamRecord, recorddata::HTSReadData) = ReadIndels(recorddata, querylength(record))

Base.IteratorSize(::Type{<:ReadIndels}) = Base.SizeUnknown()
Base.eltype(::Type{<:ReadIndels}) = IndelInfo

function Base.iterate(iter::ReadIndels)
    # Start at read position 1, with no previous genome position (-1)
    state = IndelIteratorState(1, -1)
    iterate(iter, state)
end

function Base.iterate(iter::ReadIndels, state::IndelIteratorState)
    readpos = state.readpos
    prevgenomepos = state.prevgenomepos
    alignmap = iter.recorddata.alignmap
    len = iter.len

    while readpos <= len
        genomepos = alignmap[readpos]
        
        if genomepos == 0
            start_readpos = readpos
            # Fast-forward through the insertion block
            while readpos <= len && alignmap[readpos] == 0
                readpos += 1
            end
            
            # If prevgenomepos == -1, these are leading soft-clips or unmapped.
            # If readpos > len, these are trailing soft-clips.
            # Yield true insertions that happen between mapped bases.
            if prevgenomepos != -1 && readpos <= len
                ins_len = readpos - start_readpos
                return IndelInfo(:insertion, start_readpos, prevgenomepos, ins_len), IndelIteratorState(readpos, prevgenomepos)
            end
            
        elseif prevgenomepos != -1 && genomepos > prevgenomepos + 1
            # Found a gap in the genome coordinates mapped by the read (deletion or skipped region)
            del_len = genomepos - prevgenomepos - 1
            del_genome_start = prevgenomepos + 1
            
            # Yield deletion. Next state starts at current `readpos` with `prevgenomepos = genomepos` 
            # so we don't double count the gap.
            return IndelInfo(:deletion, -1, del_genome_start, del_len), IndelIteratorState(readpos, genomepos)
            
        else
            # Normal match/mismatch base
            prevgenomepos = genomepos
            readpos += 1
        end
    end
    return nothing
end
