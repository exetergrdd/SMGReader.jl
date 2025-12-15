### functions for decoding the MM and ML fields in a BAM record

"""
    Modification

Enumerates supported DNA/RNA modification types, such as 6mA, 5mC, 5hmC, inosine, pseudouridine, and others.
Used to identify the type of modification parsed from MM/ML fields.

Currently supported mod codes for DNA and RNA:

| Enum value   | Modification        | BAM code |
|--------------|--------------------|----------|
| `mod_6mA`    | N6-methyladenine   | `a`      |
| `mod_5mC`    | 5-methylcytosine   | `m`    |
| `mod_5hmC`   | 5-hydroxymethylcytosine | `h` |
| `mod_inosine`| Inosine            | `17596`      |
| `mod_pseU`   | Pseudouridine      | `17802`      |
| `mod_4mC`    | N4-methylcytosine  | `21839`    |
| `mod_2OmeA`  | 2'-O-methyladenosine | `69426` |
| `mod_2OmeC`  | 2'-O-methylcytidine | `19228` |
| `mod_2OmeG`  | 2'-O-methylguanosine | `19229` |
| `mod_2OmeU`  | 2'-O-methyluridine | `19227` |

"""
@enum Modification begin
    mod_6mA
    mod_5mC
    mod_5hmC
    mod_inosine
    mod_pseU
    mod_4mC
    mod_2OmeA
    mod_2OmeC
    mod_2OmeG
    mod_2OmeU
    mod_C
    mod_A
    mod_T
    mod_G
end




"""
    getmodbasestrand(b::UInt8, ispos::Bool) -> UInt8

Given a base character and strand, return the encoded base for modification search.
"""
@inline function getmodbasestrand(b::UInt8, ispos::Bool)
    
    if ispos
        (b == UInt8('A')) && return 0x01
        (b == UInt8('C')) && return 0x02
        (b == UInt8('G')) && return 0x04
        (b == UInt8('T')) && return 0x08
    else
        (b == UInt8('A')) && return 0x08
        (b == UInt8('C')) && return 0x04
        (b == UInt8('G')) && return 0x02
        (b == UInt8('T')) && return 0x01
    end

    error("unrecognised base: $b")
end

"""
    parse_mod_code(data, i)

Parse the modcode at position `i` in bytearray `data`
Return `(mod, next_i)` where `mod` is a `Modifcation` and `next_i` is the next index after the modcode in the bytearray.
"""
@inline function parse_mod_code(data, i)
    ## TODO: have performant way for these not to be hardcoded
    if data[i] == UInt8('a')
        return mod_6mA, i + 1
    elseif data[i] == UInt8('m')
        return mod_5mC, i + 1
    elseif data[i] == UInt8('h')
        return mod_5hmC, i + 1
    elseif (data[i]   == UInt8('1')) && 
           (data[i+1] == UInt8('7')) && 
           (data[i+2] == UInt8('5')) && 
           (data[i+3] == UInt8('9')) && 
           (data[i+4] == UInt8('6')) && 
        return mod_inosine, i + 5
    elseif (data[i]   == UInt8('1')) && 
           (data[i+1] == UInt8('7')) && 
           (data[i+2] == UInt8('8')) && 
           (data[i+3] == UInt8('0')) && 
           (data[i+4] == UInt8('2')) && 
        return mod_pseU, i + 5
    elseif (data[i]   == UInt8('2')) && 
           (data[i+1] == UInt8('1')) && 
           (data[i+2] == UInt8('8')) && 
           (data[i+3] == UInt8('3')) && 
           (data[i+4] == UInt8('9')) && 
        return mod_4mC, i + 5
    elseif (data[i]   == UInt8('6')) && 
           (data[i+1] == UInt8('9')) && 
           (data[i+2] == UInt8('4')) && 
           (data[i+3] == UInt8('2')) && 
           (data[i+4] == UInt8('6')) && 
        return mod_2OmeA, i + 5
    elseif (data[i]   == UInt8('1')) && 
           (data[i+1] == UInt8('9')) && 
           (data[i+2] == UInt8('2')) && 
           (data[i+3] == UInt8('2')) && 
           (data[i+4] == UInt8('8')) && 
        return mod_2OmeC, i + 5
    elseif (data[i]   == UInt8('1')) && 
           (data[i+1] == UInt8('9')) && 
           (data[i+2] == UInt8('2')) && 
           (data[i+3] == UInt8('2')) && 
           (data[i+4] == UInt8('9')) && 
        return mod_2OmeG, i + 5
    elseif (data[i]   == UInt8('1')) && 
           (data[i+1] == UInt8('9')) && 
           (data[i+2] == UInt8('2')) && 
           (data[i+3] == UInt8('2')) && 
           (data[i+4] == UInt8('7')) && 
        return mod_2OmeU, i + 5

    elseif  data[i] == UInt8('C') ### invalid mod code that appears to be output but modkit when no modifications are called due to low base quality.
        return mod_C, i + 1
    elseif  data[i] == UInt8('A') ### invalid mod code that appears to be output but modkit when no modifications are called due to low base quality.
        return mod_A, i + 1
    elseif  data[i] == UInt8('T') ### invalid mod code that appears to be output but modkit when no modifications are called due to low base quality.
        return mod_T, i + 1
    elseif  data[i] == UInt8('G') ### invalid mod code that appears to be output but modkit when no modifications are called due to low base quality.
        return mod_G, i + 1
    else
        iob = IOBuffer()
        while i < length(data)
            (data[i] == UInt8(',')) && break
            write(iob, data[i])
            i += 1
        end
        
        s = String(take!(iob))
        error("Unrecognised mod code: $s")
            
    end
    
end


"""
    findnextbase(record::BamRecord, base::UInt8, start::Int, skip::Int)

Find the next occurence of `base` = (0x01, 0x02, 0x04, 0x08)` for `(A, C, G, T)` within the bam bit packed sequence,
starting from position `start` skipping over `skip` occurences of the base.

This function is designed to locate the next base for which modification data is available for parsing the MM run length encoded aux field.
"""
@inline function findnextbase(record::BamRecord, base::UInt8, start::Int, skip::Int)
    l = record.core.l_qseq
    off = record.core.l_qname + (record.core.n_cigar << 2) + 1
    
    numbytes = ((l - 1) >> 1) + 1
    bytepos = 1 + ((start - 1) >> 1)
    seqpos = start
    remaining = skip + 1
    
    @inbounds begin
        ### Uses 1 based indexing if starting at an even numbered base it means we ignore the first base of the packed byte

        if iseven(start)
            byte = record.data[off + bytepos - 1]
            base2 = byte & 0x0f
            if base2 == base
                remaining -= 1
                iszero(remaining) && return seqpos
            end
            seqpos += 1
            bytepos += 1
        end
    
        while bytepos <= numbytes 
        
            byte = record.data[off + bytepos - 1]
            base1 = byte >> 4
            base2 = byte & 0x0f

            if base1 == base
                remaining -= 1
                iszero(remaining) && return seqpos
            end

            if (seqpos + 1 <= l) && base2 == base
                remaining -= 1
                iszero(remaining) && return seqpos + 1
            end
            
            bytepos += 1
            seqpos += 2
        end
    end
        
    return 0
end

"""
    findprevbase(record::BamRecord, base::UInt8, start::Int, skip::Int)

Find the prev occurence of `base` = (0x01, 0x02, 0x04, 0x08)` for `(A, C, G, T)` within the bam bit packed sequence,
starting from position `start` skipping over `skip` occurences of the base.

This function is designed to locate the prev base for which modification data is available for parsing the MM run length encoded aux field.
    
This is the reverse version of `findnextbase` for use in reverse strand aligned reads.
"""
@inline function findprevbase(record::BamRecord, base::UInt8, start::Int, skip::Int)
    l = record.core.l_qseq
    off = record.core.l_qname + (record.core.n_cigar << 2) + 1
    
    # numbytes = ((l - 1) >> 1) + 1
    bytepos = 1 + ((start - 1) >> 1)
    seqpos = start
    remaining = skip + 1
    
    @inbounds begin
        ### Uses 1 based indexing if starting at an even numbered base it means we ignore the first base of the packed byte

        if isodd(start)
            byte = record.data[off + bytepos - 1]
            base1 = byte >> 4
            if base1 == base
                remaining -= 1
                iszero(remaining) && return seqpos
            end
            seqpos -= 1
            bytepos -= 1
        end
        
        
        while bytepos >= 1
            byte = record.data[off + bytepos - 1]
            base1 = byte >> 4
            base2 = byte & 0x0f
            
             if base2 == base
                remaining -= 1
                iszero(remaining) && return seqpos
            end

            if base1 == base
                remaining -= 1
                iszero(remaining) && return seqpos - 1
            end

           
            seqpos -= 2
            bytepos -= 1
        end
    end
        
    return 0
end

#### principle iteration logic

"""
    ModificationInfo

Holds information about a single modification call parsed from a BAM/CRAM record.

# Fields
- `mod::Modification`: The modification type.
- `strand::Bool`: `true` for forward (+), `false` for reverse (-) strand.
- `pos::Int`: 1-based position in the read sequence. 
- `prob::UInt8`: Modification probability (from ML field).
"""
struct ModificationInfo
    mod::Modification
    strand::Bool
    pos::Int
    prob::UInt8
end

"""
    struct ModIterator

An iterator over all modification calls (as `ModificationInfo`) in a BAM/CRAM record, as described by the MM and ML auxiliary fields.

# Fields
- `record::BamRecord`: The BAM/CRAM record.
- `mmaux::AuxField`: The MM auxiliary field (modification code and run length encoding of).
- `mlaux::AuxField`: The ML auxiliary field (modification probabilities).
"""
struct ModIterator
    record::BamRecord
    mmaux::AuxField
    mlaux::AuxField    
end

"""
    ModIterator(record::BamRecord)

Construct a `ModIterator` for a given BAM/CRAM record, automatically extracting the MM and ML fields.
"""
function ModIterator(record::BamRecord)
    mmaux, mlaux = getmmml(record)
    return ModIterator(record, mmaux, mlaux)
end

"""
    ModIterator(record::BamRecord, rdata::HTSReadData)

Construct a `ModIterator` for a given BAM/CRAM record, using the mm and ml fields in the auxmap
"""
ModIterator(record::BamRecord, rdata::HTSReadData) = ModIterator(record, rdata.auxmap.mm, rdata.auxmap.ml)


"""
    ModIteratorState

    Mutable internal state of ModIterator

# Fields
  - `newrecord::Bool`: Bool marking whether next iteration will return a new modification with different modcode
  - `base::UInt8`: Current base associated with modcode
  - `strand::Bool`: String of modfication
  - `mod::Modification`: Current modification
  - `mm_pos::Int`: Position in the MM string
  - `ml_pos::Int`: Position in the ML vector
  - `seq_pos::Int`: Position in the bam sequence
"""
struct ModIteratorState
    newrecord::Bool
    base::UInt8
    strand::Bool
    mod::Modification
    mm_pos::Int
    ml_pos::Int
    seq_pos::Int
end


### function to iterate over the modificatino in order that they are in the file

"""
    Base.iterate(iter::ModIterator)

    Iterate through modifications in order they are present in BAM record returning `ModInfo` for each.
"""
@inline function Base.iterate(iter::ModIterator)
    initial = ModIteratorState(true, 0x00, true, mod_6mA, iter.mmaux.start, iter.mlaux.start, 1)
    iterate(iter, initial)
end


@inline function Base.iterate(iter::ModIterator, state)
    currentrun = 0

    newrecord = state.newrecord
    base      = state.base
    strand    = state.strand
    mod       = state.mod
    mm_pos    = state.mm_pos
    ml_pos    = state.ml_pos
    seq_pos   = state.seq_pos

    
    @inbounds while mm_pos < iter.mmaux.stop ### data[iter.mmaux.stop] = null
        if newrecord
            ## mod code of the format

            ## [ACGTUN][+-][modcode][.?]
            ## initial base flip to complement if its a reverse strand read
            ## base is the modification that we need to search for in the read sequence which goes with genome strand
            base = getmodbasestrand(iter.record.data[mm_pos], ispositive(iter.record))
            mm_pos += 1
            ## strand of mod on the read
            strand = iter.record.data[mm_pos] == UInt8('+')
            mm_pos += 1
            ## mod code
            mod, mm_pos = parse_mod_code(iter.record.data, mm_pos)
            ## advance past [.?]
            mm_pos += 1
            
            ### it is possible that no modifications are called 
            ### e.g. modkit call mods appears to output C+C?;A+a
            if iter.record.data[mm_pos] == UInt8(';')
                mm_pos += 1 ### skip over the ;
                newrecord = true
                continue
            end
            seq_pos = ifelse(ispositive(iter.record), 0, querylength(iter.record) + 1)
            currentrun = 0
            newrecord = false
            mm_pos += 1
        else
            if (iter.record.data[mm_pos] == UInt8(',')) || (iter.record.data[mm_pos] == UInt8(';'))
    
                if ispositive(iter.record)
                    seq_pos = findnextbase(iter.record, base, seq_pos + 1, currentrun)
                    modpos = seq_pos
                else
                    seq_pos = findprevbase(iter.record, base, seq_pos - 1, currentrun)
                    modpos = querylength(iter.record) - seq_pos + 1
                end

                modinfo = ModificationInfo(mod, strand, modpos, iter.record.data[ml_pos])         
                newrecord = iter.record.data[mm_pos] == UInt8(';')
                mm_pos += 1
                ml_pos += 1
                return modinfo, ModIteratorState(newrecord, base, strand, mod, mm_pos, ml_pos, seq_pos)         
            else
                currentrun = currentrun*10 + (iter.record.data[mm_pos] - 48)
            end 
                
            mm_pos += 1
        end
    end
    return nothing
end

@inline Base.eltype(::Type{ModIterator}) = ModificationInfo
@inline Base.IteratorSize(::ModIterator) = Base.HasLength()
@inline Base.length(iter::ModIterator) = iter.mlaux.stop - iter.mlaux.start + 1


"""
    autodetectmods(file::String; totalreads=10)

    Detect a tuple of Modifications in `file` from first `totalreads` reads.
"""
function autodetectmods(file::String; totalreads=10)
    reader = open(HTSFileReader, file)

    modpresent = falses(length(instances(Modification)))
    t = 0
    for record in eachrecord(reader)
        validflag(record) || continue
        
        for mi in ModIterator(record)
            modpresent[Int(mi.mod) + 1] = true
        end
        
        t += 1
        (t == totalreads) && break
    end
    close(reader)

    mods = ([Modification(k-1) for (k, present) in enumerate(modpresent) if present]...,)

    mods

end