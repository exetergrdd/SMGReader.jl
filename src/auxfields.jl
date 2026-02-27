"""
    AuxFieldIter(record::BamRecord)

An iterator over the auxiliary fields of a BAM record. 
Yields `AuxField` objects describing the location and type of each aux field in the record.
"""
struct AuxFieldIter
    record::BamRecord
end

### return type for aux iteration for aux field location 
"""
    AuxField

Describes the location and type of an auxiliary field in a BAM record.

Fields:
- `tag::NTuple{2, UInt8}`: Two-character tag identifying the field.
- `typechar::UInt8`: Type character for the field (e.g., 'A', 'Z', 'B', etc.).
- `elemtypechar::UInt8`: Element type character for array fields ('B').
- `start::Int`: Start index (1-based) of the field value in the record's byte array.
- `stop::Int`: Stop index (1-based) of the field value in the record's byte array.
"""
struct AuxField
    tag::NTuple{2, UInt8}
    typechar::UInt8
    elemtypechar::UInt8
    start::Int
    stop::Int
end
AuxField() = AuxField((UInt8(0), UInt8(0)), UInt8(0), UInt8(0), 0, 0)


"""
    Base.iterate(iter::AuxFieldIter)
    Base.iterate(iter::AuxFieldIter, i)

Iterate over the auxiliary fields in a BAM record. 
Returns a tuple `(AuxField, next_index)` for each field, or `nothing` when done.
"""
@inline function Base.iterate(iter::AuxFieldIter)
    core = iter.record.core
    offset = core.l_qname + 4 * core.n_cigar + cld(core.l_qseq, 2) + core.l_qseq
    return _next_aux(iter.record.data, offset + 1)
end
@inline Base.iterate(iter::AuxFieldIter, i) = _next_aux(iter.record.data, i)
@inline Base.IteratorSize(::AuxFieldIter) = Base.SizeUnknown()
@inline Base.eltype(::Type{AuxFieldIter}) = AuxField
@inline function _next_aux(bytes, i)
    @inbounds begin
        if i > length(bytes) - 3
            return nothing
        end

        tag = (bytes[i], bytes[i+1])
        typechar = bytes[i+2]

        i += 3
        start = i
        elemtypechar = UInt8('@') ### set initial elemtypechar to a type that does not exist

        len = 0
        if typechar == UInt8('A') ||
           typechar == UInt8('c') || typechar == UInt8('C') ||
           typechar == UInt8('s') || typechar == UInt8('S') ||
           typechar == UInt8('i') || typechar == UInt8('I') ||
           typechar == UInt8('f')
            len = 1 + (typechar in (UInt8('s'), UInt8('S')) ? 1 :
                       typechar in (UInt8('i'), UInt8('I'), UInt8('f')) ? 3 : 0)

        elseif typechar == UInt8('Z') || typechar == UInt8('H')
            nullpos = findnext(==(0x00), bytes, i)
            len = nullpos === nothing ? error("unterminated aux string") : (nullpos - start + 1)
            i = nullpos

        elseif typechar == UInt8('B')
            elemtypechar = bytes[i]
            n = reinterpret(UInt32, (
                    bytes[i+1],
                    bytes[i+2],
                    bytes[i+3],
                    bytes[i+4]
                ))
            i += 4  # skip elemtypechar and 3-byte length
            start = i + 1

            elsize = elemtypechar == UInt8('c') || elemtypechar == UInt8('C') ? 1 :
                     elemtypechar == UInt8('s') || elemtypechar == UInt8('S') ? 2 :
                     elemtypechar == UInt8('i') || elemtypechar == UInt8('I') ||
                     elemtypechar == UInt8('f') ? 4 :
                     error("Unknown elemtypechar $(Char(elemtypechar))")

            len = n * elsize
        else
            error("Unsupported aux typechar $(Char(typechar))")
        end

        stop = start + len - 1
        next_i = stop + 1
        return AuxField(tag, typechar, elemtypechar, start, stop), next_i
    end
end


"""
    getmmml(record::BamRecord)

Function to retrieve location of MM and ML auxillary fields as `AuxFields`
"""
function getmmml(record::BamRecord)
    mm_aux = nothing
    ml_aux = nothing
    
    for af in AuxFieldIter(record)
        if af.tag == (UInt8('M'), UInt8('M'))
            mm_aux = af
        elseif af.tag == (UInt8('M'), UInt8('L'))
            ml_aux = af
        end
        !isnothing(mm_aux) && !isnothing(ml_aux) && break
    end
    
    (isnothing(mm_aux) || isnothing(ml_aux)) && error("MM/ML not found")
    mm_aux, ml_aux
end



#### functions for specific and common auxillary maps
"""
    autodetectaux(file::String)

Autodetect auxillary map of file. Only checks first record and assumes all records have the same aux fields. 
"""
function autodetectaux(file::String; totalreads = 10)

    reader = open(HTSFileReader, file)

    hasmm = false
    hasml = false
    hasns = false
    hasnl = false
    hasas = false
    hasal = false
    hasaq = false
    t = 0
    for record in eachrecord(reader)
        for af in AuxFieldIter(record)
            (af.tag == (UInt8('M'), UInt8('M'))) && (hasmm = true)
            (af.tag == (UInt8('M'), UInt8('L'))) && (hasml = true)
            (af.tag == (UInt8('n'), UInt8('s'))) && (hasns = true)
            (af.tag == (UInt8('n'), UInt8('l'))) && (hasnl = true)
            (af.tag == (UInt8('a'), UInt8('s'))) && (hasas = true)
            (af.tag == (UInt8('a'), UInt8('l'))) && (hasal = true)
            (af.tag == (UInt8('a'), UInt8('q'))) && (hasaq = true)
        end
        t += 1
        (t == totalreads) && break
    end
    
    if hasmm && hasml
        if hasns || hasnl || hasas || hasal || hasaq
            if hasns && hasnl && hasas && hasal && hasaq
                auxmap = AuxMapModFire()
            elseif hasns && hasnl && hasas && hasal
                auxmap = AuxMapModFiberTools()
            else
                ### Currently there is a clash between Nanopore ns and fibertools ns tag
                # https://github.com/fiberseq/fibertools-rs/issues/97
                # https://software-docs.nanoporetech.com/dorado/latest/basecaller/sam_spec/?h=ns#read-tags
                # is hasns is true assume its this clash, if not warn that FIRE/fibertools fields are incompelte
                if !hasns
                    @warn "Incomplete FIRE/fibertools fields ignoring and tracking MM/ML only"
                end
                
                auxmap = AuxMapMod()
            end
        else
            auxmap = AuxMapMod()
        end
    else
        error("No Modification MM/ML tags in $file")
    end

    close(reader)

    auxmap
end

autodetectaux(reader::HTSFileReader) = autodetectaux(reader.file)


abstract type AuxMap end

"""
    AuxMapMod()

Construct an map of the aux data for HP, MM, ML fields. For standard use in chromatin stencilling (non-fire) and direct RNA.
"""
mutable struct AuxMapMod <: AuxMap
    hp::Union{AuxField, Nothing}
    mm::AuxField
    ml::AuxField
end
AuxMapMod() = AuxMapMod(AuxField(), AuxField(), AuxField())


"""
    AuxMapModPolyA()

Construct an map of the aux data for HP, MM, ML fields. For standard use in chromatin stencilling (non-fire) and direct RNA.
"""
mutable struct AuxMapModPolyA <: AuxMap
    hp::Union{AuxField, Nothing}
    mm::AuxField
    ml::AuxField
    pt::Union{AuxField, Nothing}
end
AuxMapModPolyA() = AuxMapModPolyA(AuxField(), AuxField(), AuxField(), AuxField())

"""
    AuxMapModFiberTools()

Construct and map of the aux data for a Fiber tools BAM/CRAM mapping fields (that missing the aq fire field)

  - hp: Haplotype
  - mm: Run length encoded mod bases
  - ml: Mod Probability
  - ns: Nucleosome positions
  - nl: nucleosome length
  - as: msp positions
  - al: msp lengths
"""
mutable struct AuxMapModFiberTools <: AuxMap
    hp::Union{AuxField, Nothing}
    mm::AuxField
    ml::AuxField
    ns::AuxField
    nl::AuxField
    as::AuxField
    al::AuxField
end
AuxMapModFiberTools() = AuxMapModFiberTools(AuxField(), AuxField(), AuxField(), AuxField(), AuxField(), AuxField(), AuxField())

"""
    AuxMapModFire()

Construct and map of the aux data for a FIRE BAM/CRAM mapping fields

  - hp: Haplotype
  - mm: Run length encoded mod bases
  - ml: Mod Probability
  - ns: Nucleosome positions
  - nl: nucleosome length
  - as: msp positions
  - al: msp lengths
  - aq: msp quality score
"""
mutable struct AuxMapModFire <: AuxMap
    hp::Union{AuxField, Nothing}
    mm::AuxField
    ml::AuxField
    ns::AuxField
    nl::AuxField
    as::AuxField
    al::AuxField
    aq::AuxField
end
AuxMapModFire() = AuxMapModFire(AuxField(), AuxField(), AuxField(), AuxField(), AuxField(), AuxField(), AuxField(), AuxField())


"""
    auxmap!(record::BamRecord, auxmap<:Auxmap)

Construct map of auxillary data from `record` identifying fields specified by `auxmap`.
"""
@inline function auxmap!(record::BamRecord, auxmap::AuxMapMod)
    hp_set = false
    mm_set = false
    ml_set = false
    
    for af in AuxFieldIter(record)
        if af.tag == (UInt8('H'), UInt8('P'))
            auxmap.hp = af
            hp_set = true
        elseif af.tag == (UInt8('M'), UInt8('M'))
            auxmap.mm = af
            mm_set = true
        elseif af.tag == (UInt8('M'), UInt8('L'))
            auxmap.ml = af
            ml_set = true
        end
        mm_set && ml_set && hp_set && break
    end
    !hp_set && (auxmap.hp = nothing)
    
    (!mm_set || !ml_set) && error("MM/ML not found")

    return true
    # (isnothing(mm) || isnothing(ml)) && (success = false) # error("MM/ML not found")
    # auxmap.hp = hp
    # auxmap.mm = mm
    # auxmap.ml = ml
    # success
end


@inline function auxmap!(record::BamRecord, auxmap::AuxMapModPolyA)
    hp_set = false
    mm_set = false
    ml_set = false
    pt_set = false

    for af in AuxFieldIter(record)
        if af.tag == (UInt8('H'), UInt8('P'))
            auxmap.hp = af
            hp_set = true
        elseif af.tag == (UInt8('M'), UInt8('M'))
            auxmap.mm = af
            mm_set = true
        elseif af.tag == (UInt8('M'), UInt8('L'))
            auxmap.ml = af
            ml_set = true
        elseif af.tag == (UInt8('p'), UInt8('t'))
            auxmap.pt = af
            pt_set = true
        end
        mm_set && ml_set && hp_set && pt_set && break
    end
    !hp_set && (auxmap.hp = nothing)
    !pt_set && (auxmap.pt = nothing)
    
    (!mm_set || !ml_set) && error("MM/ML not found")

    return true
end

@inline function auxmap!(record, auxmap::AuxMapModFire)
    hp_set = false
    mm_set = false
    ml_set = false
    ns_set = false
    nl_set = false
    as_set = false
    al_set = false
    aq_set = false
    
    totalfields = 8
    totalrecords = 0
    
    for af in AuxFieldIter(record)
        if af.tag == (UInt8('H'), UInt8('P'))
            auxmap.hp = af
            hp_set = true
            totalrecords += 1
        elseif af.tag == (UInt8('M'), UInt8('M'))
            auxmap.mm = af
            mm_set = true
            totalrecords += 1
        elseif af.tag == (UInt8('M'), UInt8('L'))
            auxmap.ml = af
            ml_set = true
            totalrecords += 1
        elseif af.tag == (UInt8('n'), UInt8('s'))
            auxmap.ns = af
            ns_set = true
            totalrecords += 1
        elseif af.tag == (UInt8('n'), UInt8('l'))
            auxmap.nl = af
            nl_set = true
            totalrecords += 1
        elseif af.tag == (UInt8('a'), UInt8('s'))
            auxmap.as = af
            as_set = true
            totalrecords += 1
        elseif af.tag == (UInt8('a'), UInt8('l'))
            auxmap.al = af
            al_set = true
            totalrecords += 1
        elseif af.tag == (UInt8('a'), UInt8('q'))
            auxmap.aq = af
            aq_set = true
            totalrecords += 1
        end
        
        (totalrecords == totalfields) && break
    end


    !hp_set && (auxmap.hp = nothing)
    (!mm_set || !ml_set) && error("MM/ML not found")
    
    if ns_set && nl_set && as_set && al_set && aq_set
        return true
    else
        return false
    end
    # (isnothing(ns) || isnothing(nl) || isnothing(as) || isnothing(al) || isnothing(aq)) && (success = false)  #error("Fire Aux fields  not found")
    
    
    # auxmap.hp = hp
    # auxmap.mm = mm
    # auxmap.ml = ml
    # auxmap.ns = ns
    # auxmap.nl = nl
    # auxmap.as = as
    # auxmap.al = al
    # auxmap.aq = aq
    # auxmap
    # success
end
