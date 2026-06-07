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
    tag::NTuple{2,UInt8}
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
function autodetectaux(file::String; totalreads=10)

    reader = open(HTSFileReader, file)

    hasmm = false
    hasml = false
    hasns = false
    hasnl = false
    hasas = false
    hasal = false
    hasaq = false
    hasAQ = false
    haspt = false
    hasma = false

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
            (af.tag == (UInt8('A'), UInt8('Q'))) && (hasAQ = true)
            (af.tag == (UInt8('p'), UInt8('t'))) && (haspt = true)
            (af.tag == (UInt8('M'), UInt8('A'))) && (hasma = true)
        end
        t += 1
        (t == totalreads) && break
    end

    if hasmm && hasml
        if hasma && hasns && hasnl && hasas && hasal && hasaq
            auxmap = AuxMapModFireFiberHMM()
        elseif hasma
            auxmap = AuxMapFiberHMM()
        elseif haspt
            ### has a polyA tail length signal
            auxmap = AuxMapModPolyA()
        elseif hasns || hasnl || hasas || hasal || hasaq
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



### auxfield accessor functions

@inline getauxint(record, af) = reinterpret(Int32, (
    record.data[af.start],
    record.data[af.start+1],
    record.data[af.start+2],
    record.data[af.start+3]
))

@inline getauxuint(record, af) = reinterpret(UInt32, (
    record.data[af.start],
    record.data[af.start+1],
    record.data[af.start+2],
    record.data[af.start+3]
))



@inline function getauxint_flexible(record, af)
    if af.typechar == UInt8('c')
        return Int32(reinterpret(Int8, record.data[af.start]))
    elseif af.typechar == UInt8('C')
        return Int32(record.data[af.start])
    elseif af.typechar == UInt8('s')
        return Int32(reinterpret(Int16, (
            record.data[af.start],
            record.data[af.start+1])))
    elseif af.typechar == UInt8('S')
        return Int32(reinterpret(UInt16, (
            record.data[af.start],
            record.data[af.start+1])))
    elseif af.typechar == UInt8('i')
        return reinterpret(Int32, (
            record.data[af.start],
            record.data[af.start+1],
            record.data[af.start+2],
            record.data[af.start+3]))
    elseif af.typechar == UInt8('I')
        return Int32(reinterpret(UInt32, (
            record.data[af.start],
            record.data[af.start+1],
            record.data[af.start+2],
            record.data[af.start+3])))
    else
        error("Unsupported type")
    end
end

@inline function getauxuint_flexible(record, af::AuxField)

    if af.typechar == UInt8('C')
        return UInt32(record.data[af.start])
    elseif af.typechar == UInt8('S')
        return UInt32(reinterpret(UInt16, (
            record.data[af.start],
            record.data[af.start+1])))
    elseif af.typechar == UInt8('I')
        return reinterpret(UInt32, (
            record.data[af.start],
            record.data[af.start+1],
            record.data[af.start+2],
            record.data[af.start+3]
        ))
    else
        error("Auxfield typechar: $(af.typechar) is not UInt")
    end
end

@inline getauxfloat(record, af) = reinterpret(Float32, (
    record.data[af.start],
    record.data[af.start+1],
    record.data[af.start+2],
    record.data[af.start+3]
))



"""
    auxmap!(record::BamRecord, auxmap::AuxMap)

Construct a map of auxiliary data from `record`, identifying and extracting 
only the specific fields declared by the provided `auxmap`. 
Returns `true` if all required fields are found, `false` otherwise.
"""
function auxmap! end

"""
    @auxfieldset SetName [ExistingSets...] [require=(...)] [opt=(...)]
"""
macro auxfieldset(set_name, args...)
    all_fields = Expr[]
    
    for arg in args
        if arg isa Symbol
            # NEW: Allow sets to inherit from other sets!
            group_fields = Core.eval(__module__, arg)
            append!(all_fields, group_fields)
            
        elseif arg isa Expr && arg.head == :(=)
            wrapper = arg.args[1]
            tuple_expr = arg.args[2]
            
            (wrapper == :require || wrapper == :opt) || 
                error("Unknown keyword: $wrapper. Use 'require=' or 'opt='")
            
            items = (tuple_expr isa Expr && tuple_expr.head == :tuple) ? tuple_expr.args : [tuple_expr]
            
            for item in items
                local fname::Symbol
                local tag_str::String
                
                if item isa QuoteNode || item isa Symbol
                    sym = item isa QuoteNode ? item.value : item
                    fname = sym
                    tag_str = String(sym)
                elseif item isa String
                    fname = Symbol(item)
                    tag_str = item
                elseif item isa Expr && item.head == :call && item.args[1] == :(=>)
                    fname_raw = item.args[2]
                    tag_raw = item.args[3]
                    fname = fname_raw isa QuoteNode ? fname_raw.value : Symbol(fname_raw)
                    tag_str = tag_raw isa QuoteNode ? String(tag_raw.value) : String(tag_raw)
                else
                    error("Invalid field format: $item")
                end
                
                if wrapper == :opt
                    push!(all_fields, :($fname::Union{AuxField, Nothing} => $tag_str))
                else
                    push!(all_fields, :($fname::AuxField => $tag_str))
                end
            end
        else
            error("Invalid argument format in @auxfieldset")
        end
    end
    
    quoted_fields = [Meta.quot(f) for f in all_fields]
    
    return esc(quote
        const $set_name = [$(quoted_fields...)]
    end)
end


"""
    @auxmap StructName Set1 Set2 ...
"""
macro auxmap(struct_name, args...)
    all_fields = Expr[]
    
    # 1. Gather all the fields from the provided sets
    for arg in args
        if arg isa Symbol
            group_fields = Core.eval(__module__, arg)
            append!(all_fields, group_fields)
        else
            error("Arguments to @auxmap must be set names (e.g., ModFields)")
        end
    end
    
    # 2. Prepare the loop logic
    struct_fields = Expr[]
    inits = Any[]
    ifs = Tuple{Expr, Expr}[]
    opt_clears = Expr[]
    req_checks = Expr[]
    flags = Symbol[]
    
    for decl_pair in all_fields
        decl = decl_pair.args[2]
        tag = decl_pair.args[3]
        
        fname = decl.args[1]
        ftype = decl.args[2]
        
        t1, t2 = UInt8(tag[1]), UInt8(tag[2])
        flag = Symbol(fname, "_set")
        
        is_opt = (ftype isa Expr && ftype.head == :curly && ftype.args[1] == :Union && :Nothing in ftype.args)

        push!(struct_fields, decl)
        push!(inits, is_opt ? :nothing : :(AuxField()))
    
        push!(flags, flag)
        
        push!(ifs, (:(af.tag == ($t1, $t2)), quote
            auxmap.$fname = af
            $flag = true
            totalrecords += 1
        end))
        
        if is_opt
            push!(opt_clears, :(!$flag && (auxmap.$fname = nothing)))
        else
            push!(req_checks, :(!$flag && return false)) # Fast fail if required is missing
        end
    end
    
    total_fields = length(struct_fields)
    
    chain = :()
    for i in reverse(1:length(ifs))
        chain = Expr(:elseif, ifs[i][1], ifs[i][2], chain)
    end
    chain.head = :if 
    
    flag_inits = [:($f = false) for f in flags]
    pkg_mod = @__MODULE__ 
    
    # 3. Generate the struct and function
    return quote
        Core.@__doc__ mutable struct $(esc(struct_name)) <: AuxMap
            $(map(esc, struct_fields)...)
        end
        
        $(esc(struct_name))() = $(esc(struct_name))($(map(esc, inits)...))
        
        # Notice we removed map(esc, ...) and esc(chain). 
        # We let Julia's hygiene system handle all internal variables natively.
        @inline function $pkg_mod.auxmap!(record::BamRecord, auxmap::$(esc(struct_name)))
            $(flag_inits...)
            totalrecords = 0
            
            for af in AuxFieldIter(record)
                $chain
                totalrecords == $total_fields && break
            end
            
            $(opt_clears...)
            $(req_checks...)
            
            return true
        end
    end
end

@auxfieldset ModFields require=(mm => :MM, ml => :ML) opt=(hp => :HP)
@auxfieldset FireFields require=(:ns, :nl, :as, :al) opt=(:aq)
@auxfieldset QCFields   require=(:NM, :AS, :de) opt=(:qs)
@auxfieldset PolyAFields opt=(:pt)
@auxfieldset MAFields require=(:ma, :AQ)


abstract type AuxMap end



"""
    AuxMapMod()

Construct an map of the aux data for HP, MM, ML fields. For standard use in chromatin stencilling (non-fire) and direct RNA.
"""
@auxmap AuxMapMod ModFields


"""
    AuxMapModPolyA()

Construct an map of the aux data for HP, MM, ML fields with optional pt field for polyA tail 
"""
@auxmap AuxMapModPolyA ModFields PolyAFields


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
@auxmap AuxMapModFire ModFields FireFields



"""
    AuxMapModFireQC()

Construct and map of the aux data for a FIRE BAM/CRAM mapping fields

  - NM: edit distance
  - AS: alignment score
  - de: divergance (minimap2)
  - qs: ONT q score
  - hp: Haplotype
  - mm: Run length encoded mod bases
  - ml: Mod Probability
  - ns: Nucleosome positions
  - nl: nucleosome length
  - as: msp positions
  - al: msp lengths
  - aq: msp quality score
"""
@auxmap AuxMapModFireQC ModFields FireFields QCFields


"""
    AuxMapFiberHMM()

Construct and map of the aux data for FiberHMM BAM/CRAM mapping fields

  - hp: Haplotype
  - mm: Run length encoded mod bases
  - ml: Mod Probability
  - ma: Molecular annotation spec tag
  - AQ: Annotation quality score array
"""
@auxmap AuxMapFiberHMM ModFields MAFields

"""
    AuxMapModFireFiberHMM()

Construct an map of the aux data for both Fire and FiberHMM fields
"""
@auxmap AuxMapModFireFiberHMM ModFields FireFields MAFields


"""
    AuxMapStencilling()

Construct an map of the aux data for both Fire and FiberHMM fields
"""
@auxmap AuxMapStencilling ModFields FireFields QCFields
