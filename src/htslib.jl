using htslib_jll: libhts

function htslib_version()
    unsafe_string(ccall((:hts_version, libhts), Cstring, ()))
end


# === HTSlib C Types ===

### First htsFile layout for handling cram/bam to ensure we get the correct file pointer
struct BGZF end

struct kstring_t
    l::Csize_t
    m::Csize_t
    s::Ptr{Cchar}
end

# Main file layout
struct htsFile_layout
    flags::UInt32                     # bitfields: is_bin, is_write, is_be, is_cram, is_bgzf (1 bit each)
    lineno::Int64
    line::kstring_t
    fn::Ptr{Cchar}
    fn_aux::Ptr{Cchar}
    fp::Ptr{Cvoid}                    # union of BGZF*/cram*/hFILE* — interpret manually
    state::Ptr{Cvoid}
    format::NTuple{8, UInt8}          # enough to hold htsFormat (adjust if needed)
    idx::Ptr{Cvoid}
    fnidx::Ptr{Cchar}
    bam_header::Ptr{Cvoid}
    filter::Ptr{Cvoid}
end

@inline is_bgzf(hts::htsFile_layout) = (hts.flags & 0x10) != 0   # bit 4 (is_bgzf) — note: starts from LSB

### pointers
const htsFile_p = Ptr{htsFile_layout}
const sam_hdr_t_p = Ptr{Cvoid}
const bam1_t_p = Ptr{Cvoid}
const hts_idx_t_p = Ptr{Cvoid}
const hts_itr_t_p = Ptr{Cvoid}

# === bam1_core_t and bam1_t layout ===
struct bam1_core_t
    pos::Int64
    tid::Int32
    bin::UInt16
    qual::UInt8
    l_extranul::UInt8
    flag::UInt16
    l_qname::UInt16
    n_cigar::UInt32
    l_qseq::Int32
    mtid::Int32
    mpos::Int64
    isize::Int64
end

struct bam1_t_layout
    core::bam1_core_t       # 56 bytes
    id::UInt64              # 8 bytes
    data::Ptr{UInt8}        # 8 bytes
    l_data::Int32           # 4 bytes
    m_data::UInt32          # 4 bytes
    mempolicy::UInt32       # simulate bitfield
end



# ====== BamRecord ======
## Own wrapper types
## currently this uses an unsafe wrap of the data byte array of the ecrod
## potential performance increases to switch to pointers here

mutable struct BamRecord
    core::bam1_core_t
    data::Vector{UInt8}
end
function BamRecord(bam::bam1_t_layout)
    core = bam.core
    data = unsafe_wrap(Vector{UInt8}, bam.data, bam.l_data; own=false)
    return BamRecord(core, data)
end


# === Reader Type ===
mutable struct HTSFileReader
    file::String
    fp::htsFile_p
    hdr::sam_hdr_t_p
    idx::hts_idx_t_p
    record::bam1_t_p
    iter::hts_itr_t_p
    done::Bool
end


### automatically search for index of bam and cram files
"""
    indexfile(file)

Return a `.bai` or `.crai` file as appropriate for `bam` and `cram` files
"""
@inline function indexfile(file)
    if endswith(file, ".bam")
        return string(file, ".bai")
    elseif endswith(file, ".cram")
        return string(file, ".crai")
    else
        error("Unsupported file type: $file")
    end

end

# === Reader Setup ===
### convert to Base.open?

"""
    open(HTSFileReader, path; idx_path=indexfile(path), bamthreads=Threads.nthreads())

Use ccall to `hts_open` the file at `path`.
    - `idx_path` gives location of index and is mandatory
    - `bamthreads` specifies number of threads used by htslib
"""
@inline Base.open(::Type{HTSFileReader}, path; idx_path=indexfile(path), bamthreads=Threads.nthreads()) = openhts(path, idx_path=idx_path, bamthreads=bamthreads)
function openhts(path; idx_path=indexfile(path), bamthreads=Threads.nthreads())
   
    fp = @ccall libhts.hts_open(path::Cstring, "r"::Cstring)::htsFile_p
    fp == C_NULL && error("Failed to open $path")

    ### set reader threads
    
    @ccall libhts.hts_set_threads(fp::htsFile_p, bamthreads::Cint)::Cint

    hdr = @ccall libhts.sam_hdr_read(fp::htsFile_p)::sam_hdr_t_p
    hdr == C_NULL && error("Failed to read header")
    if isfile(idx_path)
        idx = @ccall libhts.sam_index_load2(fp::htsFile_p, path::Cstring, idx_path::Cstring)::hts_idx_t_p ### works
        idx == C_NULL && error("Failed to load index for $idx_path")
    else
        idx = C_NULL
    end

    record = @ccall libhts.bam_init1()::bam1_t_p

    HTSFileReader(path, fp, hdr, idx, record, C_NULL, false)
end

"""
    close(reader::HTSFileReader)
    
Close open `reader`.
"""
function Base.close(reader::HTSFileReader)
    reader.record != C_NULL && @ccall libhts.bam_destroy1(reader.record::bam1_t_p)::Cvoid
    reader.hdr    != C_NULL && @ccall libhts.sam_hdr_destroy(reader.hdr::sam_hdr_t_p)::Cvoid
    reader.fp     != C_NULL && @ccall libhts.hts_close(reader.fp::htsFile_p)::Cint
    reader.idx    != C_NULL && @ccall libhts.hts_idx_destroy(reader.idx::hts_idx_t_p)::Cvoid
    reader.iter   != C_NULL && @ccall libhts.hts_itr_destroy(reader.iter::hts_itr_t_p)::Cvoid
end


"""
    nrecords(reader::HTSFileReader)

Return the total number of reads in the file, using the index.
"""
function nrecords(reader::HTSFileReader)
    nref = @ccall libhts.sam_hdr_nref(reader.hdr::sam_hdr_t_p)::Cint


    total = 0
    mapped = Ref{UInt64}()
    unmapped = Ref{UInt64}()
    for tid in 0:(nref-1)
        @ccall libhts.hts_idx_get_stat(
            reader.idx::hts_idx_t_p,
            tid::Cint,
            mapped::Ref{UInt64},
            unmapped::Ref{UInt64}
        )::Cint
        total += mapped[] + unmapped[]
    end

    no_coor = @ccall libhts.hts_idx_get_n_no_coor(reader.idx::hts_idx_t_p)::UInt64
    total += no_coor

    return Int(total)
end

"""
    referencedict(reader::HTSFileReader)

Build a Int32 => String reference dict from a file

"""
function referencedict(reader::HTSFileReader)
    nrefs = @ccall libhts.sam_hdr_nref(reader.hdr::sam_hdr_t_p)::Cint
    refs = Dict{Int32,String}()
    sizehint!(refs, nrefs)
    for tid in 0:(nrefs-1)
        chrom_ptr = @ccall libhts.sam_hdr_tid2name(reader.hdr::sam_hdr_t_p, tid::Cint)::Ptr{Cchar}
        refs[tid] = unsafe_string(chrom_ptr)
    end
    return refs
end





### allocate initial bam record for use when iterating a reader
@inline function allocateinitialbam(reader::HTSFileReader)
    bam = unsafe_load(Ptr{bam1_t_layout}(reader.record))
    data = unsafe_wrap(Vector{UInt8}, bam.data, bam.l_data; own=false)
    return BamRecord(bam.core, data)
end

"""
    Base.read!(reader::HTSFileReader, record::BamRecord)

Read next record in file order from `reader` and in place store result in `record`.
"""
@inline function Base.read!(reader::HTSFileReader, record::BamRecord)
    reader.done && return false

    ret = @ccall libhts.sam_read1(reader.fp::htsFile_p, reader.hdr::sam_hdr_t_p, reader.record::bam1_t_p)::Cint
    if ret < 0
        reader.done = true
        return false
    end
    # Wrap pointer in BamRecord
    bam = unsafe_load(Ptr{bam1_t_layout}(reader.record))
    record.core = bam.core
    record.data = unsafe_wrap(Vector{UInt8}, bam.data, bam.l_data; own=false)

    return true
end

"""
    Base.read(reader::HTSFileReader)

Allocating version of `read!`
Read next record in file order from `reader` return a `record`, recommended to use `read!`.
"""
@inline function Base.read(reader::HTSFileReader)
    record = allocateinitialbam(reader)
    success = read!(reader, record)
    !success && return nothing
    record
end




########### Iteration
### itertation in file order
struct RecordIterator
    reader::HTSFileReader
    record::BamRecord
end


"""
    eachrecord(reader::HTSFileReader)

Construct an iterator of `reader` in record order. Multiple successive calls will continue iterating from current point in reader, resetting iterator requires closing and opening file.
"""
@inline function eachrecord(reader::HTSFileReader)
    record = allocateinitialbam(reader)
    return RecordIterator(reader, record)
end

@inline Base.eltype(::Type{RecordIterator}) = BamRecord
@inline Base.IteratorSize(::RecordIterator) = Base.HasLength()
@inline Base.length(iter::RecordIterator) = nrecords(iter.reader)
@inline Base.iterate(iter::RecordIterator) = iterate(iter, nothing)
@inline function Base.iterate(iter::RecordIterator, state)
    success = read!(iter.reader, iter.record)
    success || return nothing
    return (iter.record, nothing)
end


########## Region based iteration
struct RegionIterator
    reader::HTSFileReader
    chrom::String
    loc::UnitRange{Int}
    record::BamRecord
end


@inline function regstring(chrom, loc; onebased=true)
    iob = IOBuffer()
    write(iob, chrom)
    write(iob, ':')
    start = ifelse(onebased, first(loc), first(loc) + 1)
    write(iob, string(start))
    write(iob, '-')
    write(iob, string(last(loc)))
    String(take!(iob))
end

"""
    eachintersection(reader::HTSFileReader, chrom::AbstractString, loc::UnitRange{Int} ; onebased=true)

Construct an iterator of `reader` over `chrom` and `loc`.
Keyword `onebased = true` assume the interval provided is one-based

"""
function eachintersection(reader::HTSFileReader, chrom::AbstractString, loc::UnitRange{Int}; onebased=true)
    record = allocateinitialbam(reader)
    region = regstring(chrom, loc, onebased=onebased)
    set_region!(reader, region)
    return RegionIterator(reader, chrom, loc, record)
end


function set_region!(reader::HTSFileReader, region::AbstractString)
    # @assert reader.idx !== nothing "Index must be loaded to set region"
    reader.iter != C_NULL &&  @ccall libhts.hts_itr_destroy(reader.iter::hts_itr_t_p)::Cvoid
  
    reader.iter = @ccall libhts.sam_itr_querys(reader.idx::hts_idx_t_p, reader.hdr::sam_hdr_t_p, region::Cstring)::hts_itr_t_p
    reader.iter == C_NULL && error("Failed to create iterator for region: $region")
    reader.done = false
end

@inline Base.eltype(::Type{RegionIterator}) = BamRecord
@inline Base.IteratorSize(::RegionIterator) = Base.SizeUnknown()
@inline Base.iterate(iter::RegionIterator) = iterate(iter, nothing)
@inline function Base.iterate(iter::RegionIterator, state)
    reader = iter.reader
    if reader.done
        return nothing
    end

    #### hts_itr_next must be called differently for SAM and CRAM files
    #### this requires accessing load of the htsFile type
    hts = unsafe_load(reader.fp) ### TODO: this should be unsafe loaded at the start as part of the reader.
    bgzf_ptr = is_bgzf(hts) ? hts.fp : C_NULL
    ret = @ccall libhts.hts_itr_next(
            bgzf_ptr::Ptr{BGZF},              # reinterpret union safely
            reader.iter::hts_itr_t_p,
            reader.record::Ptr{bam1_t_p},
            reader.fp::htsFile_p
        )::Cint

    if ret < 0
        reader.done = true
        return nothing
    end
    
    # Wrap pointer in BamRecord
    bam = unsafe_load(Ptr{bam1_t_layout}(reader.record))
    iter.record.core = bam.core
    iter.record.data = unsafe_wrap(Vector{UInt8}, bam.data, bam.l_data; own=false)

    return (iter.record, nothing)
end



#### MultiRegionIterator

const hts_pos_t = Int64

struct hts_pair_pos_t
    start::hts_pos_t
    stop::hts_pos_t
end
 

struct hts_reglist_t
    reg::Ptr{Cchar}
    intervals::Ptr{hts_pair_pos_t}
    tid::Cint
    count::UInt32
    min_start::hts_pos_t
    max_stop::hts_pos_t
end


function constructregions(reader::HTSFileReader, chroms, locs; onebased=true)
    tid_to_chrom = referencedict(reader)
    chrom_to_tid = Dict(c => tid for (tid, c) in tid_to_chrom)
    regions = Dict{Int32, Vector{hts_pair_pos_t}}()

    for (c, l) in zip(chroms, locs)
        tid = chrom_to_tid[c]
        haskey(regions, tid) || (regions[tid] = hts_pair_pos_t[])
        push!(regions[tid], hts_pair_pos_t(ifelse(onebased, first(l) - 1, first(l)), last(l)))
    end

    n_regions = length(regions)
    c_region_list = Libc.malloc(sizeof(hts_reglist_t) * n_regions)
    hts_region_list_ptr = Ptr{hts_reglist_t}(c_region_list)

    sorted_tids = sort(collect(keys(regions)))

    for (k, tid) in enumerate(sorted_tids)
        intervals = regions[tid]
        sort!(intervals, by = x -> x.start)
        n = length(intervals)

        min_start = first(intervals).start
        max_stop  = maximum(hp -> hp.stop, intervals) ### the last interval isn't guaranteed to have the last stop position

        c_intervals = Libc.malloc(sizeof(hts_pair_pos_t) * n)
        c_intervals_ptr = Ptr{hts_pair_pos_t}(c_intervals)
        unsafe_copyto!(c_intervals_ptr, pointer(intervals), n)

        unsafe_store!(hts_region_list_ptr, hts_reglist_t(
            C_NULL,
            c_intervals_ptr,
            tid,
            UInt32(n),
            min_start,
            max_stop
        ), k)
    end

    hts_region_list_ptr, n_regions
end

struct MultiRegionIterator
    reader::HTSFileReader
    record::BamRecord
end

"""
    eachintersection(reader::HTSFileReader, chroms::Vector{<:AbstractString}, locs::Vector{UnitRange{Int}}; onebased=true)

Construct an iterator of `reader` over multiple regions `chroms` and `locs`. Keyword `onebased=true` indicates whether locs are onebased inclusive.
"""
@inline function eachintersection(reader::HTSFileReader, chroms::Vector{<:AbstractString}, locs::Vector{UnitRange{Int}}; onebased=true)
    record = allocateinitialbam(reader)
    hts_region_list_ptr, n_regions = constructregions(reader, chroms, locs, onebased=onebased)
    set_regions!(reader,  hts_region_list_ptr, n_regions)
    return MultiRegionIterator(reader, record)
end


@inline function set_regions!(reader::HTSFileReader, regions_ptr::Ptr{hts_reglist_t}, n_regions::Int)
    @assert reader.idx !== nothing "Index must be loaded to set region"
    reader.iter != C_NULL &&  @ccall libhts.hts_itr_destroy(reader.iter::hts_itr_t_p)::Cvoid
   
    reader.iter = @ccall libhts.sam_itr_regions(reader.idx::hts_idx_t_p, reader.hdr::sam_hdr_t_p, regions_ptr::Ptr{hts_reglist_t}, UInt32(n_regions)::Cuint)::hts_itr_t_p
    reader.iter == C_NULL && error("Failed to create iterator for region: $region")
    reader.done = false
end

@inline Base.eltype(::Type{MultiRegionIterator}) = BamRecord
@inline Base.IteratorSize(::MultiRegionIterator) = Base.SizeUnknown()
@inline Base.iterate(iter::MultiRegionIterator) = iterate(iter, nothing)
@inline function Base.iterate(iter::MultiRegionIterator, state)
    reader = iter.reader
    if reader.done
        return nothing
    end

    ret = @ccall libhts.hts_itr_multi_next(reader.fp::htsFile_p, reader.iter::hts_itr_t_p, reader.record::Ptr{bam1_t_p})::Cint
    if ret < 0
        reader.done = true
        return nothing
    end
    
    # Wrap pointer in BamRecord
    bam = unsafe_load(Ptr{bam1_t_layout}(reader.record))

    iter.record.core = bam.core
    iter.record.data = unsafe_wrap(Vector{UInt8}, bam.data, bam.l_data; own=false)

    return (iter.record, nothing)
end

function set_cache_size!(reader::HTSFileReader, size_mb::Int)
    # Set reference cache size (in bytes)
    ret = @ccall libhts.hts_set_opt(
        reader.fp::htsFile_p,
        0::Cint,  # HTS_OPT_CACHE_SIZE
        (size_mb * 1024 * 1024)::Cint
    )::Cint
    ret != 0 && error("Failed to set cache size")
    # @ccall libhts.hts_set_opt(reader.fp::htsFile_p, 1::Cint, 1::Cint)::Cint
    # ret != 0 && error("Failed to set cache the index")
end
