using SMGReader
using Test

bamfile = "test/data/test.bam"
cramfile = "test/data/test.cram"
samfile = "test/data/test.sam"


function iteratorlength(it)
    n = 0
    for i in it
        n += 1
    end
    n
end

if isdir("test/data")
    @testset "SMGReader.jl" begin
        reader = open(HTSFileReader, bamfile)
        close(reader)

        reader = open(HTSFileReader, cramfile)
        close(reader)

        @assert indexfile(bamfile) == "test/data/test.bam.bai"
        @assert indexfile(cramfile) == "test/data/test.cram.crai"


        reader_bam  = open(HTSFileReader, bamfile)
        reader_cram = open(HTSFileReader, bamfile)
        itbam  = eachrecord(reader_bam)
        itcram = eachrecord(reader_cram)
        nbam   = length(itbam)
        ncram  = length(itcram)

        @assert iteratorlength(itbam) == iteratorlength(itcram) == nbam == ncram
        
        regions = [("chr1", 100253156:116796142), ("chr10", 69987747:72423136), ("chr2", 1:243199373)]
        for (chrom, loc) in regions
            @assert iteratorlength(eachintersection(reader_bam, chrom, loc)) == iteratorlength(eachintersection(reader_cram, chrom, loc))
        end
        
        @test iteratorlength(eachline(samfile)) == nrecords(reader_bam) == nrecords(reader_cram)

        close(reader_bam)
        close(reader_cram)

        
        @test typeof(autodetectaux(bamfile)) == AuxMapModFire
        @test typeof(autodetectaux(cramfile)) == AuxMapModFire

    end
end

@testset "VectorBuffer Tests" begin
    include("buffertests.jl")
end

function filesequal(fileA, fileB)

    ioA = open(fileA)
    ioB = open(fileB)
    for (lA, lB) in zip(eachline(ioA), eachline(ioB))
        lA == lB || return false
    end
    bothcomplete = eof(ioA) == eof(ioB)
    close(ioA)
    close(ioB)

    bothcomplete
end

if isdir("test/data")

    @testset "SAM output tests accessor functions and aux field iteration" begin
        sambamfile = string(bamfile, ".sam")
        samcramfile = string(cramfile, ".sam")
        writesamfile(bamfile, sambamfile) ## this test auxfield iteration and hts record field access
        writesamfile(cramfile, samcramfile)
        @test filesequal(sambamfile, samfile)
        @test filesequal(samcramfile, samfile)
    end

end

if isdir("test/data")
    @testset "Aux Field Iteration" begin
        bamreader = open(HTSFileReader, bamfile)

        for record in eachrecord(bamreader)
            @test !iszero(iteratorlength(AuxFieldIter(record)))
            
            @test try SMGReader.getmmml(record)
                    true
            catch 
                false
            end
            
            for af in AuxFieldIter(record)
        
                @test 1 ≤ af.start ≤ af.stop ≤ length(record.data)
                @test af.typechar ∈ UInt8['A','c','C','s','S','i','I','f','Z','H','B']
                if af.typechar == UInt8('B')
                    @test af.elemtypechar in UInt8['c','C','s','S','i','I','f']
                else
                    @test af.elemtypechar == UInt8('@')
                end
            end
        end
        close(bamreader)
    end
end

### tests for alignmap 
if isdir("test/data")
    @testset "Process data and alignmaps" begin

        @test autodetecthtsdata(bamfile) == StencillingData
        @test autodetecthtsdata(cramfile) == StencillingData
        bamreader = open(HTSFileReader, bamfile)
        @test autodetecthtsdata(bamreader) == StencillingData
        
        recorddata = StencillingData(AuxMapMod())
        blockdata = DirectRNAAlignBlocks(AuxMapMod())
        autodata = StencillingData(autodetectaux(bamreader))
        for record in eachrecord(bamreader)
            processread!(record, recorddata)
            processread!(record, blockdata)
            processread!(record, autodata)


            @test length(recorddata.alignmap) == querylength(record) 
            ### all record in this bam are mapped
            leftposindex = findfirst(!iszero, recorddata.alignmap)
            @test !isnothing(leftposindex)
            @test recorddata.alignmap[leftposindex] == leftposition(record)

            ## check alignmap is monotonic
            pos = recorddata.alignmap[leftposindex]
            @test issorted(Iterators.filter(!iszero, recorddata.alignmap))

            ### compare alignmap and align blocks
            @test recorddata.alignmap == blockdata.alignmap == autodata.alignmap
            @test blockdata.alignblocks[1].start == leftposition(record)
            @test blockdata.alignblocks[end].stop == rightposition(recorddata)  ### both are zero based exclusive!
            ### ensure blocks are nonoverlapping
            for i in 2:length(blockdata.alignblocks)
                @test blockdata.alignblocks[i-1].stop ≤ blockdata.alignblocks[i].start
            end

        end
        close(bamreader)
    end
end

### tests for Modiciations
@testset  "Modifications" begin
    @test SMGReader.getmodbasestrand(UInt8('A'), true)  == 0x01
    @test SMGReader.getmodbasestrand(UInt8('C'), true)  == 0x02
    @test SMGReader.getmodbasestrand(UInt8('G'), true)  == 0x04
    @test SMGReader.getmodbasestrand(UInt8('T'), true)  == 0x08
    @test SMGReader.getmodbasestrand(UInt8('A'), false) == 0x08
    @test SMGReader.getmodbasestrand(UInt8('C'), false) == 0x04
    @test SMGReader.getmodbasestrand(UInt8('G'), false) == 0x02
    @test SMGReader.getmodbasestrand(UInt8('T'), false) == 0x01
    @test_throws ErrorException SMGReader.getmodbasestrand(UInt8('N'), true)
end 


@testset "Parse Mod code" begin
    s = "amhACTG175961780221839"
    data = Vector{UInt8}(s)
    i = 1

    mod, i = SMGReader.parse_mod_code(data, i)
    @test mod == mod_6mA  # 'a'
    mod, i = SMGReader.parse_mod_code(data, i)
    @test mod == mod_5mC  # 'm'
    mod, i = SMGReader.parse_mod_code(data, i)
    @test mod == mod_5hmC # 'h'
    mod, i = SMGReader.parse_mod_code(data, i)
    @test mod == SMGReader.mod_A    # 'A'
    mod, i = SMGReader.parse_mod_code(data, i)
    @test mod == SMGReader.mod_C    # 'C'
    mod, i = SMGReader.parse_mod_code(data, i)
    @test mod == SMGReader.mod_T    # 'T'
    mod, i = SMGReader.parse_mod_code(data, i)
    @test mod == SMGReader.mod_G    # 'G' 
    mod, i = SMGReader.parse_mod_code(data, i)
    @test mod == mod_inosine  # '17596'
    mod, i = SMGReader.parse_mod_code(data, i)
    @test mod == mod_pseU     # '17802'
    mod, i = SMGReader.parse_mod_code(data, i)
    @test mod == mod_4mC      # '21839'
end


if isdir("test/data")
    @testset "Modification iteration" begin
        @test autodetectmods(bamfile) == (mod_6mA, mod_5mC, mod_5hmC)
        @test autodetectmods(cramfile) == (mod_6mA, mod_5mC, mod_5hmC)
        bamreader = open(HTSFileReader, bamfile)
        recorddata = StencillingData(AuxMapMod())
        for record in eachrecord(bamreader)
            processread!(record, recorddata)
            
            modit = ModIterator(record, recorddata)
            modlength = length(modit)
            mods = collect(modit)
            @test length(mods) == modlength
            @test mods == collect(ModIterator(record)) ### confirm both ways of constructing yield same iteration
            @test all(mi -> mi.mod ∈ (mod_6mA, mod_5mC, mod_5hmC), mods)
            @test all(mi -> 1 ≤ mi.pos ≤ querylength(record), mods)

            mods_6mA  = filter(m -> m.mod == mod_6mA, mods)
            mods_5mC  = filter(m -> m.mod == mod_5mC, mods)
            mods_5hmC = filter(m -> m.mod == mod_5hmC, mods)

            @test length(mods_5mC) == length(mods_5hmC)
            @test getproperty.(mods_5mC, :pos) == getproperty.(mods_5hmC, :pos)


            ### Ensure there are not more mods than the corresponding base 
            sequence = seq(record)
            posAT = ispositive(record) ? findall(==('A'), sequence) : findall(==('T'), sequence)
            posCG = ispositive(record) ? findall(==('C'), sequence) : findall(==('G'), sequence)

            @test length(mods_6mA) ≤ length(posAT)
            @test length(mods_5mC) ≤ length(posCG) ### already confirmed that 5mC and 5hmC are identical positions
                
            @test all(m -> m.pos ∈ Set(posAT), mods_6mA)
            @test all(m -> m.pos ∈ Set(posCG), mods_5mC)
            break
        end
        
        
        close(bamreader)
    end
end

struct FireEntry
    start::Int
    stop::Int
    length::Int
    qual::Int 
    refstart::Union{Int, Nothing}
    refstop::Union{Int, Nothing}
    reflength::Union{Int, Nothing}
end

function readfiredata(file="test/data/fire_nuc_msp.test.cram.tsv")

    data = Dict{String, Dict{String, Vector{FireEntry}}}()

    for line in eachline(file)
        occursin(r"^Read", line) && continue
        fields = split(line, '\t')
        read = fields[1]
        feature = fields[2]
        start  = parse(Int, fields[3])
        stop   = parse(Int, fields[4])
        length = parse(Int, fields[5])
        qual   = parse(Int, fields[6])
        refstart  = isempty(fields[7]) ? 0 : parse(Int, fields[7])
        refstop   = isempty(fields[8]) ? 0 : parse(Int, fields[8])
        reflength = isempty(fields[9]) ? 0 : parse(Int, fields[9])

        haskey(data, read) || (data[read] = Dict{String, Vector{FireEntry}}())
        haskey(data[read], feature) || (data[read][feature] = FireEntry[])

        push!(data[read][feature], FireEntry(start, stop, length, qual, refstart, refstop, reflength))
        
    end
    data

end


if isdir("test/data")
    @testset "FIRE fields" begin
        bamreader = open(HTSFileReader, bamfile)
        recorddata = StencillingData(AuxMapModFire())
        firedata = readfiredata()
        for record in eachrecord(bamreader)
            processread!(record, recorddata)

            ### nucleoseoms
            ftnucs = firedata[qname(record)]["nucleosome"]
            fn = collect(firenucs(record, recorddata))

            if ispositive(record)
                fn_start = first.(fn)
                fn_stop = first.(fn) .+ last.(fn)
            else
                reverse!(fn)
                fn_start = querylength(record) .- (first.(fn) .+ last.(fn))
                fn_stop  = querylength(record) .- (first.(fn))
            end
            gn = firegenomecoords.(first.(fn), last.(fn), Ref(record), Ref(recorddata), onebased=false, ftstop=true)

            @test fn_start   == getproperty.(ftnucs, :start)
            @test fn_stop    == getproperty.(ftnucs, :stop)
            @test last.(fn)  == getproperty.(ftnucs, :length)
            @test first.(gn) == getproperty.(ftnucs, :refstart)
            @test last.(gn)  == getproperty.(ftnucs, :refstop)


            ### msps
            ftmsps = firedata[qname(record)]["msp"]
            fm = collect(firemsps(record, recorddata))

            if ispositive(record)
                fm_start = first.(fm)
                fm_stop  = first.(fm) .+ getindex.(fm, 2)
            else
                reverse!(fm)
                fm_start = querylength(record) .- (first.(fm) .+ getindex.(fm, 2))
                fm_stop  = querylength(record) .- (first.(fm))
            end
            gm = firegenomecoords.(first.(fm), getindex.(fm, 2), Ref(record), Ref(recorddata), onebased=false, ftstop=true)

            @test fm_start         == getproperty.(ftmsps, :start)
            @test fm_stop          == getproperty.(ftmsps, :stop)
            @test getindex.(fm, 2) == getproperty.(ftmsps, :length)
            @test last.(fm)        == getproperty.(ftmsps, :qual)
            @test first.(gm)       == getproperty.(ftmsps, :refstart)
            @test last.(gm)        == getproperty.(ftmsps, :refstop)
        end
    end
end
