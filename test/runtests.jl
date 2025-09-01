using SMGReader
using Test
using BenchmarkTools

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

@testset "SMGReader.jl" begin
    reader = open(HTSFileReader, bamfile)
    close(reader)

    reader = open(HTSFileReader, cramfile)
    close(reader)

    @assert indexfile(bamfile) == "test/data/test.bam.bai"
    @assert indexfile(cramfile) == "test/data/test.cram.crai"


    reader_bam  = open(HTSFileReader, bamfile)
    reader_cram = open(HTSFileReader, bamfile)
    @assert iteratorlength(eachrecord(reader_bam)) == iteratorlength(eachrecord(reader_cram))
    
    regions = [("chr1", 100253156:116796142), ("chr10", 69987747:72423136), ("chr2", 1:243199373)]
    for (chrom, loc) in regions
        @assert iteratorlength(eachintersection(reader_bam, chrom, loc)) == iteratorlength(eachintersection(reader_cram, chrom, loc))
    end
    
    close(reader_bam)
    close(reader_cram)

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
        bamreader = open(HTSFileReader, bamfile)
        recorddata = StencillingData(AuxMapMod())
        blockdata = DirectRNAAlignBlocks(AuxMapMod())
        for record in eachrecord(bamreader)
            processread!(record, recorddata)
            processread!(record, blockdata)
            @test length(recorddata.alignmap) == querylength(record) 
            ### all record in this bam are mapped
            leftposindex = findfirst(!iszero, recorddata.alignmap)
            @test !isnothing(leftposindex)
            @test recorddata.alignmap[leftposindex] == leftposition(record)

            ## check alignmap is monotonic
            pos = recorddata.alignmap[leftposindex]
            @test issorted(Iterators.filter(!iszero, recorddata.alignmap))

            ### compare alignmap and align blocks
            @test recorddata.alignmap == blockdata.alignmap
            @test blockdata.alignblocks[1].start == leftposition(record)
            @test blockdata.alignblocks[end].stop - 1 == rightposition(recorddata)  ### alignblocks are zero based exclusive!
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
        bamreader = open(HTSFileReader, bamfile)
        recorddata = StencillingData(AuxMapMod())

    end
end







