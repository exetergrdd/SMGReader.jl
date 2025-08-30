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

@testset "SAM output tests accessor functions and aux field iteration" begin
    sambamfile = string(bamfile, ".sam")
    samcramfile = string(cramfile, ".sam")
    writesamfile(bamfile, sambamfile)
    writesamfile(cramfile, samcramfile)
    @test filesequal(sambamfile, samfile)
    @test filesequal(samcramfile, samfile)
end




reader = open(HTSFileReader, bamfile)
record = read(reader)

samrecord(record, reader)

   sambamfile = string(bamfile, ".sam")
    samcramfile = string(cramfile, ".sam")
    writesamfile(bamfile, sambamfile)
    writesamfile(cramfile, samcramfile)
# s1 = SMGReader.seq(record)
# s2 = SMGReader.seq2(record)
# s3 = SMGReader.seq3(record)
# @assert s1 == s2 == s3
# SMGReader.seq2(record) == SMGReader.seq(record)
# sb, l = SMGReader.seq_byte(record)
# @btime SMGReader.seq(record);
# @btime SMGReader.seq2(record);
# @btime SMGReader.seq_byte(record);
# @btime SMGReader.seq3(record);

