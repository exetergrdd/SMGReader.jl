using SMGReader
using Test

bamfile = "test/data/test.bam"
cramfile = "test/data/test.cram"
samfile = "test/data/test.sam"

directrna_bam = "test/data/pdx1.bam"
directrna_sam = "test/data/pdx1.sam"

hg38_file = "/Users/ndlo201/projects/smf/fire/hap/0307_H1_60min/fire/0307_H1_60min.fire.cram"

# autodetecthtsdata(hg38_file)

function iteratorlength(it)
    n = 0
    for i in it
        n += 1
    end
    n
end


if isdir("test/data")
    @testset "SMGReader.jl" begin
        ### stencilling data checks
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

        @test iteratorlength(itbam) == iteratorlength(itcram) == nbam == ncram
        
        regions = [("chr1", 100253156:116796142), ("chr10", 69987747:72423136), ("chr2", 1:243199373)]
        for (chrom, loc) in regions
            @test iteratorlength(eachintersection(reader_bam, chrom, loc)) == iteratorlength(eachintersection(reader_cram, chrom, loc))
        end
        
        @test iteratorlength(eachline(samfile)) == nrecords(reader_bam) == nrecords(reader_cram)

        close(reader_bam)
        close(reader_cram)
        
        @test typeof(autodetectaux(bamfile)) == AuxMapModFire
        @test typeof(autodetectaux(cramfile)) == AuxMapModFire

        ### direct rna test
        reader = open(HTSFileReader, directrna_bam)
        close(reader)
        @test indexfile(directrna_bam) == "test/data/pdx1.bam.bai"

        reader_bam = open(HTSFileReader, directrna_bam)
        itbam  = eachrecord(reader_bam)

        @test iteratorlength(itbam) == length(collect(eachline(directrna_sam)))
        close(reader_bam)
        @test typeof(autodetectaux(directrna_bam)) == AuxMapModPolyA


    end
end

@testset "VectorBuffer Tests" begin
    include("buffertests.jl")
end

function filesequal(fileA, fileB)

    ioA = open(fileA)
    ioB = open(fileB)
    for (lA, lB) in zip(eachline(ioA), eachline(ioB))
        if lA != lB 
            println("Line mismatch:", "\n", lA, "\n", lB)
            for (k, (fA, fB)) in enumerate(zip(eachsplit(lA, '\t'), eachsplit(lB, '\t')))
                if fA != fB
                    println("Field mismatch:", "\t", k, "\t", fA, "\t", fB)
                end
            end
            return false
        end
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

        samdirectrna_bam = string(directrna_bam, ".sam")
        writesamfile(directrna_bam, samdirectrna_bam)
        @test filesequal(samdirectrna_bam, directrna_sam)
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


if isdir("test/data") 
    @testset "polyA stats" begin


        ref_polya_lengths = Dict{String, Int}()
        for line in  eachline("test/data/pdx1.sam")
            readname = ""
            polyAfound = false
            for (k, field) in enumerate(eachsplit(line, '\t'))
                (k == 1) && (readname = field)

                m = match(r"pt:i:(\d+)", field)
                if !isnothing(m)
                    ref_polya_lengths[readname] = parse(Int, m[1])
                    polyAfound = true
                    break
                end
            end
        end

        reader = open(HTSFileReader, directrna_bam)
        recorddata = DirectRNA(AuxMapModPolyA())

        for record in eachrecord(reader)
            validflag(record) || continue
            processread!(record, recorddata)
        
            @test polyAtaillength(record, recorddata) == ref_polya_lengths[qname(record)]

        end
        close(reader)
    end
end


function parse_dorado_line(line)

    # @show line
    fields = split(line, "\t")

    

    read_id = fields[2]
    sequence_length_template = parse(Int, fields[10])
    mean_qscore_template = parse(Float64, fields[11])
    alignment_genome = fields[13]
    alignment_genome_start = parse(Int, fields[14])
    alignment_genome_end = parse(Int, fields[15])
    alignment_strand_start = parse(Int, fields[16])
    alignment_strand_end = parse(Int, fields[17])
    alignment_direction = fields[18]
    alignment_length = parse(Int, fields[19])
    alignment_num_aligned = parse(Int, fields[20])
    alignment_num_correct = parse(Int, fields[21])
    alignment_num_insertions = parse(Int, fields[22])
    alignment_num_deletions = parse(Int, fields[23])
    alignment_num_substitutions = parse(Int, fields[24])
    alignment_mapq = parse(Int, fields[25])
    alignment_strand_coverage = parse(Float64, fields[26])
    alignment_identity = parse(Float64, fields[27])
    alignment_accuracy = parse(Float64, fields[28])

   

    (; read_id,
    sequence_length_template,
    mean_qscore_template,
    alignment_genome,
    alignment_genome_start,
    alignment_genome_end,
    alignment_strand_start,
    alignment_strand_end,
    alignment_direction,
    alignment_length,
    alignment_num_aligned,
    alignment_num_correct,
    alignment_num_insertions,
    alignment_num_deletions,
    alignment_num_substitutions,
    alignment_mapq,
    alignment_strand_coverage,
    alignment_identity,
    alignment_accuracy)


end

parse_dorado_summary(file="test/data/test.bam.doradosummary.tsv") = map(((k, l),) -> parse_dorado_line(l), Iterators.filter(((k, line),) -> k > 1, enumerate(eachline(file))))


if isdir("test/data")
    @testset "Alignment and QC fields " begin

        dorado_summary = parse_dorado_summary()

        bamreader = open(HTSFileReader, bamfile)
        recorddata = StencillingData(AuxMapModFireQC())
        t = 0
        for (record, summary) in zip(eachrecord(bamreader), dorado_summary)
            @test processread!(record, recorddata)


            ### invariants on sequence alignment
            # ed = editdistance(record, recorddata)
            aligned_bases = alignedbases(record, recorddata)
            reference_span = referencespan(record, recorddata)
            numins = numinsertions(record)
            alignment_length = alignmentlength(record, recorddata)
            alignment_num_correct = alignmentnumcorrect(record, recorddata)

            alignment_accuracy = alignmentaccuracy(record, recorddata)
            alignment_identity = alignmentidentity(record, recorddata) 
        


            @test qname(record) == summary.read_id
            @test refname(bamreader, record) == summary.alignment_genome
            @test querylength(record) == summary.sequence_length_template
            @test leftposition(record) == summary.alignment_genome_start
            @test rightposition(recorddata) == summary.alignment_genome_end
            
 
            @test summary.alignment_strand_start == (findfirst(!iszero, recorddata.alignmap) - 1)
            @test summary.alignment_strand_end == findlast(!iszero, recorddata.alignmap)
            @test ifelse(ispositive(record), "+", "-") == summary.alignment_direction

         

            @test summary.alignment_length == alignment_length
            @test summary.alignment_num_aligned == aligned_bases
            @test summary.alignment_num_insertions == numins
            @test summary.alignment_num_correct == alignment_num_correct


            @test abs(summary.alignment_accuracy - alignment_accuracy) < 1e-6
            @test abs(summary.alignment_identity - alignment_identity) < 1e-6
            

            ######## AS and de fields
            as = alignmentscore(record, recorddata)
            de = gapcompresseddivergence(record, recorddata)               

            @test as isa Int32
            @test de isa Float32
            @test 0 <= de <= 1


            ### nanopore qscore qs fields
            qs = basecallqscore(record, recorddata)
            if !isnothing(qs) ### pac bio won't have the same qs field
                @test qs isa Float32
                @test 0 <= qs <= 60 ### strictly dorado is qscore max is 50
            end

            t += 1
        end
        
        close(bamreader)

    end
end


@testset "Remote test" begin
    htsfiles = ["http://penrose.ex.ac.uk/data/ont/data/P2/110326_stage4_Hia5/110326_stage4_Hia5/20260311_1456_P2S-02817-A_PBK54013_70d428b2/dorado_v13/fire_hg38/results/110326_stage4_Hia5/fire/110326_stage4_Hia5.fire.cram",
               "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam"]

    for h in htsfiles
        reader = open(HTSFileReader, h)
        @test !isnothing(reader)
        @test reader.idx != C_NULL
        close(reader)
    end
end

