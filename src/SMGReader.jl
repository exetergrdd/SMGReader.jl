module SMGReader

export HTSFileReader, indexfile, BamRecord, nrecords, referencedict,
    eachintersection, eachrecord,
    qname, refname, cigarvec, cigarstring, qual, ispositive, validflag, seq, flag, leftposition, mappingquality, matepos, materefname, querylength, templatelength,
    @auxfieldset, @auxmap,
    AuxField, AuxFieldIter, AuxMapMod, AuxMapModPolyA, AuxMapModFire, AuxMapModFiberTools, AuxMapModFireQC, AuxMapFiberHMM, AuxMapModFireFiberHMM, AuxMapStencilling,
    autodetectaux, Modification, ModificationInfo, ModIterator, autodetectmods, kmer, kmer_index, kmer3_index, kmer5_index, kmer7_index, index_to_kmer, index_to_kmer3, index_to_kmer5, index_to_kmer7,
    kmer_mod_index, kmer3_mod_index, kmer5_mod_index, kmer7_mod_index,
    haplotype, hashaplotype, polyAtaillength, haspolyAtail, editdistance, alignmentscore, gapcompresseddivergence, basecallqscore, rightposition,
    numinsertions, referencespan, alignedbases, alignmentnumcorrect, alignmentlength, alignmentidentity, alignmentaccuracy,
    mod_6mA, mod_5mC, mod_5hmC, mod_4mC, mod_inosine, mod_pseU, mod_2OmeA, mod_2OmeC, mod_2OmeG, mod_2OmeU,
    StencillingData, StencillingDataKmer, DirectRNA, DirectRNAAlignBlocks, FiberHMMData, processread!, genomecoords, firegenomecoords, autodetecthtsdata,
    firemsps, firenucs, firegenome, fiberhmm_nucs, fiberhmm_msps, fiberhmm_tfs,
    samrecord, writesamfile,
    VectorBuffer, setlength!

include("buffer.jl")
include("htslib.jl")
include("htsrecord.jl")
include("alignmentmap.jl")
include("auxfields.jl")
include("processedread.jl")
include("modifications.jl")
include("kmers.jl")
include("sam.jl")
include("alignmentaccessors.jl")



end
