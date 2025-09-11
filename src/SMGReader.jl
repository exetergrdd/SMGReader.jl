module SMGReader

export  HTSFileReader, indexfile, BamRecord,
        eachintersection, eachrecord,
        qname, refname, cigarvec, cigarstring, qual, ispositive, validflag, seq, flag, leftposition, mappingquality, matepos, materefname, querylength, templatelength, haplotype, hashaplotype, rightposition,
        AuxField, AuxFieldIter, AuxMapMod, AuxMapModFire, AuxMapModFiberTools, autodetectaux, Modification, ModificationInfo, ModIterator, mod_6mA, mod_5mC, mod_5hmC, mod_4mC, mod_inosine, mod_pseU,
        StencillingData, DirectRNA, DirectRNAAlignBlocks, processread!, genomecoords, firegenomecoords,
        firemsps, firenucs, firegenome,
        samrecord, writesamfile,
        VectorBuffer, setlength!

include("buffer.jl")
include("htslib.jl")
include("htsrecord.jl")
include("alignmentmap.jl")
include("auxfields.jl")
include("processedread.jl")
include("modifications.jl")
include("sam.jl")




end
