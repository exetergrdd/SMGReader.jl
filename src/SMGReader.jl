module SMGReader

export  HTSFileReader, indexfile, BamRecord,
        eachintersection, eachrecord,
        qname, refname, cigarvec, cigarstring, qual, ispositive, validflag, seq, flag, leftposition, mappingquality, matepos, materefname, querylength, templatelength, haplotype, hashaplotype, rightpos,
        AuxField, AuxMapMod, AuxMapModFire, Modification, ModificationInfo, ModIterator,
        StencillingData, DirectRNA, DirectRNAAlignBlocks, processread!,
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
