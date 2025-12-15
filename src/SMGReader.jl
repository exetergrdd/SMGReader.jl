module SMGReader

export  HTSFileReader, indexfile, BamRecord, nrecords, referencedict,
        eachintersection, eachrecord,
        qname, refname, cigarvec, cigarstring, qual, ispositive, validflag, seq, flag, leftposition, mappingquality, matepos, materefname, querylength, templatelength, haplotype, hashaplotype, rightposition,
        AuxField, AuxFieldIter, AuxMapMod, AuxMapModFire, AuxMapModFiberTools, autodetectaux, Modification, ModificationInfo, ModIterator, autodetectmods,
        mod_6mA, mod_5mC, mod_5hmC, mod_4mC, mod_inosine, mod_pseU, mod_2OmeA, mod_2OmeC, mod_2OmeG, mod_2OmeU,
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
