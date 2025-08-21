module SMGReader

export  HTSFileReader, indexfile, BamRecord,
        eachintersection, eachrecord,
        qname, refname, cigarvec, cigarstring, qual, ispositive, validflag, seq, flag, leftposition, mappingquality, matepos, materefname, querylength, templatelength,
        AuxField,
        Modification, getmodbasestrand,
        samrecord, writesamfile,
        VectorBuffer, setlength!

include("buffer.jl")
include("htslib.jl")
include("htsrecord.jl")
include("auxfields.jl")
include("modifications.jl")
include("sam.jl")


end
