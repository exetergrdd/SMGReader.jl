module SMGReader

export  HTSFileReader, indexfile, BamRecord,
        eachintersection, eachrecord,
        qname, refname, cigarvec, cigarstring, qual, ispositive, validflag, seq

include("bufferl.jl")
include("htslib.jl")
include("htsrecord.jl")

end
