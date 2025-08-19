module SMGReader

export  HTSFileReader, indexfile, BamRecord,
        eachintersection, eachrecord,
        qname, refname, cigarvec, cigarstring, qual, ispositive, validflag, seq,
        samrecord,
        VectorBuffer, setlength!

include("buffer.jl")
include("htslib.jl")
include("htsrecord.jl")
include("sam.jl")

end
