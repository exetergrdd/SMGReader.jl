# SMGReader.jl

This package is a wrapper around `htslib` C library for BAM/CRAM file reading focussed on single molecule genomics.


```@contents
```

## Functions

```@docs
indexfile(file)
Base.open(::Type{HTSFileReader}, path; idx_path=indexfile(path), bamthreads=Threads.nthreads())
Base.close(reader::HTSFileReader)
Base.read!(reader::HTSFileReader, record::BamRecord)
Base.read(reader::HTSFileReader)
eachrecord(reader::HTSFileReader)
eachintersection(reader::HTSFileReader, chrom::AbstractString, loc::UnitRange{Int})
qname(record::BamRecord)
refname(reader::HTSFileReader, record::BamRecord)
cigarvec(record)
cigarstring(record)
qual(record)
ispositive(record)
validflag(record)
```

## Index

```@index
```