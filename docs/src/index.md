# SMGReader.jl

This package is a wrapper around `htslib` C library for BAM/CRAM file reading focussed on single molecule genomics.

This package has been designed for fast iteration over BAM/CRAM files containing modification calls. Primarily tested and used on those output by [dorado](https://github.com/nanoporetech/dorado) from Oxford Nanopore.

The package is built for two principle applications:
  1. Chromatin stencilling - 6mA induced at accessible base - also known as Fiber-seq. This package integrates our FIRE annotated BAM/CRAM files [FIRE](https://github.com/fiberseq/FIRE).
  2. Direct RNA-seq - detection of modification on RNA

## Installation
To install:
```julia
] # package mode
add https://github.com/exetergrdd/SMGReader.jl
```

Note installation will include `htslib` binary included in https://github.com/JuliaBinaryWrappers/htslib_jll.jl. This may not be the most recent version of hts lib.


## Basic Usage
To iterate over all records in file order:
```julia
using SMGReader

file = "chromsten.cram"
reader = open(HTSFileReader, file)
for record in eachrecord(reader)
    ### do something
end
close(reader)
```

To iterate over record in a given region:
```julia
using SMGReader

chrom = "chr1"
loc   = 1:10000 

file = "chromsten.cram"
reader = open(HTSFileReader, file)
for record in eachintersection(reader, chrom, loc)
    ### do something
end
close(reader)
```

Note that if interested in a large proportion of reads of a file it may be faster to iterate over entire file.

## Contents
```@contents
```