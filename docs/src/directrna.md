# Direct RNA

This documentation is a work in progress. However, the advice on the chromatin stencilling page is relevent here.

However, the main differnce is that the processed read data should:
  - `DirectRNAA(AuxMapMod())` - for just processing modified bases
  - `DirectRNAAlignBlocks(AuxMapMod())` - for processing modified bases with spliced genomic alignments

```julia
using SMGReader
using  DataStructures
file = "chromsten.cram"
reader = open(HTSFileReader, file)

modcounts = counter(Modification)

recorddata = DirectRNAAlignBlocks(AuxMapMod())

for record in eachrecord(reader)
    processread!(record, recorddata)
    for mi in ModIterator(record, recorddata)
        push!(modcounts, mi.mod)
    end
end
close(reader)

modcounts
```


## Supported Modifications