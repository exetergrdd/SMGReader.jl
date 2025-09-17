# Chromatin Stencilling

This library is designed for streaming through a file handling a single read at a time avoiding excessive memory allocation. While it is possible to collect multiple reads into memory, see example below, this should be discouraged currently. 

## Iteration of modification

Below is an example of calculating the total number of modification calls, note this is a simple example for all other uses you will want to preallocate `StencillingData` to help manage memory used during iteration:
```julia
using SMGReader
using  DataStructures
file = "chromsten.cram"
reader = open(HTSFileReader, file)

modcounts = counter(Modification)

for record in eachrecord(reader)
    for mi in ModIterator(record)
        push!(modcounts, mi.mod)
    end
end
close(reader)

modcounts
```
The `ModIterator` returns `ModificationInfo` which has the following fields
```julia
    mi.mod::Modification  ### the modification
    strand::Bool          ### the strand of the modification
    pos::Int              ### the position of the mod on the read in the direction of the read
    prob::UInt8           ### mod probability [0,255], divide by 255 to get the probability.
```
Supported modifications that `mi.mod` can take:
Currently supported mod codes:

| Enum value   | Modification        | BAM code |
|--------------|--------------------|----------|
| `mod_6mA`    | N6-methyladenine   | `a`      |
| `mod_5mC`    | 5-methylcytosine   | `m`    |
| `mod_5hmC`   | 5-hydroxymethylcytosine | `h` |
| `mod_inosine`| Inosine            | `17596`      |
| `mod_pseU`   | Pseudouridine      | `17802`      |
| `mod_4mC`    | N4-methylcytosine  | `21839`    |
| `mod_2OmeA`  | 2'-O-methyladenosine | `69426` |
| `mod_2OmeC`  | 2'-O-methylcytidine | `19228` |
| `mod_2OmeG`  | 2'-O-methylguanosine | `19229` |
| `mod_2OmeU`  | 2'-O-methyluridine | `19227` |

### Mapping Read Coordinates to the Genome
In BAM format the modifications are stored with respect to the strand of the read, we need to map these to the genome using the alignment information of the read. Below is an example for making a 6mA coverage track over a genomic interval. For this we need to preallocate a `StencillingData` type, this allocates memory to store information on how the read maps to the genome and index where the modification information is stored for each read:

```julia
using SMGReader
using  DataStructures
file = "chromsten.cram"
reader = open(HTSFileReader, file)

chrom, loc = "chr10", 71027756:71163638 ### coordinates of the HK1 locus ## 1-based coords
modthresh = 0.9 ## anything with mod prob > 0.9 will be counted a modified
n = length(loc)
modcounts = zeros(n)
basecounts = zeros(n)

### Preallocate memory to store information about the read
### AuxMapMod() specifies the MM/ML fields in the auxillary read data
recorddata = StencillingData(AuxMapMod()) 

for record in eachintersection(reader, chrom, loc)
    processread!(record, recorddata)                 ## update reccorddata with new read 
    for mi in ModIterator(record, recorddata) 
        (mi.mod == mod_6mA) || continue              ## require that the mod is 6mA
        genomepos = recorddata.alignmap[mi.pos]      ## 0-based coordinate
        if !iszero(genomepos)                        ## if 0 position not in genome alignment
            locpos = (genomepos + 1) - loc.start + 1 ## 1-based coords in interval
            if mi.prob > 0.9*255                     ## check its modified
                modcounts[locpos] += 1               ## total modified counts at position
            end
            basecounts[locpos] += 1                  ## total bases called at position
        end
    end
end
close(reader)

proportionmod = modcounts./basecount ### calculate the proportion modified
```

## Access to FIRE fields
There are a number of convenience functions to access FIRE fields. There are two key sets of fields:
  - Nucleosomes: auxillary field `ns` and `nl` specify their position and length in read coordinates
  - Methylation sensitive patches (MSPs): `as`, `al`, and `aq` specifies their position, length and quality
    - MSP quality <  0.9*255 intepreted as a accessible linker region
    - MSP quality >= 0.9*255 interpreted as FIRE element. 

To access this information, the key difference is a different auxillary map is needed `AuxMapModFire`.

```julia
...
recorddata = StencillingData(AuxMapModFire())
for record in eachrecord(reader)
    processrecord!(record, reccorddata)
    ## iterate over the nucleosomes in read coordinates
    for (pos, len) in firenucs(record, recorddata)
        ### do something 
    end
    ## iterate over the MSPs in read coordinates
    for (pos, len, qual) in firemsps(record, recorddata)
        ### do something 
    end
end
...
```
There are functions to access the fire positions, lengths, qualities etc. as vectors rather than iterating over them.
Currently this is an internal API:
```julia
nucpos  = SMGReader.firenucpos(record, recorddata) 
nuclen  = SMGReader.firenuclen(record, recorddata)

msppos  = SMGReader.firemsppos(record, recorddata)  
msplen  = SMGReader.firemsplen(record, recorddata)  
mspqual = SMGReader.firemspqual(record, recorddata) 
```


## Mapping FIRE fields to genome coordinates
When mapping a FIRE field with a position and length to the genome, it is possible that either or both ends may be absent from the alignment. 

This is currently down by constructing a convertor for a give record with `firegenome`. This API is a bit clunky and likely will be refined. However, it will skip over any fire field that maps incompletely to the genome. 

```julia
...
recorddata = StencillingData(AuxMapModFire())
for record in eachrecord(reader)
    processrecord!(record, reccorddata)
    firegenomecoords = firegenome(record, recorddata)

    ## iterate over the nucleosomes in genome coordinates
    for (pos, len) in firegenomcoords(firenucs(record, recorddata))
        ### do something 
    end
    ## iterate over the MSPs in genome coordinates
    for (pos, len, qual) in firegenomcoords(firemsps(record, recorddata))
        ### do something with genome coordinates
    end
end
...
```