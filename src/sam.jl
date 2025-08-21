
using Printf

function writesamfile(infile, samfile, numrecords=-1)
    reader = open(HTSFileReader, infile)

    nr = 0
    io = open(samfile, "w")
    for record in eachrecord(reader)
        samrecord(record, reader, io=io)
        nr += 1
        (nr == numrecords) && break
    end
    close(io)
    close(reader)
    nr
end

function samrecord(record, reader; io=stdout)
    print(io, qname(record), '\t')
    print(io, flag(record), '\t')
    print(io, refname(reader, record), '\t')
    print(io, leftposition(record) + 1, '\t')
    print(io, mappingquality(record), '\t')
    print(io, cigarstring(record), '\t')
    print(io, materefname(reader, record), '\t')
    print(io, matepos(record) + 1, '\t')
    print(io, templatelength(record), '\t')
    print(io, seq(record), '\t')
    print(io, qualstring(record))
    samauxstring(record, io=io);
    println(io)
end



function samauxstring(record; io=stdout)
    
    for af in AuxFieldIter(record)

        write(io, '\t')
        write(io, Char(af.tag[1]))
        write(io, Char(af.tag[2]))
        write(io, ':')

        ### in sam format type chars C, c, s, S, I are all i
        samtypechar = af.typechar ∈ (UInt8('C'), UInt8('c'), UInt8('s'), UInt8('S'), UInt8('I')) ? UInt8('i') : af.typechar
        write(io, samtypechar)
        write(io, ':')
        

        if af.typechar == UInt8('A')
            write(io, Char(record.data[af.start]))
        elseif af.typechar == UInt8('c')
            print(io, Int8(record.data[af.start]))
        elseif af.typechar == UInt8('C')
            print(io, record.data[af.start])
        elseif af.typechar == UInt8('s')
            val = reinterpret(Int16, (record.data[af.start], record.data[af.start+1]))
            print(io, val)
        elseif af.typechar == UInt8('S')
            val = reinterpret(UInt16, (record.data[af.start], record.data[af.start+1]))
            print(io, val)
        elseif af.typechar == UInt8('i')
            val = reinterpret(Int32, (
                record.data[af.start],
                record.data[af.start+1],
                record.data[af.start+2],
                record.data[af.start+3]
            ))
            print(io, val)
        elseif af.typechar == UInt8('I')
            val = reinterpret(UInt32, (
                record.data[af.start],
                record.data[af.start+1],
                record.data[af.start+2],
                record.data[af.start+3]
            ))
            print(io, val)
        elseif af.typechar == UInt8('f')
            val = reinterpret(Float32, (
                record.data[af.start],
                record.data[af.start+1],
                record.data[af.start+2],
                record.data[af.start+3]
            ))
            @printf(io, "%g", val)
        elseif af.typechar in (UInt8('Z'), UInt8('H'))
            write(io, view(record.data, af.start:af.stop-1))  # Drop null
        elseif af.typechar == UInt8('B')
            T = af.elemtypechar == UInt8('c') ? Int8 :
                af.elemtypechar == UInt8('C') ? UInt8 :
                af.elemtypechar == UInt8('s') ? Int16 :
                af.elemtypechar == UInt8('S') ? UInt16 :
                af.elemtypechar == UInt8('i') ? Int32 :
                af.elemtypechar == UInt8('I') ? UInt32 :
                af.elemtypechar == UInt8('f') ? Float32 :
                error("Invalid B subtype $subtype")

            write(io, af.elemtypechar)
            write(io, ',')

            # sizeof_T = sizeof(T)
            # n = (stop - start + 1) ÷ sizeof_T
            vals = reinterpret(T, view(record.data, af.start:af.stop))

            for j in eachindex(vals)
                print(io, vals[j])
                j != length(vals) && write(io, ',')
            end
        else
            error("Unknown aux field type $typechar")
        end
       
    end
    io
end
