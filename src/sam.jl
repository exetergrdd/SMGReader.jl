

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
    # sam_aux_string(record, io=io);
    println(io)
end
