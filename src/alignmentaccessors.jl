
###### alignment accession 
"""
    numinsertions(record::BamRecord)

    Count number of inserted bases using CIGAR information in `record`
"""
@inline function numinsertions(record::BamRecord)
    offset = record.core.l_qname + 1
    n_cigar = record.core.n_cigar

    cigarview = @view record.data[offset:(offset + 4 * n_cigar - 1)]
    cigar = reinterpret(UInt32, cigarview)

    num_ins = 0
    @inbounds @simd for c in cigar
        op = c & 0xF
        len = c >> 4
        if op == 1  # I: insertion
            num_ins += len
        end
    end
    return num_ins
end

"""
    alignedbases(record::BamRecord, recorddata::HTSReadData)

Get total number of aligned bases in an alignment
"""
@inline alignedbases(record::BamRecord, recorddata::HTSReadData) = sum(!iszero, recorddata.alignmap)

"""
    alignmentlength(record::BamRecord, recorddata::HTSReadData)

Get alignment length = reference span + number of insertions
"""
@inline alignmentlength(record::BamRecord, recorddata::HTSReadData) = referencespan(record, recorddata) + numinsertions(record)


"""
    alignmentnumcorrect(record::BamRecord, rec)

Get the number of correct bases in an alignment
"""
@inline alignmentnumcorrect(record, recorddata::StencillingData{AuxMapModFireQC}) = alignmentlength(record, recorddata) - editdistance(record, recorddata)


"""
    alignmentaccuracy(record::BamRecord, recorddata::StencillingData{AuxMapModFireQC})

Get the alignment accuracy = num_correct/alignmentlength
"""
@inline function alignmentaccuracy(record, recorddata::StencillingData{AuxMapModFireQC})
    al = alignmentlength(record, recorddata)
    ed = editdistance(record, recorddata)
    
    (al - ed)/al
end


"""
    alignmentidentity(record::BamRecord, recorddata::StencillingData{AuxMapModFireQC})

Get the alignment identity = num_correct/alignedbases
"""
@inline function alignmentidentity(record, recorddata::StencillingData{AuxMapModFireQC})
    ab = alignedbases(record, recorddata)
    al = alignmentlength(record, recorddata)
    ed = editdistance(record, recorddata)
    
    (al - ed)/ab
end
