BEGIN {
    OFS="\t"
}
{
    contig = $5
    start = $6 - 1
    end = $7
    strand = ($9 == "C") ? "-" : "+"
    split($10, a, "#")
    repeat_class = a[2]
    print contig, start, end, repeat_class, ".", strand
}