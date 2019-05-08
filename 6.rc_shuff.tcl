source testcds_ab.tcl
set flist [glob -directory New *tfa]
###############
set testlist [lrange $flist 0 25]
###############
set rout_plus [open rplusC1.csv w]
set rout_min [open rmin1C.csv w]

set header "CDS,Raw Seq"
set raw_data_plus ""
set raw_data_min ""
set gn_cnt 0
foreach fl $flist {
    incr gn_cnt
    set fname [file tail $fl]
    set accn [file rootname $fname]
    if {[expr $gn_cnt%100]==0} {puts "Analyzing gene $gn_cnt of bare..."}

    set in [open $fl r]
    set data [read $in]
    close $in

    set newline [string first "\n" $data]
    set rawseq [string range $data $newline end]
    regsub -all {[^A-Za-z]} $rawseq {} seq

    set seqlen [string length $seq]

    append raw_data_plus "\n$accn"
    append raw_data_min "\n$accn"

    set plus [string range $seq 1 end]
    set min "x$seq"

    set codons_p [regexp -all -inline {.{1,3}} $plus]
    set codons_m [regexp -all -inline {.{1,3}} $min]

    set trans_p ""
    set trans_m ""

    foreach codon $codons_p {
        if {$codon=="tag" | $codon=="taa" | $codon=="tga"} {append trans_p "X"} else {append trans_p "N"}
    }

    foreach codon $codons_m {
        if {$codon=="tag" | $codon=="taa" | $codon=="tga"} {append trans_m"X"} else {append trans_m "N"}
    }

    set num_plus_stops [findstops $trans_p]
    set num_min_stops [findstops $trans_m]

    append raw_data_plus ",$num_plus_stops"
    append raw_data_min ",$num_min_stops"

}

puts $rout_plus "$header $raw_data_plus"
puts $rout_min "$header $raw_data_min"

close $rout_plus
close $rout_min

puts FINISHED