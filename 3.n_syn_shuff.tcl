source testcds_ab.tcl
source shuff_test.tcl
set flist [glob *tfa]

set out_plus [open plus1N.csv w]
set out_min [open min1N.csv w]
set out_zero [open zeroN.csv w]

set nums ""
for {set i 1} {$i < 1001} {incr i} {
    append nums ",$i"
}

set header "NCDS,Length$nums"

set raw_data_plus ""
set raw_data_min ""
set raw_data_zero ""

set gn_cnt 0
set total_genes [llength $flist]
puts "Analyzing gene 1 of $total_genes..."

set testl [string range $flist 0 20]

foreach fl $flist {
         
    incr gn_cnt
    set fname [file tail $fl]
    set accn [file rootname $fname]
    if {[expr $gn_cnt%100]==0} {puts "Analyzing gene $gn_cnt of $total_genes..."}

    #open file
    set in [open $fl r]
    set data [read $in]
    close $in

    #get and clean seq
    set newline [string first "\n" $data]
    set rawseq [string range $data $newline end]
    regsub -all {[^A-Za-z]} $rawseq {} seq

    set seqlen [string length $seq]

    for {set i 0} {$i < 1000} {incr i} {

        append raw_data_plus "\n$accn,$seqlen"
        append raw_data_min "\n$accn,$seqlen"
        append raw_data_zero "\n$accn,$seqlen"

        set randseq shuffle $seq

        set plus [string range $randseq 1 end]
        set min "x$randseq"
        set zero $randseq
        set plus_c [regexp -all -inline {.{1,3}} $plus]
        set min_c [regexp -all -inline {.{1,3}} $min
        set zero_c [regexp -all -inline {.{1,3}} $zero]
        set trans_p ""
        set trans_m ""
        set trans_z ""

        foreach codon $plus_c {
            if {$codon=="tag" | $codon=="taa" | $codon=="tga"} {append trans_p "X"} else {append trans_p "N"}
        }
        foreach codon $min_c {
            if {$codon=="tag" | $codon=="taa" | $codon=="tga"} {append trans_m "X"} else {append trans_p "N"}
        }
        foreach codon $zero_c {
            if {$codon=="tag" | $codon=="taa" | $codon=="tga"} {append trans_z "X"} else {append trans_p "N"}
        }

        set num_plus_stops [findstops [translate $trans_p]]
        set num_min_stops [findstops [translate $trans_m]]
        set num_zero_stops [findstops [translate $trans_z]]
        append raw_data_plus ",$num_plus_stops"
        append raw_data_min ",$num_min_stops"
        append raw_data_zero ",$num_zero_stops"
    }
    
}

puts $out_plus "$header $raw_data_plus"
puts $out_min "$header $raw_data_min"
puts $out_zero "$header $raw_data_zero"

close $out_plus
close $out_min
close $out_zero

puts FINISHED