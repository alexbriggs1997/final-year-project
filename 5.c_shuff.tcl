source testcds_ab.tcl
set flist [glob *tfa]

set out_plus [open plus1C.csv w]
set out_min [open min1C.csv w]

set nums ""
for {set i 1} {$i < 1001} {incr i} {
    append nums ",$i"
}

set header "CDS,Length$nums"

set raw_data_plus ""
set raw_data_min ""

set gn_cnt 0
set total_genes [llength $flist]
puts "Analyzing gene 1 of $total_genes..."

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

    append raw_data_plus "\n$accn,$seqlen"
    append raw_data_min "\n$accn,$seqlen"

    #do 1000 syn randomisations
    for {set i 0} {$i < 1000} {incr i} {
        
        translate $seq
        randomiser $transseq

        set plus_rand [string range $randomseq 1 end]
        set min_rand "x$randomseq"

        set num_plus_stops [findstops [translate $plus_rand]]
        set num_min_stops [findstops [translate $min_rand]]

        append raw_data_plus ",$num_plus_stops"
        append raw_data_min ",$num_min_stops"
    }
}

puts $out_plus "$header $raw_data_plus"
puts $out_min "$header $raw_data_min"

close $out_plus
close $out_min

puts FINISHED