set starttime [clock seconds]
#GET NON-CODING SEQUENCES

if {![file exists NCDS]} {
    file mkdir NCDS
}

#open file
set in [open bac_anno.txt r]
set s_anno [read $in]
close $in

set s_anno [split $s_anno \n]

#open genome and sort it
set in_g [open bac_seq.txt r]
set full_g [read $in_g]
close $in_g

#find ORIGIN 
set origin [string first "ORIGIN " $full_g]
set ln_break [string first "\n" $full_g [expr $origin +10]]
set raw_seq [string range $full_g $ln_break end]

regsub -all {[^A-Za-z]} $raw_seq {} clean_seq
set seq [string tolower $clean_seq]



#output each non-coding sequence to a file
set sRNA_cnt 0
foreach sRNA $s_anno {
    
    if {[regexp {GenBank	sRNA} $sRNA]} {
    incr sRNA_cnt
    #get numbers
    regexp {GenBank	sRNA	([0-9]+?)[^0-9]+?([0-9]+?)[^0-9]} $sRNA all s_frst s_lst

    #get 'gene'
    if {[regexp {;gene=(.+?);} $sRNA]} {
    regexp {;gene=(.+?);} $sRNA all s_num} else {set s_num S$sRNA_cnt}

    #get product name
    if {[regexp {;product=(.+?);} $sRNA]} {
    regexp {;product=(.+?);} $sRNA all s_product} else {set s_product NA}

    #get locus tag
    if {[regexp {;LocusTag=(.+?);} $sRNA]} {
    regexp {;LocusTag=(.+?);} $sRNA all s_loc_tag} else {set s_loc_tag NA}

    #get non-coding sequence
    set srna_seq [string range $seq [expr $s_frst -1] [expr $s_lst -1]]

    if {[regexp {[^atcg]} $srna_seq]} {
        set okay 1
        puts $out_error "$accn\t$loc_tag\tactg"
    }

    #write to file
    if {![file exists NCDS/$s_num.tfa]} {
    set s_out [open NCDS/$s_num.tfa w]
    puts $s_out ">$s_num:product=$s_product\n$srna_seq"
    close $s_out
    }
    if {[expr $sRNA_cnt%100] == 0} {puts "Analysing $sRNA_cnt of 1548 NCDS..."}
    }
    
}

puts "NON-CODING SEQUENCES DONE"
