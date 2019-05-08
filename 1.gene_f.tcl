set starttime [clock seconds]
if {![file exists CDS]} {
file mkdir CDS}

set out_error [open errors.xls w]
puts $out_error "Accn\tloc_tag\terror"

set flist [glob -directory test_genomes *embl]

set num_genomes [llength $flist]
set gn_cnt 0

foreach fl $flist {

    set gn_cnt [expr $gn_cnt +1]
    
    set fname [file tail $fl]
    set accn [file rootname $fname]
    
    #see if you can put in a % thing...
    puts "Analysing $accn: $gn_cnt of $num_genomes"
    
    #careful with the slashes make sure you know which one wait on code from hurst
    if {![file exists CDS/$accn]} {
    file mkdir CDS/$accn}
    
    #open file
    set in [open $fl r]
    set data [read $in]
    close $in
    
    #get sequence
    set firstX [string first "XX\nSQ" $data]
    set ln_break [string first "\n" $data [expr $firstX +10]]
    set raw_seq [string range $data $ln_break end]
    
    regsub -all {[^A-Za-z]} $raw_seq {} clean_seq
    set seq [string tolower $clean_seq]
    
    #get annotation
    set top_bit [string range $data 0 [expr $firstX -1]]
    unset data
    
    set numFTCDS [regsub -all "FT   CDS" $top_bit "£FT   CDS" ntop_bit]
    
    if {$numFTCDS < 20} {
    
    set numFTgene [regsub -all "FT   gene" $top_bit "£FT   gene" top_bit]
    } else {set top_bit $ntop_bit}
    
    set gene_list [split $top_bit £]
    
    set real_top [lindex $gene_list 0]
    
    set gene_list [lrange $gene_list 1 end]
    
    set cnt 0
    
    #go get them
    foreach gn $gene_list {
    
    set cnt [expr $cnt +1]
    
    #does it say complement on the first line
    regexp {(.*?\n)} $gn all firstl
    set comp [regexp {complement} $firstl]
    
    #get the numbers
    set bits [split $firstl ,]
    set g_seq ""
    
    foreach exon $bits {
    append exon " "
    regexp {[^0-9]([0-9]+?)\.\.([0-9]+?)[^0-9]} $exon all frst lst
    set exon_seq [string range $seq [expr $frst -1] [expr $lst -1]]
    append g_seq $exon_seq
    }
    
    #reverse if complement
    if {$comp == 1} {
    set g_seq [string reverse $g_seq]
    regsub -all g $g_seq 1 g_seq
    regsub -all c $g_seq g g_seq
    regsub -all 1 $g_seq c g_seq
    regsub -all t $g_seq 1 g_seq
    regsub -all a $g_seq t g_seq
    regsub -all 1 $g_seq a g_seq
    }
    
    #locus tag
    if {[regexp {/locus_tag"(.+?)"} $gn]} {
    regexp {/locus_tag"(.+?)"} $gn all loc_tag} else {set loc_tag $cnt}
    
    #gene name
    if {[regexp {/gene="(.+?)"} $gn]} {regexp {/gene="(.+?)"} $gn all gname} elseif {
    [regexp {/product="(.+?)"} $gn]} {regexp {/product="(.+?)"} $gn all gname} else {
    set gname NA}
    
    #translation table
    if {[regexp {/transl_table=([0-9]+?)[^0-9]} $gn]} {
    regexp {/transl_table=([0-9]+?)[^0-9]} $gn all trans_tab} else {set trans_tab NA}
    
    #check for errors
    #mult of 3
    set okay 0
    set glen [string length $g_seq]
    set rem [expr $glen % 3]
    
    if {$rem !=0} {
    set okay 1
    puts $out_error "$accn\t$loc_tag\tx3"
    }
    
    #ends with a stop
    set last3 [string range $g_seq end-2 end]
    if {![regexp {tga|tag|taa} $last3]} {
    set okay 1
    puts $out_error "$accn\t$loc_tag\tnstop:$last3"
    }
    
    #funny bases
    if {[regexp {[^atcg]} $g_seq]} {
    set okay 1
    puts $out_error "$accn\t$loc_tag\tactg"
    }
    
    #write to file (careful on slash again) 
    set out [open CDS/$accn/$loc_tag.tfa w]
    puts $out ">$loc_tag:tt=$trans_tab:gen=$gname:comp=$comp:$accn\n$g_seq"
    close $out




}




}

close $out_error
set endtime [clock seconds]

puts "$num_genomes genomes done in [expr $endtime - $starttime] seconds!"







