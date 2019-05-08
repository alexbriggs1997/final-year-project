#Determine Processes

#frameshift procedure
proc frameshift {seq n} {

    set tempname [string range $seq $n end]
    for {set i 0} {$i < $n } {incr i} {
        append tempname "x"
    }


    return $tempname
}


#Declaring Global variables - codons used
set c_used {

}

#Runs through the sequence in incrs of 3 and translates into aa sequence via a switch
proc translate {seq} {
    global c_used ;
    set c_used {} 
    global transseq
    set transseq ""

    #global codon lists
    global cods_F
    global cods_L
    global cods_I
    global cods_M
    global cods_V
    global cods_S
    global cods_P
    global cods_T
    global cods_A
    global cods_Y
    global cods_X
    global cods_H
    global cods_Q
    global cods_N
    global cods_K
    global cods_D
    global cods_E
    global cods_C
    global cods_W
    global cods_R
    global cods_S
    global cods_G
    
    #Create lists for codons used
    set cods_F [list ]
    set cods_L [list ]
    set cods_I [list ]
    set cods_M [list ]
    set cods_V [list ]
    set cods_S [list ]
    set cods_P [list ]
    set cods_T [list ]
    set cods_A [list ]
    set cods_Y [list ]
    set cods_X [list ]
    set cods_H [list ]
    set cods_Q [list ]
    set cods_N [list ]
    set cods_K [list ]
    set cods_D [list ]
    set cods_E [list ]
    set cods_C [list ]
    set cods_W [list ]
    set cods_R [list ]
    set cods_S [list ]
    set cods_G [list ]

    for {set t 0} {$t < [string length $seq] } {set t [expr $t + 3]} {
        set codon [string range $seq $t [expr $t +2]]
        set iscodon true 
        switch $codon {
            ttt {append transseq "F"
                 lappend cods_F ttt}   
            ttc {append transseq "F"
                 lappend cods_F ttc}
            tta {append transseq "L"
                 lappend cods_L tta}
            ttg {append transseq "L"
                 lappend cods_L ttg}
            ctt {append transseq "L"
                 lappend cods_L cct}
            ctc {append transseq "L"
                 lappend cods_L ctc}
            cta {append transseq "L"
                 lappend cods_L cta}
            ctg {append transseq "L"
                 lappend cods_L ctg}
            att {append transseq "I"
                 lappend cods_I att}
            atc {append transseq "I"
                 lappend cods_I atc}
            ata {append transseq "I"
                 lappend cods_I ata}
            atg {append transseq "M"
                 lappend cods_M atg}
            gtt {append transseq "V"
                 lappend cods_V gtt}
            gtc {append transseq "V"
                 lappend cods_V gtc}
            gta {append transseq "V"
                 lappend cods_V gta}
            gtg {append transseq "V"
                 lappend cods_V gtg}
            tct {append transseq "S"
                 lappend cods_S tct}
            tcc {append transseq "S"
                 lappend cods_S tcc}
            tca {append transseq "S"
                 lappend cods_S tca}
            tcg {append transseq "S"
                 lappend cods_S tcg}
            cct {append transseq "P"
                 lappend cods_P cct}
            ccc {append transseq "P"
                 lappend cods_P ccc}
            cca {append transseq "P"
                 lappend cods_P cca}
            ccg {append transseq "P"
                 lappend cods_P ccg}
            act {append transseq "T"
                 lappend cods_T act}
            acc {append transseq "T"
                 lappend cods_T acc}
            aca {append transseq "T"
                 lappend cods_T aca}
            acg {append transseq "T"
                 lappend cods_T acg}
            gct {append transseq "A"
                 lappend cods_A gct}
            gcc {append transseq "A"
                 lappend cods_A gcc}
            gca {append transseq "A"
                 lappend cods_A gca}
            gcg {append transseq "A"
                 lappend cods_A gcg}
            tat {append transseq "Y"
                 lappend cods_Y tat}
            tac {append transseq "Y"
                 lappend cods_Y tac}
            taa {append transseq "X"
                 lappend cods_X taa}
            tag {append transseq "X"
                 lappend cods_X tag}
            cat {append transseq "H"
                 lappend cods_H cat}
            cac {append transseq "H"
                 lappend cods_H cac}
            caa {append transseq "Q"
                 lappend cods_Q caa}
            cag {append transseq "Q"
                 lappend cods_Q cag}
            aat {append transseq "N"
                 lappend cods_N aat}
            aac {append transseq "N"
                 lappend cods_N aac}
            aaa {append transseq "K"
                 lappend cods_K aaa}
            aag {append transseq "K"
                 lappend cods_K aag}
            gat {append transseq "D"
                 lappend cods_D gat}
            gac {append transseq "D"
                 lappend cods_D gac}
            gaa {append transseq "E"
                 lappend cods_E gaa}
            gag {append transseq "E"
                 lappend cods_E gag}
            tgt {append transseq "C"
                 lappend cods_C tgt}
            tgc {append transseq "C"
                 lappend cods_C tgc}
            tga {append transseq "X"
                 lappend cods_X tga}
            tgg {append transseq "W"
                 lappend cods_W tgg}
            cgt {append transseq "R"
                 lappend cods_R cgt}
            cgc {append transseq "R"
                 lappend cods_R cgc}
            cga {append transseq "R"
                 lappend cods_R cga}
            cgg {append transseq "R"
                 lappend cods_R cgg}
            agt {append transseq "S"
                 lappend cods_S agt}
            agc {append transseq "S"
                 lappend cods_S agc}
            aga {append transseq "R"
                 lappend cods_R aga}
            agg {append transseq "R"
                 lappend cods_R agg}
            ggt {append transseq "G"
                 lappend cods_G ggt}
            ggc {append transseq "G"
                 lappend cods_G ggc}
            gga {append transseq "G"
                 lappend cods_G gga}
            ggg {append transseq "G"
                 lappend cods_G ggg}
            default {set iscodon false}
        }
        if {"$iscodon"} {append c_used $codon} 
    } 

    return $transseq 
}

#Finds the stops by searching for X's and then incr the nstops
proc findstops {seq} {
    set nstops 0
    for {set s 0} {$s < [string length $seq]} {incr s} {
        if {[string index $seq $s] == "X"} {
            incr nstops 
        }
    } 
    return $nstops 
}

#Randomiser - Take the translated sequence and produce a randomised version at syn sites
#why is rawcodons an input? or is seq rawcodons one thing
#what do i need to globalise?
proc randomiser {seq} {

     global randomseq

    global cods_F
    global cods_L
    global cods_I
    global cods_M
    global cods_V
    global cods_S
    global cods_P
    global cods_T
    global cods_A
    global cods_Y
    global cods_X
    global cods_H
    global cods_Q
    global cods_N
    global cods_K
    global cods_D
    global cods_E
    global cods_C
    global cods_W
    global cods_R
    global cods_S
    global cods_G

    set randomseq ""
    for {set c 0} {$c < [string length $seq]} {incr c} {
        set aa [string index $seq $c] 
        set poscodons {}
        switch $aa {
            F {set poscodons  $cods_F}
            L {set poscodons  $cods_L}
            I {set poscodons  $cods_I}
            M {set poscodons  $cods_M}
            V {set poscodons  $cods_V}
            S {set poscodons  $cods_S}
            P {set poscodons  $cods_P}
            T {set poscodons  $cods_T}
            A {set poscodons  $cods_A}
            Y {set poscodons  $cods_Y}
            X {set poscodons  $cods_X}
            H {set poscodons  $cods_H}
            Q {set poscodons  $cods_Q}
            N {set poscodons  $cods_N}
            K {set poscodons  $cods_K}
            D {set poscodons  $cods_D}
            E {set poscodons  $cods_E}
            C {set poscodons  $cods_C}
            W {set poscodons  $cods_W}
            R {set poscodons  $cods_R}
            S {set poscodons  $cods_S}
            G {set poscodons  $cods_G}
            default {set poscodons [list "na"]}

        }
        set rancodonno [expr {int(rand() * [expr [llength $poscodons]-1])}]
        set rancodon [lindex $poscodons $rancodonno]
        append randomseq $rancodon

        #remove codon used from bin
        set idx [lsearch [set cods_$aa] $rancodon]
        set cods_$aa [lreplace [set cods_$aa] $idx $idx]
        
        
    }
    return $randomseq 
}


#Find the Genes and Open them
#set genelist [glob -directory Genes *.tfa]
#
#foreach filepath $genelist {
#    set fname [file tail $filepath]
#    set accn [file rootname $fname]
#
#    set file [open $filepath r]
#    set data [read $file]
#    close $file
#
#    #this is splitting the data by new line breaks then finding the 3rd from last element which is the seq
#
#    set splitdata [split $data \n] 
#    set seq [lindex $splitdata [expr [llength $splitdata] -3]]
#
##FRAMESHIFT BY +1 AND +2
#    set frameshift1 [frameshift $seq 1]
#    set frameshift2 [frameshift $seq 2]
#
#    set prot [translate $seq]
#    set prot_codons $c_used
#    set prot1 [translate $frameshift1]
#    set prot2 [translate $frameshift2]
#    set p1stops [findstops $prot1]
#    set p2stops [findstops $prot2]
#    set randomboi [randomiser $prot {}]
##Set out destination and puts all of the variables we want
#    set out [open E:/Tests/testcds/Genes/output/$fname w]
#    close $out
#
#    set out [open E:/Tests/testcds/Genes/output/$fname a]
#    puts $out "Original seq\n $seq"
#    puts $out "Frameshift +1\n $frameshift1"
#    puts $out "Frameshift +2\n $frameshift2"
#    puts $out "Translation\n $prot"
#    puts $out "Translation +1\n $prot1"
#    puts $out "Translation +2\n $prot2"
#    puts $out "Number of Stops +1 = $p1stops"
#    puts $out "Number of Stops +2 = $p2stops"
#    puts $out "Randomised sequence\n $randomboi"
#    close $out
#
#}
#


#to do
# FIX RANDOMISER done
# Once the randomiser works... make a bin that for each gene only selects the codons that exist within the original sequence.
#Create a for loop around randomboi that creates x number of randomboi and reports each stop it finds to a list that then makes average/min/max
#Change the output to a csv instead of my stupid test file
#Change the input directory to the actual genes 
#profit

#call lists codon_Q/G/Fetc then change codon_$aa