proc shuffle {str} {
set lst [split $str {}]

     
     set len [llength $lst]
     set len2 $len
     for {set i 0} {$i < $len-1} {incr i} {
         set n [expr {int($i + $len2 * rand())}]
         incr len2 -1

         # Swap elements at i & n
         set temp [lindex $lst $i]
         lset lst $i [lindex $lst $n]
         lset lst $n $temp
        }
        set lst [join $lst ""]
         return $lst
}