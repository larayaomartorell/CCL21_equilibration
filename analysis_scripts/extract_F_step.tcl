### Open the log file for reading and the output .txt file for writing
set file [open "../../input_9/equilNpT_smd.log" r]
set output [open "../data/F_step_9.txt" w]

### Loop over all lines of the log file
while { [gets $file line] != -1 } {

   ### Determine if a line contains SMD output. If so, write the
   ### step number + coordinates + force vector
   
   if {[string range $line 0 3] == "SMD "} {
      puts $output "[expr [lindex $line 1]] [expr [lindex $line 2]] [expr [lindex $line 3]] [expr [lindex $line 4]] [expr [lindex $line 5]] [expr [lindex $line 6]] [expr [lindex $line 7]]"
   }

}

### Close the log file and the output file
close $file
close $output
