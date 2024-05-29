# 
# Use with: vmd -dispdev text -e <nameofscript.tcl>
# Original by Jordi Faraudo March 2017 modified by David 2020

#
# PROCEDURE 1: BIG DCD
#First of all I include here the BigDCD command
#

proc bigdcd { script type args } {
    global bigdcd_frame bigdcd_proc bigdcd_firstframe vmd_frame bigdcd_running
  
    set bigdcd_running 1
    set bigdcd_frame 0
    set bigdcd_firstframe [molinfo top get numframes]
    set bigdcd_proc $script

    # backwards "compatibility". type flag is omitted.
    if {[file exists $type]} { 
        set args [linsert $args 0 $type] 
        set type auto
    }
  
    uplevel #0 trace variable vmd_frame w bigdcd_callback
    foreach dcd $args {
        if { $type == "auto" } {
            mol addfile $dcd waitfor 0
        } else {
            mol addfile $dcd type $type waitfor 0
        }
    }
    after idle bigdcd_wait
}

proc bigdcd_callback { tracedvar mol op } {
    global bigdcd_frame bigdcd_proc bigdcd_firstframe vmd_frame
    set msg {}
 
    # If we're out of frames, we're also done 
    # AK: (can this happen at all these days???). XXX
    set thisframe $vmd_frame($mol)
    if { $thisframe < $bigdcd_firstframe } {
        puts "end of frames"
        bigdcd_done
        return
    }
 
    incr bigdcd_frame
    if { [catch {uplevel #0 $bigdcd_proc $bigdcd_frame} msg] } { 
        puts stderr "bigdcd aborting at frame $bigdcd_frame\n$msg"
        bigdcd_done
        return
    }
    animate delete beg $thisframe end $thisframe $mol
    return $msg
}

proc bigdcd_done { } {
    global bigdcd_running
    
    if {$bigdcd_running > 0} then {
        uplevel #0 trace vdelete vmd_frame w bigdcd_callback
        puts "bigdcd_done"
        set bigdcd_running 0
    }
}

proc bigdcd_wait { } {
    global bigdcd_running bigdcd_frame
    while {$bigdcd_running > 0} {
        global bigdcd_oldframe
        set bigdcd_oldframe $bigdcd_frame
        # run global processing hooks (including loading of scheduled frames)
        display update ui
        # if we have read a new frame during then the two should be different.
        if { $bigdcd_oldframe == $bigdcd_frame } {bigdcd_done}
    }
}

proc average L {
    expr ([join $L +])/[llength $L].
}

#
# Define procedure for counting aminoacids near PEG (resname PEG)
#
proc contaramino { frame } {
  global fp mol

  #aminoacids 3 letter cod    
  set amino {ARG HIS LYS ASP GLU SER THR ASN GLN CYS GLY PRO ALA VAL ILE LEU MET PHE TYR TRP}
    
  puts "$frame"

#  puts $fp "$frame" 

	
   set sel1 [atomselect $mol {resname PEG and type OG301 and within 3.5 of resname ARG}]
   set sel2 [atomselect $mol {resname PEG and type OG301 and within 3.5 of resname HIS}]
   set sel3 [atomselect $mol {resname PEG and type OG301 and within 3.5 of resname LYS}]
   set sel4 [atomselect $mol {resname PEG and type OG301 and within 3.5 of resname ASP}]
   set sel5 [atomselect $mol {resname PEG and type OG301 and within 3.5 of resname GLU}]
   set sel6 [atomselect $mol {resname PEG and type OG301 and within 3.5 of resname SER}]
   set sel7 [atomselect $mol {resname PEG and type OG301 and within 3.5 of resname THR}]
   set sel8 [atomselect $mol {resname PEG and type OG301 and within 3.5 of resname ASN}]
   set sel9 [atomselect $mol {resname PEG and type OG301 and within 3.5 of resname GLN}]
   set sel10 [atomselect $mol {resname PEG and type OG301 and within 3.5 of resname CYS}]
   set sel12 [atomselect $mol {resname PEG and type OG301 and within 3.5 of resname GLY}]
   set sel13 [atomselect $mol {resname PEG and type OG301 and within 3.5 of resname PRO}]
   set sel14 [atomselect $mol {resname PEG and type OG301 and within 3.5 of resname ALA}] 
   set sel15 [atomselect $mol {resname PEG and type OG301 and within 3.5 of resname VAL}]
   set sel16 [atomselect $mol {resname PEG and type OG301 and within 3.5 of resname ILE}]
   set sel17 [atomselect $mol {resname PEG and type OG301 and within 3.5 of resname LEU}]
   set sel18 [atomselect $mol {resname PEG and type OG301 and within 3.5 of resname MET}]
   set sel19 [atomselect $mol {resname PEG and type OG301 and within 3.5 of resname PHE}]
   set sel20 [atomselect $mol {resname PEG and type OG301 and within 3.5 of resname TYR}]
   set sel21 [atomselect $mol {resname PEG and type OG301 and within 3.5 of resname TRP}]
    
    
   $sel1 frame $frame
   $sel2 frame $frame
   $sel3 frame $frame
   $sel4 frame $frame
   $sel5 frame $frame
   $sel6 frame $frame
   $sel7 frame $frame
   $sel8 frame $frame
   $sel9 frame $frame
   $sel10 frame $frame
   $sel12 frame $frame
   $sel13 frame $frame
   $sel14 frame $frame
   $sel15 frame $frame
   $sel16 frame $frame
   $sel17 frame $frame
   $sel18 frame $frame
   $sel19 frame $frame
   $sel20 frame $frame
   $sel21 frame $frame
    
   $sel1 update
   $sel2 update
   $sel3 update
   $sel4 update
   $sel5 update
   $sel6 update
   $sel7 update
   $sel8 update
   $sel9 update
   $sel10 update
   $sel12 update
   $sel13 update
   $sel14 update
   $sel15 update
   $sel16 update
   $sel17 update
   $sel18 update
   $sel19 update
   $sel20 update
   $sel21 update  
    
   set n1 [$sel1 num]
   set n2 [$sel2 num] 
   set n3 [$sel3 num]
   set n4 [$sel4 num] 
   set n5 [$sel5 num]
   set n6 [$sel6 num] 
   set n7 [$sel7 num]
   set n8 [$sel8 num] 
   set n9 [$sel9 num]
   set n10 [$sel10 num] 
   set n12 [$sel12 num] 
   set n13 [$sel13 num]
   set n14 [$sel14 num] 
   set n15 [$sel15 num]
   set n16 [$sel16 num]
   set n17 [$sel17 num]
   set n18 [$sel18 num] 
   set n19 [$sel19 num]
   set n20 [$sel20 num] 
   set n21 [$sel21 num] 
    
   puts $fp "$frame $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 $n9 $n10 $n12 $n13 $n14 $n15 $n16 $n17 $n18 $n19 $n20 $n21"
      
   $sel1 delete
   $sel2 delete
   $sel3 delete
   $sel4 delete
   $sel5 delete
   $sel6 delete
   $sel7 delete
   $sel8 delete
   $sel9 delete
   $sel10 delete
   $sel12 delete
   $sel13 delete
   $sel14 delete
   $sel15 delete
   $sel16 delete
   $sel17 delete
   $sel18 delete
   $sel19 delete
   $sel20 delete
   $sel21 delete       
 
    
}

#
#Main program
#


# output file: PLEASE MODIFY NAME
set fp [ open "../data/contacts_aas_O.txt" w ] 

# open structure file of simulated system PLEASE MODIFY
set mol [mol new ../../input/ionized.psf type psf waitfor all]

#set list n1 {}

#perform calculation over each frame using BigDCD
puts "Please wait. Calculating..."
bigdcd contaramino auto ../../output3/equilNpT_stride10.dcd
bigdcd_wait

#set an1 [average $n1]
#puts "$an1"

#close output file          
close $fp 

exit

