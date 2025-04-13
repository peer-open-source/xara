set inFilename ExtEleForcesFloorsRIgravityForceForce.out

set numFloor 21
set numCline 6


if [catch {open $inFilename r} inFileID] {;	            # Open the input file and check for error
	puts stderr "Cannot open $inFilename for reading";	# output error statement
} else {

	foreach line [split [read $inFileID] \n] {; 	    # Look at each line in the file
	
	if {[llength $line] == 0} {;	                    # Blank line --> do nothing
	continue;
	} else {
        set timestep [expr [lindex $line 0]]
		
		if {$timestep == 1.0} {;                        # read the last step of gravity analysis result
		
			set counter 0
		
			for {set colLine 1} {$colLine <= $numCline} {incr colLine [expr $numCline-1]} {
				for {set floor1 1; set floor2 2} {$floor1 < $numFloor} {incr floor1 1; incr floor2 1} {
				
					set Ngrav [expr [lindex $line [expr 2+$counter*6]]]
								
					set counter [expr $counter+1]
				
					lappend HSSgrav $Ngrav        ;     # create list axial forces in external columns
					
				}
			}
		} else {
		}
		
	}
	}

close $inFileID; ;	# Close the input file
}

puts "axial list = $HSSgrav"


for {set colLine 1} {$colLine <= $numCline} {incr colLine 1} {
	    for {set floor1 1; set floor2 2} {$floor1 < $numFloor} {incr floor1 1; incr floor2 1} {
		
			if {$colLine == 1 || $colLine == $numCline} {
				
				set axial [expr [lindex $HSSgrav [expr $floor1 -1]]] 
				puts $axial
			}

		}
	}





