# read external data file 
set inFilename Test_xq1/FBeam.out

if [catch {open $inFilename r} inFileID] {;	# Open the input file and check for error
puts stderr "Cannot open $inFilename for reading";	# output error statement
} else {

	set outFileID [open output.txt w]

	foreach line [split [read $inFileID] \n] { ;	# Look at each line in the file
	#puts "line = $line"
	if {[llength $line] == 0} {;	# Blank line --> do nothing
	continue;
	} else {
	puts $outFileID $line
	set Xvalues $line;	# execute operation on read data
	set col1 [expr [lindex $line 2]]
	}
	}
close $outFileID;
close $inFileID; ;	# Close the input file
}

puts $Xvalues
puts $col1
