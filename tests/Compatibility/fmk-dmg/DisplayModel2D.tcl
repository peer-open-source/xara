proc DisplayModel2D { {ShapeType nill} {dAmp 5}  {xLoc 10} {yLoc 10} {xPixels 512} {yPixels 384} {nEigen 1} } {
	######################################################################################
	## DisplayModel2D $ShapeType $dAmp $xLoc $yLoc $xPixels $yPixels $nEigen
	######################################################################################
	## display Node Numbers, Deformed or Mode Shape in 2D problem
	##			Silvia Mazzoni & Frank McKenna, 2006
	##
	## 	ShapeType : 	type of shape to display. # options: ModeShape , NodeNumbers , DeformedShape 
	## 	dAmp : 		relative amplification factor for deformations
	## 	xLoc,yLoc  : 	horizontal & vertical location in pixels of graphical window (0,0=upper left-most corner)
	## 	xPixels,yPixels :	width & height of graphical window in pixels
	## 	nEigen : 		if nEigen not=0, show mode shape for nEigen eigenvalue
	##	
	#######################################################################################
	global TunitTXT;					# load time-unit text
	global ScreenResolutionX ScreenResolutionY;	# read global values for screen resolution

	if {  [info exists TunitTXT] != 1} {set TunitTXT ""};		# set blank if it has not been defined previously.

	if {  [info exists ScreenResolutionX] != 1} {set ScreenResolutionX 1024};		# set default if it has not been defined previously.
	if {  [info exists ScreenResolutionY] != 1} {set ScreenResolutionY 768};		# set default if it has not been defined previously.

	if {$xPixels == 0} {
		set xPixels [expr int($ScreenResolutionX/2)];		
		set yPixels [expr int($ScreenResolutionY/2)]
		set xLoc 10
		set yLoc 10
	}
	if {$ShapeType == "nill"} {
		puts ""; puts ""; puts "------------------"
		puts "View the Model? (N)odes, (D)eformedShape, anyMode(1),(2),(#). Press enter for NO."
		gets stdin answer
		if {[llength $answer]>0 } { 
			if {$answer != "N" & $answer != "n"} {
				puts "Modify View Scaling Factor=$dAmp? Type factor, or press enter for NO."
				gets stdin answerdAmp
				if {[llength $answerdAmp]>0 } { 
					set dAmp $answerdAmp
				}
			}
			if {[string index $answer 0] == "N" || [string index $answer 0] == "n"} {
				set ShapeType NodeNumbers
			} elseif {[string index $answer 0] == "D" ||[string index $answer 0] == "d" } {
				set ShapeType DeformedShape
			} else {
				set ShapeType ModeShape
				set nEigen $answer
			}
		} else {
			return
		}
	}

	if {$ShapeType ==  "ModeShape" } {
		set lambdaN [eigen $nEigen];		# perform eigenvalue analysis for ModeShape
		set lambda [lindex $lambdaN [expr $nEigen-1]];
		set omega [expr pow($lambda,0.5)]
		set PI 	[expr 2*asin(1.0)];		# define constant
		set Tperiod [expr 2*$PI/$omega];	   	# period (sec.) 
		set fmt1 "Mode Shape, Mode=%.1i Period=%.3f %s  "
		set windowTitle [format $fmt1 $nEigen $Tperiod  $TunitTXT]
	} elseif  {$ShapeType ==  "NodeNumbers" } {
		set windowTitle "Node Numbers"
	} elseif  {$ShapeType ==  "DeformedShape" } {
		set windowTitle "Deformed Shape"
	}

	set viewPlane XY
	recorder display $windowTitle $xLoc $yLoc $xPixels $yPixels  -wipe ; # display recorder
	DisplayPlane $ShapeType $dAmp $viewPlane $nEigen 0
}
