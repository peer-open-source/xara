
### Example 1 a -  Nonlinear dynamic analysis using Portal Frame 

### Example 1.1 - Basic Truss 

- 2d 3 Element Elastic Truss
- Single Nodal Load, Static Analysis

Objectives:
- Simple Introduction to OpenSees


### Example 2.1 - Moment-Curvature 

- Zero length element with fiber section
- Single Nodal Load, Static Analysis

Objectives:
- Moment-Curvature Analysis in OpenSees


### Example 3.1 - Portal Frame 

- Reinforced concrete one-bay, one-story frame
- Distributed vertical load on girder


Objectives:
- Nonlinear beam-column elements
- Gravity load analysis and eigenvalue analysis



### Example 3.2 - Portal Frame 

- Reinforced concrete one-bay, one-story frame
- Distributed vertical load on girder
- Lateral Load at top of frame


### Example 3.1 -  Nonlinear pushover analysis using Portal Frame  as starting point


### Example 3.3 - Portal Frame 

- Reinforced concrete one-bay, one-story frame
- Distributed vertical load on girder
- Uniform excitation acting at fixed nodes in horizontal direction




### Example 4.1 - 2 Story Multi Bay Frame 

- Reinforced concrete multi-bay, two-story frame with distributed load on girder

Objectives:
- Nonlinear beam-column elements
- Gravity load analysis followed by pushover analysis
- Demonstrate scripting for the algorithmic level


### Example 5.1 - 3 Story One-by-One Bay Frame 

- Reinforced concrete one-bay, three-story frame with distributed load on girder

Objectives:
- 3D building with rigid diaphragms
- Nonlinear beam-column elements
- Gravity load analysis followed by transient analysis


### Example 6.1 - Solid Simply Supported Beam 

- Simply supported beam modeled with two dimensional solid elements

Objectives:
- test different quad elements
- free vibration analysis starting from static deflection


### Example 7.1 - 3D Shell Structure 

- Shell roof modeled with three dimensional linear shell elements

Objectives:

- test linear-elastic shell element
- free vibration analysis starting from static deflection


### Example 8.1 - Brick Cantilever Beam

- Cantilever beam modeled with three dimensional brick elements
  
Objectives:
-  test different brick elements
-  free vibration analysis starting from static deflection


## Other

The example subdirectories include:

        Example1:  contains a main() C++ routine which creates and
	    performs the analysis of the example1 problem outlined
	    in 'Adding an Element to G3'

                main.C - example C++ main driver (could be simpler)


               typing 'make' builds the executable example1
               typing 'make wipe' removes example1 and .o files

	ExampleScripts:  contains TCL scripts which define and analyze
	    R/C frames. to run these examples enter the interpreter
            and source in these files with the command:
		              source <fileName>

		RCFrame1.tcl - pushover analysis
		RCFrame2.tcl - pushover analysis with a simple
			TCL procedure that performs the analysis
		RCFrame3.tcl - linear elastic analysis
		RCFrame4.tcl - linear and nonlinear EQ analysis
		RCFrame5.tcl - pushover analysis with frame defined
			using a node and element generating TCL
			procedure
                RCFrameDisplay.tcl - commands for opening a window
	                used to display the RCFrame1-5 examples
	        RigidFrame3D.tcl - a nonlinear EQ analysis of a 3d frame
                        with elastic columns and nonlinear beams with rigid
                        diaphragm constraints
		genPlaneFrame.tcl - defines node and element generating
			TCL procedure


        PlaneFrame:  contains a main() C++ routine which uses
		the PlaneFrame model builder and performs a linear
	        static analysis.

                main.C - example C++ main driver (could be simpler)
	        quick.in - an example PlaneFrame input file.

               typing 'make' builds the executable PlaneFrame
               typing 'make wipe' removes PlaneFrame and .o files


        TclModelBuilder: contains the files outlined in the document
		'Adding an Element to g3'.

                myG3 - an executable.
                myCommands.C - contains procedure for adding new commands to g3
                TclPlaneTruss.h/.C - an example of a tcl model builder
                MyTruss.h/.C - files to add the new truss member
		example1.tcl, example2.tcl - the two tcl scripts
		     presented in the document 'Adding an Element to g3'

		examples - a directory containing some more scripts

                to load a script in the running interpreter type
                    source testX.tcl

                typing 'make' builds the executable myG3
                typing 'make wipe' remove myG3 and .o files


