from math import cos,sin,sqrt,pi
import opensees as ops
# Usage: 
#   OpenSees main.tcl <case>
#     where <case> is in {hwd1.1, hwd2.1, ...}
#
#   OpenSees main.tcl <case> - [<motion.zip>]
#     where <motion.zip> is an optional path to a zip file
#     with ground motion data.
#

#
# SETTINGS
#
# auto_path.append("/home/claudio/brace/Scripts")
config = (csmip_dir) "/home/claudio/brace/CSMIP/painter/"
default_motion = "config(csmip_dir)/RioDell_Petrolia_Processed_Data.zip"
# 
default_case = "hwd1.1"

#
# COMMAND LINE
#
if {argc > 0} {
  case = [lindex argv 0]  
else:
  case = default_case


if {argc > 1} {
  brace2_main = [lindex argv 1]


if {argc > 2} {
  motion_file = [lindex argv 2]
else:
  motion_file = default_motion


#
# SETUP
#
# package require brace2
source {C:\Users\16507\Documents\GitHub\Scripts\OpenSeesScripts\brace2.tcl}


# turn off case-specific analysis routines
# brace2_main = 1
source "{case}.tcl"


if {[info 'exists', brace2_main]} {
  print -JSON "dataDir/modelDetails.json"

  # Eigenvalue analysis                          output file              num modes
  [brace2::new 'EigenvalueAnalysis]', ana.analyze(file=dataDir, )/ModeShapes.yaml    3

  ana.algorithm('Newton')
  ana.constraints('Plain')
  ana.system('ProfileSPD')


  ana.rayleigh(){*}[brace2::rayleigh_alpha {1 0.03} {2 0.03}]

  #wipeAnalysis
  #brace2::new GravityAnalysis ga
  #ga analyze 10
  # Eigenvalue analysis                          output file              num modes
  #[brace2::new EigenvalueAnalysis] analyze -file dataDir/ModeShapes.yaml    3   -v

  #
  # RESPONSE HISTORY
  #
  #wipeAnalysis

  # RECORDERS
  # recorder Element -ele 4020 -xml dataDir/SectionDeformationHist.xml section 4 deformations
  # recorder Node -file dataDir/node705accel.txt -time -node 705 -dof 1 accel

  # EXCITATION
  brace2::new 'ResponseHistory', 'rh'
  rh 'patterns', {
    #            DOF (Ground motion zip file) Chn   (scale from cm to inch/s^2)   (rotation)
    UniformQuake 1      "motion_file"        20       -s 0.3937007874015748   --   -r 20
    UniformQuake 2      "motion_file"        20       -s 0.3937007874015748



  rh 'print'
  #         [steps]  [dt]
  #rh analyze 25; #1000; #10000 ; # 0.005

