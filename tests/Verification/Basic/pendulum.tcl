model basic 2 3

set L 1
set mass 10
set g 9.8
set dt 0.001

node 1 0  0
# node 2 0 -$L -mass $mass $mass $mass

fix 1 1 1 0


# rigidLink beam 2 1


timeSeries Constant 1 -factor $g

pattern UniformExcitation 1 2 -accel 1

constraints Transformation
numberer Plain
algorithm Linear
integrator Newmark 0.5 [expr 1.0/6.0]
system BandGeneral
analysis Transient

foreach t [linspace 0 10 10000] {
  analyze 1 $dt
  puts "$t\t[nodeDisp 2]"
}

