
set N 10 ; # Number of elements
set L 10 ; # Cantilever length

#               A,E,G,J,Iy,Iz
set properties {1 1 1 1  1  1}

model basic 3 6

transform Linear 1 {0 0 1}

foreach i [range 1 [= 2+$N]] x [linspace 0 $L [= 1+$N]]  {
  node $i $x 0 0;
}

foreach i [range 1 [= 1+$N]] {
  element ElasticBeamColumn $i $i [expr int($i+1)] {*}$properties 1
}

fix 1 ; # 1 1 1 1 1 1
print -json


