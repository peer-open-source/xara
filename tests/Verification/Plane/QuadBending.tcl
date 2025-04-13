
set thk 0.01
set L 5.0
set H 1.0

set P 100.
set E 200.e6
set nu 0.3
set nodes {
  {  1 0. 0. }
  {  2 5. 0. }
  {  3 5. 1. }
  {  4 0. 1. }
  {  5 2.5 0. }
  {  6 5. .5 }
  {  7 2.5 1. }
  {  8 0. .5 }
  {  9 2.5 .5}
};



foreach elem {{quad8n 8} {quad9n 9} {quad 4}} {
  puts "$elem"
  set eleType [lindex $elem 0];
  set numNode [lindex $elem 1];

  wipe
  model Basic -ndm 2 -ndf 2
  foreach node [range $numNode] {node {*}[lindex $nodes $node]}

  nDMaterial ElasticIsotropic 1 $E $nu

# element quad9n 1 1 2 3 4 5 6 7 8 9 $thk "PlaneStress" 1
# element quad8n 1 1 2 3 4 5 6 7 8 $thk "PlaneStress" 1
  element $eleType 1 {*}[getNodeTags] $thk "PlaneStress" 1

  fix 1 1 1
  fix 4 1 0
# fix 8 1 0

  timeSeries Linear 1

  pattern Plain 1 1 {
          load 2 $P 0.
          load 3 -$P 0.
  }

  analysis Static

  analyze 1

  print -node 2 3

  puts {
  # verification:
  #   tip vertical displacement (node 2 and 3) = 0.0075
  #   bottom Gauss Point stress_xx = 46475.8
  #   bottom extrem stress_xx (extrapolated) = 60000.0
  }
}
