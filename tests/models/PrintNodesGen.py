from math import cos,sin,sqrt,pi
import opensees as ops
NodeCor =             [open, dataDir/ModeShape/Nodes.txt, w];
foreach, 'tag', [getNodeTags] {
  puts NodeCor [
    format "%d %lf %lf %lf" tag [nodeCoord tag 1] [nodeCoord tag 2] [nodeCoord tag 3]
  ]

close, NodeCor