# written: fmk
# date: 02/99
#
# purpose: example1 in OpenSeesIntro.tex

#create the ModelBuilder object
model BasicBuilder -ndm 2 -ndf 2

# build the model 

# add nodes - command: node nodeId xCrd yCrd
node 1   0.0  0.0
node 2 144.0  0.0
node 3 168.0  0.0
node 4  72.0 96.0

# add material - command: uniaxialMaterial <matType> matID <matArgs>
uniaxialMaterial Elastic 1 3000

# add truss elements - command: element truss trussID node1 node2 A matID
element truss 1 1 4 10.0 1
element truss 2 2 4 5.0 1
element truss 3 3 4 5.0 1

# set the boundary conditions - command: fix nodeID xResrnt? yRestrnt?
fix 1 1 1 
fix 2 1 1
fix 3 1 1

pattern Plain 1 "Linear" {
    # apply the load - command: load nodeID xForce yForce
    load 4 100 -50
}

