"""
This script creates two distinct models (`system_1` and `system_2`)
using the `opensees.openseespy.Model` class, and analyzes them independently.
"""

import opensees.openseespy as ops


def make_truss(system):
    # create nodes & add to Domain - command: node nodeId xCrd yCrd
    system.node(1,   0.0,  0.0)
    system.node(2, 144.0,  0.0)
    system.node(3, 168.0,  0.0)
    system.node(4,  72.0, 96.0)

    # set the boundary conditions - command: fix nodeID xRestrnt? yRestrnt?
    system.fix(1, 1, 1)
    system.fix(2, 1, 1)
    system.fix(3, 1, 1)

    # Define materials for truss elements
    # -----------------------------------
    # Create Elastic material prototype - command: uniaxialMaterial Elastic matID E
    system.uniaxialMaterial("Elastic", 1, 3000.0)

    # Define elements
    # ---------------
    system.element("truss", 1, 1, 4, 10.0, 1)
    system.element("truss", 2, 2, 4,  5.0, 1)
    system.element("truss", 3, 3, 4,  5.0, 1)

    # Define loads
    # ------------
    # create a Linear TimeSeries (load factor varies linearly with time)
    system.timeSeries("Linear", 1)

    # create a Plain load pattern
    system.pattern("Plain", 1, 1, fact=1.0)
    # create the nodal load: load nodeID xForce yForce
    system.load(4, 100.0, -50.0)


def analyze(system):
    # ------------------------------
    # Analysis generation
    # ------------------------------
    # create the system of equation, a SPD using a band storage scheme
    system.system("BandSPD")

    # create the DOF numberer, the reverse Cuthill-McKee algorithm
    system.numberer("RCM")

    # create the constraint handler, a Plain handler is used as homo constraints
    system.constraints("Plain")

    # create the solution algorithm, a Linear algorithm is created
    system.algorithm("Linear")

    # create the integration scheme, the LoadControl scheme using steps of 1.0
    system.integrator("LoadControl", 1.0)

    # create the analysis object 
    system.analysis("Static")

    # perform the analysis
    system.analyze(1)


# create first independent model (with two-dimensions and 2 DOF/node)
system_1 = ops.Model("basic", "-ndm", 2, "-ndf", 2)

# create second independent model. Note the optional improved syntax 
# for `ndm` and `ndf`.
system_2 = ops.Model("basic", ndm=2, ndf=2)

# create truss models in each system
make_truss(system_1)
make_truss(system_2)

# Analyze system_1 and print its state
analyze(system_1)
system_1.printModel("node", 4)
system_1.printModel("ele")

# Print the state of system_2 to show that it has not
# been analyzed
system_2.printModel("node", 4)
system_2.printModel("ele")


analyze(system_2)
system_2.printModel("node", 4)
system_2.printModel("ele")

