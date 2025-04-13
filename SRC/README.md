
All files are in the following subdirectories:

<dl>

<dt><a href="./api"><code>api</code></a></dt>
<dd>
- `internal/` internal C API for extending
- `external/` implementation of Tcl commands
</dd>

<dt><a href="./element"><code>element</code></a></dt>
<dd>Contains the element classes, all in subdirectories, i.e. the
`Truss.h` and `Truss.C` files are in a subdirectory `truss/`.
</dd>

<dt><a href="./domain"><code>domain</code></a></dt>
<dd>contains the classes for the domain. all in subdirectores:

- `domain/` the Domain class and abstract Iters for the objects of
           the Domain, two subdirectories single- contains iters  
           for Domain, partitioned contains the PartitionedDomain
           class and the iters for this class.
- `node/` the `Node` and `NodalLoad` classes.
- `constraints/` the `SP_Constraint` and `MP_Constraint` classes.
- `subdomain/` the `Subdomain`, `ActorSubdomain` and `ShadowSubdoma` in classes.
- `pattern/` the abstact classes `LoadPattern` and `TimeSeries` and some concrete classes.
- `partitioner/` the DomainPartitioner class.
- `loadBalancer/` the LoadBalancer class and some concrete subclasses.

</dd>

<dt><a href="./analysis"><code>analysis</code></a></dt>
<dd>contains the classes for the analysis. again the classes in subdirectories:

 - `analysis/` the Analysis, StaticAnalysis, and some others
 - `algorithm/` SolutionAlgorithm class and 2 subdirectories:
           equiSolnAlgo - Linear, NewtonRaphson and ModifiedNewton
           domainDecompAlgo - DomainDecompAlgo
 - `dof_grp/` DOF_Group
 - `fe_ele/` `FE_Element` class and subdirectory penalty containing
           `PenaltySP_FE` and PenaltyMP_FE
 - `penalty/` ConstraintHandler, PenaltyConstraintHandler and
           PlainHandler
 - `integrator/` Integrator, IncrementalIntegrator, LoadControl,
           ArcLength, StaticIntegrator, TransientIntegrator and
           Newmark
 - `model/` AnalysisModel and its iters.
 - `numberer/` DOF_Numberer
           
<dt><a href="./graph"><code>graph</code></a></dt>
<dd>contains the graph classes, again in subdirectories

- `graph/` - `Graph`, `Vertex`, and others
- `partitioner/` - `GraphPartitioner` and `Metis`, also in metis-2.0
          contains the code downloaded to build metis.
- `numberer/` - `GraphNumberer`, `RCM` and some others

</dd>

<dt><a href="./system_of_eqn"><code>system_of_eqn</code></a></dt>
<dd>contains the SystemofEqn class and linearSOE subclass:
        linearSOE: contains LinearSOE, LinearSolver and DomainSolver
        classes and a bunch of subdirectories for each soe and solver:
        
- `fullGEN/` for a full general solver
- `bandGEN/` for a banded general solver
- `profileSPD/` for my profile solver and a solver using skypack,
           developed by O.Marques, now working with Jim Demmel
           (the code for skypack in a subdirectory)
- `petsc/` for the petsc solver
- `sparseGEN/` for superLU and thraeded superLU.
- `symSparse/` for Kincho Law's symmetric sparse solver

</dd>

<dt><a href="./tagged"><code>tagged</code></a></dt>
<dd>

`TagggedObject` and it's container classes (in subdirectory `storage/`).
This includes the abstract `TaggedObjectStorage`, the implementation
`ArrayOfTaggedObjects` and their iterators
These store the objects in `Domain` and `AnalysisModel`.

</dd>
 


<dt><a href="./runtime"><code>runtime</code></a></dt>
<dd>contains stuff for the interpreter.
</dd>


<dt><a href="./material"><code>material</code></a></dt>
<dd>Material, UniaxialMaterial and some concrete classes need by
        [FMK]'s Truss ele and Filip's element.
        
<dt><a href="./actor"><code>actor</code></a></dt>
<dd>The classes for [FMK]'s parallel stuff, again in subdirectories:

 - `shadow/` the Shadow class.  
 - `actor/` the Actor class
 - `channel/` the `Channel`, `TCP_Socket` and `MPI_Channel` classes
 - `machineBroker/` the `MachineBroker` and some classes for the machines over here.
 - `objectBroker/` the `FEM_ObjectBroker` class.
 - `message/` the `Message` class
 - `address/` the `Address` class.

<dt><a href="./matrix"><code>matrix</code></a><dd>contains the Matrix, Vector and ID classes.</dd>

</dl>

