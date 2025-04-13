# TripleFrictionPendulum

This command is used to construct a Triple Friction Pendulum Bearing
(TPB) (Figure 1) element object. The element is a 3-dimensional element
with variable <a
href="http://opensees.berkeley.edu/wiki/index.php/FrictionModel_Command">friction
coefficient models</a> [1] [2]. The element accounts for the
vertical-horizontal coupling and the bidirectional coupling in
horizontal behavior. The friction coefficient model is a general model
that accounts for the variation of friction coefficient on velocity and
vertical force. Other simplified friction coefficient models such as
vertical-force-independent friction coefficient, velocity-independent
friction coefficient and constant friction coefficient can also be
defined. The element can also be used for modeling single friction
pendulum bearings or double friction pendulum bearings by simplifying
the general backbone curve of the TPB.

```tcl
element TripleFrictionPendulum $eleTag $iNode $jNode
        $frnTag1 $frnTag2 $frnTag3 $vertMatTag $rotZMatTag $rotXMatTag
        $rotYMatTag $L1 $L2 $L3 $d1 $d2 $d3 $W $uy $kvt $minFv
        $tol
```

<hr />

<table>
<tbody>
<tr class="odd">
<td><code class="parameter-table-variable">eleTag</code></td>
<td><p>= unique element object tag</p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">iNode jNode</code></p></td>
<td><p>= end nodes</p></td>
</tr>
<tr class="odd">
<td><p><strong>frnTag1, frnTag2, frnTag2</code></p></td>
<td><p>= tags associated with previously-defined <a
href="http://opensees.berkeley.edu/wiki/index.php/FrictionModel_Command">FrictionModels</a>
at the three sliding interfaces</p></td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">vertMatTag</code></td>
<td><p>= Pre-defined material tag for COMPRESSION behavior of the
bearing</p></td>
</tr>
<tr class="odd">
<td><p><code>rotZMatTag, rotXMatTag, rotYMatTag</code></p></td>
<td><p>= Pre-defined material tags for rotational behavior about 3-axis,
1-axis and 2-axis, respectively.</p></td>
</tr>
<tr class="even">
<td><p><code>L1, L2, L3</code></p></td>
<td><p>= effective radii. Li = R_i - h_i (see Figure 1)</p></td>
</tr>
<tr class="odd">
<td><p><code>d1, d2, d3</code></p></td>
<td><p>= displacement limits of pendulums (Figure 1). Displacement limit
of the bearing is 2d1+d2+d3+L1.d3/L3-L1.d2/L2</p></td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">W</code></td>
<td><p>= axial force used for the first trial of the first analysis
step.</p></td>
</tr>
<tr class="odd">
<td><code class="parameter-table-variable">uy</code></td>
<td><p>= lateral displacement where sliding of the bearing starts.
Recommended value = 0.25 to 1 mm. A smaller value may cause convergence
problem.</p></td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">kvt</code></td>
<td><p>= Tension stiffness k_vt of the bearing.</p></td>
</tr>
<tr class="odd">
<td><p><code>minFv</code> (&gt;=0)</p></td>
<td><p>= minimum vertical compression force in the bearing used for
computing the horizontal tangent stiffness matrix from the normalized
tangent stiffness matrix of the element. `minFv` is substituted for the
actual compressive force when it is less than `minFv`, and prevents the
element from using a negative stiffness matrix in the horizontal
direction when uplift occurs. The vertical nodal force returned to nodes
is always computed from `kvc` (or `kvt`) and vertical deformation, and
thus is not affected by `minFv`.</p></td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">tol</code></td>
<td><p>= relative tolerance for checking the convergence of the element.
Recommended value = 1.e-10 to 1.e-3.</p></td>
</tr>
</tbody>
</table>

<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Fig1.jpg" title="TPB_Nhan_Fig1.jpg" width="500"
alt="TPB_Nhan_Fig1.jpg" />
<figcaption aria-hidden="true">TPB_Nhan_Fig1.jpg</figcaption>
</figure>

<p>The horizontal normalized behavior of the element is an extension of
the unidirectional behavior proposed by Fenz and Constantinou [3] and
Morgan and Mahin [4]. The envelope normalized backbone curve for
unidirectional behavior is in Figure 2 where effective radii are
calculated from L_i=R_i- h_i, based on the bearing geometry in Figure 1.
Displacements u_i^* and normalized forces f_i^* are evaluated according
to [2], [3] or [4]. <span style="color:blue"> The Excel file for
generating a backbone curve of a TPB can be downloaded
here:</span><a href="Media:TPB_Nhan_BackboneCurve.xls"
title="wikilink">Media:TPB_Nhan_BackboneCurve.xls</a></p>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Fig2.jpg" title="TPB_Nhan_Fig2.jpg" width="400"
alt="TPB_Nhan_Fig2.jpg" />
<figcaption aria-hidden="true">TPB_Nhan_Fig2.jpg</figcaption>
</figure>
<p>Overturning moment and torsion due to the eccentricity of internal
forces are equally distributed to the 2 nodes of the element.</p>

In the vertical direction, the element is multi-linear elastic with
different stiffnesses k_vc and k_vt in compression and tension,
respectively. Even though a TPB has no resistance in tension, a small
nonzero stiffness should be provided in tension for stability of the
numerical procedure. The reasonable vertical stiffness in tension
depends on the stiffness of the superstructure, but a value between 10
N/m to 100 N/m should work well in most cases. Very small ratio of
k_vt/k_vc may cause the convergence difficulty when the superstructure
is uplifted.

<p><strong>SPECIAL CASES:</strong></p>
<dl>
<dt>Vertical-horizontal uncoupled bearing:</dt>
<dd>
  Restrain vertical displacement of the 2 nodes and set `minFv` = static
  vertical reaction of the bearing.
</dd>

<dt>Neglecting rotational stiffness:</dt>
<dd>
  Rotational stiffness of the bearing can be neglected by defining a very
  small value for the `rotXMatTag`, `rotYMatTag` and `rotZMatTag`. However,
  using a too small number may cause a numerical convergence problem.
</dd>
</dl>

<p><strong>NOTES:</strong></p>
- The current element requires 6 degree of freedoms at each node and
  defines the local coordinate system to be the same as the global
  coordinate system, where the vertical axis must be 3.

- Since the element accounts for the vertical-horizontal coupling
  behavior of TPBs, the time step size in analysis of vertically stiff
  structures should be small enough so that the high frequency components
  in responses can be captured.

- Rayleigh damping is automatically included when using Rayleigh command.

- The height of the bearing (for computing overturning moment from
  horizontal force) is computed from the vertical distance between the two
  end nodes.


## Examples

```tcl
frictionModel VelNormalFrcDep 1 [expr 0.012/pow($W,0.8-1.0)] 0.8 \
        [expr 0.018/pow($W,0.7-1.0)] 0.7 25.0 0.0 0.0 3.0

frictionModel VelNormalFrcDep 2 [expr 0.052/pow($W,0.8-1.0)] 0.8 \
        [expr 0.075/pow($W,0.7-1.0)] 0.7 25.0 0.0 0.0 3.0

frictionModel VelNormalFrcDep 3 [expr 0.12/pow($W,0.8-1.0)] 0.8 \
        [expr 0.16/pow($W,0.7-1.0)] 0.7 25.0 0.0 0.0 3.0

uniaxialMaterial Elastic 1 1.e6
uniaxialMaterial Elastic 2 100.;
uniaxialMaterial Elastic 3 100.;
uniaxialMaterial Elastic 4 10.;
element TripleFrictionPendulum 1 1 2 1 2 3 1 4 2 3 0.36 1.25 1.25 \
       0.1 0.2 0.2 1000.0 0.0005 1.0 0.1 1.E-5;
```

<hr />

Download these 3 ground motion files and change the extension ".tcl"
to ".ATH" for EXAMPLES 2 to 7: <a href="Media:TCU065-E.tcl"
title="wikilink">Media:TCU065-E.tcl</a>, <a href="Media:TCU065-N.tcl"
title="wikilink">Media:TCU065-N.tcl</a>, <a href="Media:TCU065-V.tcl"
title="wikilink">Media:TCU065-V.tcl</a>

### Example 1
<p>Unidirectional static pushover of a TPB element with constant
friction coefficients.</p>
<p>Tcl code: <a
href="Media:CheckingTripleFrictionPendulum_StaticCyclicPushover.tcl"
title="wikilink">Media:CheckingTripleFrictionPendulum_StaticCyclicPushover.tcl</a></p>

<p>Results:</p>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex1a.png" title="TPB_Nhan_Ex1a.png" width="400"
alt="TPB_Nhan_Ex1a.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex1a.png</figcaption>
</figure>

### Example 2
Unidirectional dynamic seismic analysis of a single mass supported by
a TPB element with constant friction coefficients.
<p>Tcl code: <a
href="Media:CheckingTripleFrictionPendulum_1DDynamics_ConstFriction.tcl"
title="wikilink">Media:CheckingTripleFrictionPendulum_1DDynamics_ConstFriction.tcl</a></p>
<p>Results:</p>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex2_1Disp.png" title="TPB_Nhan_Ex2_1Disp.png"
width="500" alt="TPB_Nhan_Ex2_1Disp.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex2_1Disp.png</figcaption>
</figure>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex2_2DispTrace.png"
title="TPB_Nhan_Ex2_2DispTrace.png" width="250"
alt="TPB_Nhan_Ex2_2DispTrace.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex2_2DispTrace.png</figcaption>
</figure>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex2_3Reaction.png" title="TPB_Nhan_Ex2_3Reaction.png"
width="500" alt="TPB_Nhan_Ex2_3Reaction.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex2_3Reaction.png</figcaption>
</figure>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex2_4Loop.png" title="TPB_Nhan_Ex2_4Loop.png"
width="500" alt="TPB_Nhan_Ex2_4Loop.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex2_4Loop.png</figcaption>
</figure>

```tcl
EXAMPLE 3:
```
<p>Unidirectional dynamic seismic analysis of a single mass supported by
a TPB element with variable friction coefficients.</p>
<p>Tcl code: <a
href="Media:CheckingTripleFrictionPendulum_1DDynamics_VariableFriction.tcl"
title="wikilink">Media:CheckingTripleFrictionPendulum_1DDynamics_VariableFriction.tcl</a></p>
<p>Results:</p>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex3_1Disp.png" title="TPB_Nhan_Ex3_1Disp.png"
width="500" alt="TPB_Nhan_Ex3_1Disp.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex3_1Disp.png</figcaption>
</figure>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex3_2DispTrace.png"
title="TPB_Nhan_Ex3_2DispTrace.png" width="250"
alt="TPB_Nhan_Ex3_2DispTrace.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex3_2DispTrace.png</figcaption>
</figure>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex3_3Reaction.png" title="TPB_Nhan_Ex3_3Reaction.png"
width="500" alt="TPB_Nhan_Ex3_3Reaction.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex3_3Reaction.png</figcaption>
</figure>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex3_4Loop.png" title="TPB_Nhan_Ex3_4Loop.png"
width="500" alt="TPB_Nhan_Ex3_4Loop.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex3_4Loop.png</figcaption>
</figure>

```tcl
EXAMPLE 4:
```
<p>Two-dimensional dynamic seismic analysis of a single mass supported
by a TPB element with constant friction coefficients.</p>
<p>Tcl code: <a
href="Media:CheckingTripleFrictionPendulum_2DDynamic_ConstFriction.tcl"
title="wikilink">Media:CheckingTripleFrictionPendulum_2DDynamic_ConstFriction.tcl</a></p>
<p>Results:</p>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex4_1Disp.png" title="TPB_Nhan_Ex4_1Disp.png"
width="500" alt="TPB_Nhan_Ex4_1Disp.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex4_1Disp.png</figcaption>
</figure>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex4_2DispTrace.png"
title="TPB_Nhan_Ex4_2DispTrace.png" width="250"
alt="TPB_Nhan_Ex4_2DispTrace.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex4_2DispTrace.png</figcaption>
</figure>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex4_3aReaction.png"
title="TPB_Nhan_Ex4_3aReaction.png" width="500"
alt="TPB_Nhan_Ex4_3aReaction.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex4_3aReaction.png</figcaption>
</figure>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex4_4Loop.png" title="TPB_Nhan_Ex4_4Loop.png"
width="500" alt="TPB_Nhan_Ex4_4Loop.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex4_4Loop.png</figcaption>
</figure>

```tcl
EXAMPLE 5:
```
<p>Two-dimensional dynamic seismic analysis of a single mass supported
by a TPB element with variable friction coefficients.</p>
<p>Tcl code: <a
href="Media:CheckingTripleFrictionPendulum_2DDynamic_VariableFriction.tcl"
title="wikilink">Media:CheckingTripleFrictionPendulum_2DDynamic_VariableFriction.tcl</a></p>
<p>Results:</p>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex5_1Disp.png" title="TPB_Nhan_Ex5_1Disp.png"
width="500" alt="TPB_Nhan_Ex5_1Disp.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex5_1Disp.png</figcaption>
</figure>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex5_2Trace.png" title="TPB_Nhan_Ex5_2Trace.png"
width="250" alt="TPB_Nhan_Ex5_2Trace.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex5_2Trace.png</figcaption>
</figure>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex5_3Reaction.png" title="TPB_Nhan_Ex5_3Reaction.png"
width="500" alt="TPB_Nhan_Ex5_3Reaction.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex5_3Reaction.png</figcaption>
</figure>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex5_4aLoop.png" title="TPB_Nhan_Ex5_4aLoop.png"
width="500" alt="TPB_Nhan_Ex5_4aLoop.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex5_4aLoop.png</figcaption>
</figure>

```tcl
EXAMPLE 6:
```
<p>Three-dimensional dynamic seismic analysis of a single mass supported
by a TPB element with constant friction coefficients.</p>
<p>Tcl code: <a
href="Media:CheckingTripleFrictionPendulum_3DDynamic_ConstFriction.tcl"
title="wikilink">Media:CheckingTripleFrictionPendulum_3DDynamic_ConstFriction.tcl</a></p>
<p>Results:</p>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex6_1Disp.png" title="TPB_Nhan_Ex6_1Disp.png"
width="500" alt="TPB_Nhan_Ex6_1Disp.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex6_1Disp.png</figcaption>
</figure>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex6_2DispTrace.png"
title="TPB_Nhan_Ex6_2DispTrace.png" width="250"
alt="TPB_Nhan_Ex6_2DispTrace.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex6_2DispTrace.png</figcaption>
</figure>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex6_3Reaction.png" title="TPB_Nhan_Ex6_3Reaction.png"
width="500" alt="TPB_Nhan_Ex6_3Reaction.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex6_3Reaction.png</figcaption>
</figure>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex6_4Loop.png" title="TPB_Nhan_Ex6_4Loop.png"
width="500" alt="TPB_Nhan_Ex6_4Loop.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex6_4Loop.png</figcaption>
</figure>

```tcl
EXAMPLE 7:
```
<p>Three-dimensional dynamic seismic analysis of a single mass supported
by a TPB element with variable friction coefficients.</p>
<p>Tcl code: <a
href="Media:CheckingTripleFrictionPendulum_3DDynamic_VariableFriction.tcl"
title="wikilink">Media:CheckingTripleFrictionPendulum_3DDynamic_VariableFriction.tcl</a></p>
<p>Results:</p>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex7_1Disp.png" title="TPB_Nhan_Ex7_1Disp.png"
width="500" alt="TPB_Nhan_Ex7_1Disp.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex7_1Disp.png</figcaption>
</figure>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex7_2DispTrace.png"
title="TPB_Nhan_Ex7_2DispTrace.png" width="250"
alt="TPB_Nhan_Ex7_2DispTrace.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex7_2DispTrace.png</figcaption>
</figure>
<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex7_3Reaction.png" title="TPB_Nhan_Ex7_3Reaction.png"
width="500" alt="TPB_Nhan_Ex7_3Reaction.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex7_3Reaction.png</figcaption>
</figure>

<figure>
<img src="/OpenSeesRT/contrib/static/TPB_Nhan_Ex7_4Loop.png" title="TPB_Nhan_Ex7_4Loop.png"
width="500" alt="TPB_Nhan_Ex7_4Loop.png" />
<figcaption aria-hidden="true">TPB_Nhan_Ex7_4Loop.png</figcaption>
</figure>

<hr />

## References
<p>[1] Dao ND, Ryan KL, Sato E, Sasaki T. Predicting the displacement of
triple pendulum™ bearings in a full-scale shaking experiment using a
three-dimensional element. Earthquake Engineering and Structural
Dynamics, 2013.</p>
<p>[2] Dao ND. Seismic Response of a Full-scale 5-story Steel Frame
Building Isolated by Triple Pendulum Bearings under 3D Excitations.
Dissertation, University of Nevada - Reno, 2012.</p>
<p>[3] Fenz DM, Constantinou MC. Spherical sliding isolation bearings
with adaptive behavior: Theory. Earthquake Engineering and Structural
Dynamics 2008; 37(2):163-183.</p>
<p>[4] Morgan TA, Mahin SA. The use of innovative base isolation systems
to achieve complex seismic performance objectives. PEER-2011/06
2011.</p>

<hr />
<p>Code developed by: <span style="color:blue"> Nhan D. Dao,
University of Nevada - Reno. E-mail: nhan.unr@gmail.com
</span>.</p>

