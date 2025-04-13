# Triple Friction Pendulum Bearing

<p>This command is used to construct a Triple Friction Pendulum Bearing
element object, which is defined by two nodes. The element can have zero
length or the appropriate bearing height. The bearing has unidirectional
(2D) or coupled (3D) friction properties (with post-yield stiffening due
to the concave sliding surface) for the shear deformations, and
force-deformation behaviors defined by `UniaxialMaterials` in the
remaining two (2D) or four (3D) directions. To capture the uplift
behavior of the bearing, the user-specified UniaxialMaterial in the
axial direction is modified for no-tension behavior. P-Delta moments are
entirely transferred to the concave sliding surface (iNode). It is
important to note that rotations of the concave sliding surface
(rotations at the iNode) affect the shear behavior of the bearing. If
the element has non-zero length, the local x-axis is determined from the
nodal geometry unless the optional x-axis vector is specified in which
case the nodal geometry is ignored and the user-defined orientation is
utilized.</p>
<figure>
<img src="/OpenSeesRT/contrib/static/TFP_backbone.gif" title="TFP_backbone.gif"
alt="TFP_backbone.gif" />
<figcaption aria-hidden="true">TFP_backbone.gif</figcaption>
</figure>

```tcl
element TFP $eleTag $iNode $jNode $R1 $R2 $R3 $R4 $D1 $D2 $D3 $D4 \
        $d1 $d2 $d3 $d4 $mu1 $mu2 $mu3 $mu4 $h1 $h2 $h3 $h4 $H0 $colLoad < $K >
```

<hr />
<table>
<tbody>
<tr class="odd">
<td><code class="parameter-table-variable">eleTag</code></td>
<td><p>unique element object tag</p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">iNode jNode</code></p></td>
<td><p>end nodes</p></td>
</tr>
<tr class="odd">
<td><p><code class="parameter-table-variable">R1</code></p></td>
<td><p>Radius of inner bottom sliding surface</p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">R2</code></p></td>
<td><p>Radius of inner top sliding surface</p></td>
</tr>
<tr class="odd">
<td><p><code class="parameter-table-variable">R3</code></p></td>
<td><p>Radius of outer bottom sliding surface</p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">R4</code></p></td>
<td><p>Radius of outer top sliding surface</p></td>
</tr>
<tr class="odd">
<td><p><code class="parameter-table-variable">D1</code></p></td>
<td><p>Diameter of inner bottom sliding surface</p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">D2</code></p></td>
<td><p>Diameter of inner top sliding surface</p></td>
</tr>
<tr class="odd">
<td><p><code class="parameter-table-variable">D3</code></p></td>
<td><p>Diameter of outer bottom sliding surface</p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">D4</code></p></td>
<td><p>Diameter of outer top sliding surface</p></td>
</tr>
<tr class="odd">
<td><p><code class="parameter-table-variable">d1</code></p></td>
<td><p>diameter of inner slider</p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">d2</code></p></td>
<td><p>diameter of inner slider</p></td>
</tr>
<tr class="odd">
<td><p><code class="parameter-table-variable">d3</code></p></td>
<td><p>diameter of outer bottom slider</p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">d4</code></p></td>
<td><p>diameter of outer top slider</p></td>
</tr>
<tr class="odd">
<td><p><code class="parameter-table-variable">mu1</code></p></td>
<td><p>friction coefficient of inner bottom sliding surface</p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">mu2</code></p></td>
<td><p>friction coefficient of inner top sliding surface</p></td>
</tr>
<tr class="odd">
<td><p><code class="parameter-table-variable">mu3</code></p></td>
<td><p>friction coefficient of outer bottom sliding surface</p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">mu4</code></p></td>
<td><p>friction coefficient of outer top sliding surface</p></td>
</tr>
<tr class="odd">
<td><p><code class="parameter-table-variable">h1</code></p></td>
<td><p>height from inner bottom sliding surface to center of
bearing</p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">h2</code></p></td>
<td><p>height from inner top sliding surface to center of
bearing</p></td>
</tr>
<tr class="odd">
<td><p><code class="parameter-table-variable">h3</code></p></td>
<td><p>height from outer bottom sliding surface to center of
bearing</p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">h4</code></p></td>
<td><p>height from inner top sliding surface to center of
bearing</p></td>
</tr>
<tr class="odd">
<td><code class="parameter-table-variable">H0</code></td>
<td><p>total height of bearing</p></td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">colLoad</code></td>
<td><p>initial axial load on bearing (only used for first time step then
load come from model)</p></td>
</tr>
<tr class="odd">
<td><code class="parameter-table-variable">K</code></td>
<td><p>optional, stiffness of spring in vertical dirn (dof 2 if ndm= 2,
dof 3 if ndm = 3) (default=1.0e15)</p></td>
</tr>
</tbody>
</table>
<figure>
<img src="/OpenSeesRT/contrib/static/TFPwH0.gif" title="TFPwH0.gif" alt="TFPwH0.gif" />
<figcaption aria-hidden="true">TFPwH0.gif</figcaption>
</figure>
<figure>
<img src="/OpenSeesRT/contrib/static/TFP_displaced.gif" title="TFP_displaced.gif"
alt="TFP_displaced.gif" />
<figcaption aria-hidden="true">TFP_displaced.gif</figcaption>
</figure>
<p>NOTE:</p>
<ol>
<li>If the element has zero length and optional orientation vectors are
not specified, the local element axes coincide with the global axes.
Otherwise the local z-axis is defined by the cross product between the
$x$- and $y$-vectors specified on the command line.</li>

<li>The valid queries to a triple friction pendulum bearing element when
creating an ElementRecorder object are `force`, `localForce`,
`basicForce`, `localDisplacement`, `basicDisplacement`, `relativeDisp`,
`plasticDisp`, and `material $matNum matArg1 matArg2 ...` Where `matNum`
is the number associated with the material whose data is to be output.
  <ol>
  <li><strong>relativeDisp</strong> returns relative displacements between
  the sliding components in the bearing. Relative displacements is the
  rotation (as shown in the figure above) multiplied by the respective
  radii. For each time step it returns 8 values; 4 for each horizontal
  direction.</li>
  <li><strong>plasticDisp</strong> returns plastic displacements
  associated with relativeDisp</li>
  </ol></li>
</ol>

## Examples

```tcl
element TFP 1 1 2 12.0 12.0 88.0 88.0 12.0 12.0 44.0 44.0 8.0 8.0 \
    12.5 12.5 0.02 0.02 0.09 0.12 3.0 3.0 4.5 4.5 12.5 45.0;
```

## References
<p>Becker, TC, Mahin, SA. "Experimental and analytical study of the
bi-directional behavior of the triple friction pendulum isolator,"
Earthquake Engineering and Structural Dynamics. (accepted for
publication 03/11)<a
href="http://onlinelibrary.wiley.com/doi/10.1002/eqe.1133/pdf">read the
paper</a></p>
<hr />
<p>Code Developed by: <span style="color:blue"> Tracy Becker,
University of California, Berkeley. </span></p>
