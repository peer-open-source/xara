To run multiple models simultaneously, use `openseespy.Model(...)` (capital `M`)
instead of the regular `openseespy.model(...)` (lowercase m) function, and invoke
all subsequent modeling functions (e.g. `node(...)`, `element(...)`, `fix(...)`, etc)
as methods on the object returned from `Model()` instead of the `openseespy` submodule
directly. For example, instead of:

```python
ops.model("basic", "-ndm", 2, "-pdf", 3)
ops.node(1, 2, 3)
``` 
do
```python
model = ops.Model("basic", "-ndm", 2, "-ndf", 3)
model.node(1, 2, 3)
```
where the `openseespy` is assumed to have been imported 
with the `ops` alias:
```python
import opensees.openseespy as ops
```
