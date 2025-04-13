
To add a new command, implement a Tcl command function. For example,
the peridynamics commands are in `SRC/runtime/commands/modeling/peridynamics.cpp`.
These functions can be added to the command table in  `SRC/runtime/commands/modeling/commands.h`

```bash
python -m opensees tests/peridynamics.tcl
```
