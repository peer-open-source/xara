
- improved `printA`:
  - does not require `FullGen` solver; now all one needs
    to do is call `model.getTangent()`. This is possible due to the
    new consolidated memory ownership model.
  - `-mck` options

- improved `printA` command

- improved `print` command
  - new `-registry` option.
  - more reliable JSON printing
  - includes `MP_Constraint` information

- add `-det` option to static integrators:
  - ArcLength, ...

- add `-det` capability to solvers:
  - `FullGenLapack`, `Umfpack`, `BandGenLapack`

- Verbosity control

- new `export` command in Python and Tcl, when run through Python
- new `getTangent()` method in Python
- new `invoke` Tcl command and Python constructs
- new `progress` command in Tcl
- new `=` command, fixes vexing operator precedence in `expr`

