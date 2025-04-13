# `opensees`

<img align="left" src="https://raw.githubusercontent.com/claudioperez/sdof/master/docs/assets/peer-black-300.png" width="150px" alt="PEER Logo">

Nonlinear finite element analysis.

<br>

<div style="align:center">

[![Latest PyPI version](https://img.shields.io/pypi/v/opensees?logo=pypi)](https://pypi.python.org/pypi/opensees)
[![](https://img.shields.io/conda/v/opensees/opensees?color=%23660505)](https://anaconda.org/opensees/opensees)
[![PyPI Downloads](https://img.shields.io/pypi/dm/opensees)](https://pypi.org/project/opensees)

</div>


`opensees` is a Python package that provides an intuitive API for nonlinear
finite element analysis, implemented in C++ through the OpenSees framework. 
OpenSees features state-of-the-art finite element formulations and solution 
algorithms, including mixed formulations for beams and solids, over 200 material models, and an
extensive collection of continuation algorithms to solve highly nonlinear
problems. 

The `opensees` package supports interactive post processing via the
[`sees`](https://pypi.org/project/sees) package.


The package may be used as a drop-in replacement for both `OpenSees.exe` and
OpenSeesPy (see *Getting Started* below), and generally provides a substantial performance boost.

<p style="text-align: center;">
<b>This package is <i>experimental</i> and not yet intended for public use.</b>
</p>


## Features

- **Performance** Switching Python scripts to use `opensees` typically results in a 4x to 5x performance boost.
- **Interactive Tasks**: Easily return stiffness, mass, and damping matrices as NumPy arrays and join meshes without duplicate nodes and constraints.
- **Extensive Modeling Library**: State-of-the-art element formulations with over 200 material models to choose from.
- **Continuation Algorithms**: Robust algorithms for solving highly nonlinear problems.
- **Intuitive and Reliable** The core OpenSees runtime has been redesigned so that all program 
  state is encapsulated in user-instantiated classes,
  and global variables/singletons are avoided. 
  This eliminates several preexisting vulnerabilities to inadvertent state corruption.


<!-- 
- **Semantics** Unlike interfaces which rely on global state, this package can be used 
  with true library semantics. 
-->

Additional features include:

- Convert OpenSeesPy scripts into equivalent Tcl files that can be used
  for faster processing or serialization. Unlike most conversion utilities,
  this conversion is done *exactly* and does not rely on hand-rolled parsing.

- The package can be installed with `pip` for Python versions 3.7 - 3.12 on Linux, MacOS and
  Windows, but eigenvalue analysis is currently broken on Windows.

> [!NOTE]
> This package is independent of the [`openseespy`](https://pypi.org/project/openseespy)
> library, which is documented in the OpenSees [documentation](https://opensees.github.io/OpenSeesDocumentation)
> website. OpenSeesPy can be installed by running the following command:
>
> ```shell
> pip install opensees[py]
> ```



### Getting Started

The `opensees` package can be installed into a Python environment
in the standard manner. For example, using `pip`:

```shell
pip install opensees
```

There are several ways to use the `opensees` package:

- To execute Tcl procedures from a Python script, just create an instance
  of the `opensees.tcl.Interpreter` class and call its `eval()` method:
  ```python
  interp = opensees.tcl.Interpreter()
  interp.eval("model Basic -ndm 2")
  interp.eval("print -json")
  ```

- To start an interactive interpreter run the shell command:

  ```bash
  python -m opensees
  ```
  To quit the interpreter, just run `exit`:
  ```tcl
  opensees > exit
  ```

- The `opensees` package exposes a compatibility layer that exactly reproduces
  the *OpenSeesPy* functions, but does so without mandating a single
  global program state. To run OpenSeesPy scripts, just change the import:

  ```python
  import openseespy.opensees
  ```
  to
  ```python
  import opensees.openseesrt
  ```
  For true stateless modeling, the `Model` class should be used instead of the legacy
  `model` function; documentation is under development.


## Development

To compile the project see [help/compiling](https://github.com/claudioperez/opensees/blob/master/help/compiling.md)

<!-- Badge links -->

[pypi-d-image]: https://img.shields.io/pypi/dm/opensees.svg
[license-badge]: https://img.shields.io/pypi/l/opensees.svg
[pypi-d-link]: https://pypi.org/project/opensees
[pypi-v-image]: https://img.shields.io/pypi/v/opensees.svg
[pypi-v-link]: https://pypi.org/project/opensees


## See also

- [`osmg`](https://pypi.org/project/osmg) OpenSees Model Generator
- [`sees`](https://pypi.org/project/sees) Modern rendering library
- [`mdof`](https://pypi.org/project/mdof) Optimized system identification library
- [`sdof`](https://pypi.org/project/sdof) Optimized integration for single degree of freedom systems

For more projects by the STAIRlab, visit https://github.com/STAIRlab .

## Support

<table align="center" style="border: 0;">

 <tr style="background-color:rgba(0, 0, 0, 0);">
  <td style="background-color:rgba(0, 0, 0, 0);" >
    <a href="https://peer.berkeley.edu">
    <img src="https://raw.githubusercontent.com/claudioperez/sdof/master/docs/assets/peer-black-300.png"
         alt="PEER Logo" width="200"/>
    </a>
  </td>

  <td>
    <a href="https://dot.ca.gov/">
    <img src="https://raw.githubusercontent.com/claudioperez/sdof/master/docs/assets/Caltrans.svg.png"
         alt="Caltrans Logo" width="200"/>
    </a>
  </td>

  <td>
    <!-- <a href="https://brace2.herokuapp.com"> -->
    <a href="https://stairlab.berkeley.edu">
    <img src="https://raw.githubusercontent.com/claudioperez/sdof/master/docs/assets/stairlab.svg"
         alt="STAIRlab Logo" width="200"/>
    </a>
  </td>
 
 </tr>
</table>

