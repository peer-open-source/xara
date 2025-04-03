
<img align="left" src="https://raw.githubusercontent.com/claudioperez/sdof/master/docs/assets/peer-black-300.png" width="150px" alt="PEER Logo">

Nonlinear finite element analysis.

<br>

<div style="align:center">

[![Latest PyPI version](https://img.shields.io/pypi/v/opensees?logo=pypi)](https://pypi.python.org/pypi/opensees)
[![](https://img.shields.io/conda/v/opensees/opensees?color=%23660505)](https://anaconda.org/opensees/opensees)
[![PyPI Downloads](https://img.shields.io/pypi/dm/opensees)](https://pypi.org/project/opensees)

</div>

*xara* is a Python package that provides an intuitive and performant API for nonlinear
finite element analysis, implemented in C++ through the OpenSeesRT framework. 
OpenSees features state-of-the-art finite element formulations and solution 
algorithms, including mixed formulations for beams and solids, over 200 material models, and an
extensive collection of continuation algorithms to solve highly nonlinear
problems. 


This package may be used as a drop-in replacement for both `OpenSees.exe` and
OpenSeesPy (see *Getting Started* below), and generally provides a substantial performance boost.

<p style="text-align: center;">
<b>This package is <i>experimental</i> and not yet intended for public use.</b>
</p>


> [!NOTE]
> This package is independent of the [`openseespy`](https://pypi.org/project/openseespy)
> library, which is documented in the OpenSees [documentation](https://opensees.github.io/OpenSeesDocumentation)
> website. 


### Getting Started

The `xara` package can be installed into a Python environment
in the standard manner. For example, using `pip`:

```shell
pip install xara
```

There are several ways to use the `xara` package:

- To execute Tcl procedures from a Python script, just create an instance
  of the `xara.Model` class and call its `eval()` method:
  ```python
  model = xara.Model()
  model.eval("model Basic -ndm 2")
  model.eval("print -json")
  ```

- To start an interactive interpreter run the shell command:

  ```bash
  python -m opensees
  ```
  To quit the interpreter, just run `exit`:
  ```tcl
  opensees > exit
  ```

- The `xara` package exposes a compatibility layer that exactly reproduces
  the *OpenSeesPy* functions, but does so without mandating a single
  global program state. To run OpenSeesPy scripts, just change the import:

  ```python
  import openseespy.opensees
  ```
  to
  ```python
  import opensees.openseespy
  ```
  For true stateless modeling, the `Model` class should be used instead of the legacy
  `model` function; see the documentation [here](https://xara.so/user/manual/model/model_class.html).


## Development

To compile the project see [about/compiling](https://github.com/claudioperez/OpenSeesRT/blob/main/about/compiling.md)

<!-- Badge links -->

[pypi-d-image]: https://img.shields.io/pypi/dm/opensees.svg
[license-badge]: https://img.shields.io/pypi/l/opensees.svg
[pypi-d-link]: https://pypi.org/project/opensees
[pypi-v-image]: https://img.shields.io/pypi/v/opensees.svg
[pypi-v-link]: https://pypi.org/project/opensees


## Support

<table align="center" style="border: 0;">
<tr>
  <td>
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

 </tr>
</table>

