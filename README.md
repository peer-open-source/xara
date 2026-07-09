

<p align="center">
  <a href="https://xara.so/">
    <img src="https://raw.githubusercontent.com/STAIRLab/xara-docs/master/source/_static/images/xara-chi.png" alt="xara logo" width="150" >
  </a>
</p>

<hr>

<!--
<h3 align="center">xara</h3>
-->


<p align="center">
  Nonlinear finite element analysis.
  <br>
</p>



<img align="left" src="https://raw.githubusercontent.com/claudioperez/sdof/master/docs/assets/peer-black-300.png" width="150px" alt="PEER Logo">


<br>

<div style="align:center">

[![Latest PyPI version](https://img.shields.io/pypi/v/xara?logo=pypi)](https://pypi.python.org/pypi/xara)
[![](https://img.shields.io/conda/v/opensees/opensees?color=%23660505)](https://anaconda.org/opensees/opensees)
[![PyPI Downloads](https://img.shields.io/pypi/dm/opensees)](https://pypi.org/project/opensees)
[![Mentioned in Awesome Finite Element Method](https://awesome.re/mentioned-badge.svg)](https://github.com/tkoyama010/awesome-finite-elements)

</div>

*xara* is a Python package that provides an intuitive and performant API for nonlinear
finite element analysis, implemented in C++ through the OpenSeesRT framework. 
OpenSees features state-of-the-art finite element formulations and solution 
algorithms, including mixed formulations for beams and solids, over 200 material models, and an
extensive collection of continuation algorithms to solve highly nonlinear
problems. 


This package may be used as a drop-in replacement for both `OpenSees.exe` and
OpenSeesPy (see *Getting Started* below), and generally provides a substantial performance boost.


### Getting Started

The `xara` package can be installed into a Python environment
in the standard manner. For example, using `pip`:

```shell
pip install xara
```

There are several ways to use the `xara` package:

- The standard way to use `xara` from Python is to create a `xara.Model`,
  and invoke its methods:
  ```python
  model = xara.Model(ndm=2, ndf=2) # Create a 2D model with 2 DOFs per node
  model.node(1, 2.0, 3.0)
  ...
  ```
  Most of the functions from the OpenSeesPy library can be invoked directly
  as methods of a `xara.Model` without any changes in syntax, although it is
  generally encouraged to use the safer variants supported by `xara`. For
  example:
  ```python
  # BAD
  model.pattern("Plain", 1, "Linear")
  model.load(1, 2.0, 3.0)
  # GOOD
  model.pattern("Plain", 1, "Linear")
  model.load(1, 2.0, 3.0, pattern=1)
  ```


- To start an interactive Tcl interpreter run the shell command:

  ```bash
  python -m opensees
  ```
  To quit the interpreter, just run `exit`:
  ```tcl
  opensees > exit
  ```

- To execute Tcl procedures programmatically from Python, create an instance
  of the `xara.Model` class and call its `eval()` method to evaluate Tcl code:
  ```python
  model = xara.Model()
  model.eval("model Basic -ndm 2")
  model.eval("print -json")
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

To compile the project see [about/compiling](https://xara.so/user/guides/compile.html)


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

