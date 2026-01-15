<a id="top"></a>

# SplineOP

## A Dynamic Programming Knot Selection Algorithm for Time-Series Compression With Quadratic Splines

### Nicolás Cecchi, Vincent Runge, Charles Truong and Laurent Oudre

- This project is the implementation of a paper submitted for peer-review. The paper is currently being handled by the editors of the journal. 

- There are plans for refactoring this code and integrating the method to the [ruptures library](https://pypi.org/project/ruptures/), stay tuned. 

> [Quick start](#start)

> [Models And Data Generators](#Models)

> [splineOP](#splineOP)

<a id="start"></a>

## Quick start

### Introduction

SplineOP finds changes in acceleration in time-series by using a Dynamic Programming algorithm. 


### Installing SplineOP for Python

**REQUIREMENTS:** - pybind11 scikit-build-core build ;

Currently we distribute the code for local compilation. As usual, it is recommended to work on a virtual environment.

You can create one by executing `python3 -m venv /path/to/env` and then activating it with `source /path/to/env/bin/activate`. (see [venv](https://docs.python.org/3/library/venv.html) for details)

Once there, you can install the dependencies with `pip install pybind11 scikit-build-core build`.

Finally, to build the package you git clone the repo and navigate to the /splineOPRCpp folder, build with `python3 -m build` and install with `pip install .` (see [scikit-build-core](https://scikit-build-core.readthedocs.io/en/latest/guide/getting_started.html), for details).

Once you have done that, you can import the package as following:

    import splineop_cpp as spop

### Example

Simple examples are provided in [this notebook](https://github.com/nicolascecchi/splineOPRcpp/blob/main/example.ipynb).

[(Back to Top)](#top)
