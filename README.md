<a id="top"></a>

# splineOP

## A Dynamic Programming Knot Selection Algorithm for Time-Series Compression With Quadratic Splines

### Nicolás Cecchi, Vincent Runge, Charles Truong and Laurent Oudre

> [Quick start](#start)

> [Models And Data Generators](#Models)

> [splineOP](#splineOP)

<a id="start"></a>

## Quick start

### Introduction

SplineOP finds changes in acceleration in time-series. 


### Installing SplineOP for Python

pip install 
python3 -m build
pip install .

**REQUIREMENTS:** - pybind11 scikit-build-core build ;

Currently we distribute the code for local compilation. As usual, it is recommended to work on a virtual environment.

You can create one by executing `pip -m venv /path/to/env` and then activating it with `source /path/to/env/bin/activate`.

Once there, you can install the dependencies with `pip install pybind11 scikit-build-core build`

Finally, to build the package you git clone the repo and navigate to the /splineOPRCpp folder, build with `python3 -m build` and install with `pip install .`.

Once you have done that, you can import the package as following:

    import splineop_cpp as spop

### A simple example

[(Back to Top)](#top)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

<a id="Models"></a>

## Models And Data Generators

### Quadratic Splines

### Data Generators

[(Back to Top)](#top)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

<a id="splineOP"></a>

## splineOP

[(Back to Top)](#top)
