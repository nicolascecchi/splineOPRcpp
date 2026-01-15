import sys
from setuptools import setup, Extension
import setuptools.command.build_ext as _build_ext

# Check if the pybind11 module can be imported to get its include path
try:
    import pybind11
except ImportError:
    pass 

# Custom build_ext command to handle compiler flags
class build_ext(_build_ext.build_ext):
    def finalize_options(self):
        super().finalize_options()
        # Check if compiler supports C++17, which is often required for modern C++ and Pybind11 features
        for ext in self.extensions:
            ext.extra_compile_args = ['-std=c++17', '-Wall', '-O3']
            # Important: Since your Rcpp code uses Eigen, we need to ensure Eigen headers are found.
            # You may need to manually add the Eigen include path here if it's not global.
            # E.g., ext.include_dirs.append('/path/to/Eigen/include/')

# Define the C++ extension module
splineop_module = Extension(
    'splineop_py', # The name of the resulting Python module (you will 'import splineop_py')
    sources=['py_bindings.cpp'], # Your binding source file
    include_dirs=[
        # If using Method 1 (pip install), uncomment the line below:
        # pybind11.get_include(),
        # Add any directories containing your core C++ headers (SplineOP.h, etc.)
        '.', 
        # You may need to add the Eigen include path here
        # Example: '/usr/include/eigen3', 
    ],
    language='c++',
    # The modules need to be linked against the core libraries from your Rcpp build (if any).
    # Since SplineOP likely only relies on Eigen, usually only Eigen is needed here.
    libraries=[], 
)

setup(
    name='splineop_cpp',
    version='0.1.0',
    description='Python interface for SplineOP change point detection engine',
    author='Nicolás Cecchi',
    ext_modules=[splineop_module],
    cmdclass={'build_ext': build_ext},
    zip_safe=False,
)
