# src/splineOPRcpp/__init__.py

# Import the compiled extension module (the .so/.pyd file) 
# and expose its functions/classes.
try:
    from .splineOPRcpp import *
except ImportError as e:
    raise ImportError(f"Could not load C++ extension: {e}")

__version__ = "0.1.0"
