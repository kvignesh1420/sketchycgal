"""
_libsketchycgal.so loader for importing into python.
"""
import os
import ctypes
import sys
import inspect
import importlib

def _load_libsketchycgal():
    """_load_library"""
    # f = inspect.getfile(sys._getframe(1))  # pylint: disable=protected-access

    # # Construct filename
    # f = os.path.join(os.path.dirname(f), "_libsketchycgal.so")
    filenames = []

    # Add datapath to load if en var is set, used for running tests where shared
    # libraries are built in a different path
    datapath = os.environ.get("DATAPATH", "bazel-bin")
    if datapath is not None:
        f = os.path.join(
            datapath,
            "sketchycgal",
        )
        filenames.append(f)
    # Function to load the library, return True if file system library is loaded

    load_fn = lambda f: ctypes.CDLL(f, mode=ctypes.RTLD_GLOBAL)

    # Try to load all paths for file, fail if none succeed
    for f in filenames:
        sys.path.append(f)

_load_libsketchycgal()
# print(sys.path)
_pywrap_libsketchycgal = importlib.import_module("_libsketchycgal")
