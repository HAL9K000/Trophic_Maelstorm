"""
GPU Acceleration Module with Automatic CPU Fallback
Provides drop-in replacements for numpy/scipy with automatic GPU acceleration where possible
"""

import os
import sys
import warnings
from functools import wraps
import importlib

# Determine GPU usage from environment variable
USE_CUDA = os.getenv("USE_GPU", "0") == "1"
GPU_AVAILABLE = False

# Standard CPU imports (always available)
import numpy as _cpu_np
import scipy as _cpu_scipy
import scipy.fft as _cpu_fft
import scipy.signal as _cpu_signal
import scipy.interpolate as _cpu_interpolate
import scipy.stats as _cpu_stats

# Try GPU imports
if USE_CUDA:
    try:
        if sys.platform.startswith("darwin"):
            raise ImportError("CUDA libraries are not supported on macOS, use Linux machine for GPU support.")

        import cupy as _cupy
        import cupyx.scipy.fft as _gpu_fft
        import cupyx.scipy.signal as _gpu_signal

        GPU_AVAILABLE = True
        # RAPIDs libraries (WORK ONLY ON LINUX WITH NVIDIA GPUS)
        if sys.platform.startswith("linux"):
            try:
                import cudf.pandas
                cudf.pandas.install()
                import cuml.accel
                cuml.accel.install()
                print("RAPIDs acceleration enabled")
            except ImportError as e:
                warnings.warn(f"RAPIDs libraries not found: {e}. Using CuPy only for GPU acceleration...")
        elif sys.platform.startswith("win"):
            warnings.warn("RAPIDs libraries are not supported on Windows, using CuPy only for GPU support.")
    except ImportError as e:
        warnings.warn(f"Error importing GPU libraries: {e}. Falling back to CPU usage.")
        USE_CUDA = False
        GPU_AVAILABLE = False

# Utlity functions and wrapper modules
"""Recursively convert GPU arrays to CPU arrays."""
def to_cpu(obj):
    
    if not GPU_AVAILABLE:
        return obj
    if hasattr(_cupy, 'ndarray') and isinstance(obj, _cupy.ndarray):
        return obj.get()
    elif isinstance(obj, dict):
        return {k: to_cpu(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return type(obj)(to_cpu(v) for v in obj)
    elif isinstance(obj, set):
        return {to_cpu(v) for v in obj}
    return obj

"""Convert CPU arrays (iterators) to GPU arrays where possible."""
def to_gpu(obj):
    
    if not GPU_AVAILABLE:
        return obj
    if hasattr(_cupy, 'ndarray') and isinstance(obj, _cupy.ndarray):
        return obj
    if isinstance(obj, _cpu_np.ndarray):
        try:
            return _cupy.asarray(obj)
        except Exception:
            return obj
    elif isinstance(obj, dict):
        return {k: to_gpu(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return type(obj)(to_gpu(v) for v in obj)
    elif isinstance(obj, set):
        return {to_gpu(v) for v in obj}
    return obj

def cpu_fallback(func_name, cpu_module, gpu_module=None):
    """Decorator to provide CPU fallback for GPU functions."""
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            if GPU_AVAILABLE and gpu_module:
                try:
                    # Try GPU version first
                    gpu_args = to_gpu(args)
                    gpu_kwargs = to_gpu(kwargs)
                    result = func(*gpu_args, **gpu_kwargs)
                    return to_cpu(result)
                except Exception as e:
                    warnings.warn(f"GPU operation failed for {func_name}: {e}. Falling back to CPU.")
            
            # Use CPU version
            cpu_args = to_cpu(args)
            cpu_kwargs = to_cpu(kwargs)
            cpu_func = getattr(cpu_module, func_name)
            return cpu_func(*cpu_args, **cpu_kwargs)
        return wrapper
    return decorator

class SmartModule:
    """A module wrapper that tries GPU first, falls back to CPU."""
    def __init__(self, cpu_module, gpu_module=None, module_name=""):
        self._cpu_module = cpu_module
        self._gpu_module = gpu_module
        self._module_name = module_name

        # Define attributes that should always come from CPU module
        # (types, constants, etc.)
        self._cpu_only_attrs = {
            'ndarray', 'dtype', 'generic',
            'float64', 'float32', 'float16', 'int64', 'int32', 'int16', 'int8',
            'uint64', 'uint32', 'uint16', 'uint8', 'bool_', 'complex64', 'complex128',
            #Constants
            'nan', 'inf', 'pi', 'e', 'euler_gamma', 'newaxis', 'NINF', 'NZERO', 'PZERO' }
            # Add other constants and types as needed
        
        
    def __getattr__(self, name):
        # For types and constants, always use CPU module
        if name in self._cpu_only_attrs:
            return getattr(self._cpu_module, name)
        # Check if attribute exists in GPU module first
        if GPU_AVAILABLE and self._gpu_module and hasattr(self._gpu_module, name):
            gpu_attr = getattr(self._gpu_module, name)
            if callable(gpu_attr):
                @wraps(gpu_attr)
                def smart_func(*args, **kwargs):
                    try:
                        # Convert args to GPU
                        gpu_args = to_gpu(args)
                        gpu_kwargs = to_gpu(kwargs)
                        result = gpu_attr(*gpu_args, **gpu_kwargs)
                        # ✅ Keep result on GPU! Don't convert to CPU automatically
                        return result  # or to_cpu(result) if you want to ensure always returning
                        #CPU arrays for compatibility (MASSIVE performance hit)
                    except Exception as e:
                        # Fall back to CPU
                        warnings.warn(f"GPU operation failed for {self._module_name}.{name}: {e}")
                        cpu_attr = getattr(self._cpu_module, name)
                        cpu_args = to_cpu(args)
                        cpu_kwargs = to_cpu(kwargs)
                        return cpu_attr(*cpu_args, **cpu_kwargs)
                return smart_func
            else:
                return gpu_attr
        
        # Fall back to CPU module
        cpu_attr = getattr(self._cpu_module, name)
        if callable(cpu_attr):
            @wraps(cpu_attr)
            def cpu_func(*args, **kwargs):
                cpu_args = to_cpu(args)
                cpu_kwargs = to_cpu(kwargs)
                return cpu_attr(*cpu_args, **cpu_kwargs)
            return cpu_func
        return cpu_attr

class CPUOnlyModule:
    """Wrapper for modules that don't have GPU equivalents."""
    def __init__(self, cpu_module, module_name=""):
        self._cpu_module = cpu_module
        self._module_name = module_name
        
    def __getattr__(self, name):
        attr = getattr(self._cpu_module, name)
        if callable(attr):
            @wraps(attr)
            def cpu_only_func(*args, **kwargs):
                # Convert all inputs to CPU
                cpu_args = to_cpu(args)
                cpu_kwargs = to_cpu(kwargs)
                return attr(*cpu_args, **cpu_kwargs)
            return cpu_only_func
        return attr

""" IMPORTANT: Decorator to ensure function ALWAYS returns CPU arrays - USE SPARINGLY DUE TO PERFORMANCE HIT! """     
def ensure_cpu(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        result = func(*args, **kwargs)
        return to_cpu(result)
    return wrapper

def cpu_safe_import(module_path):
    """Import and wrap a CPU-only module."""
    try:
        mod = importlib.import_module(module_path)
        return CPUOnlyModule(mod, module_path)
    except ImportError as e:
        raise ImportError(f"Could not import {module_path}: {e}")

# Create smart modules
if GPU_AVAILABLE:
    numpy = SmartModule(_cpu_np, _cupy, "numpy")
    fft = SmartModule(_cpu_fft, _gpu_fft, "scipy.fft")
    signal = SmartModule(_cpu_signal, _gpu_signal, "scipy.signal")
    # Pre-wrapped common CPU-only modules for convenience
    interpolate = cpu_safe_import("scipy.interpolate")
    stats = cpu_safe_import("scipy.stats")
else:
    numpy = SmartModule(_cpu_np, None, "numpy")
    fft = SmartModule(_cpu_fft, None, "scipy.fft")
    signal = SmartModule(_cpu_signal, None, "scipy.signal")
    interpolate = _cpu_scipy.interpolate
    stats = _cpu_scipy.stats

# Create a numpy alias for convenience
np = numpy


    
# OTHER USEFUL UTILITIES FOR EXPLICIT TO_CPU() CONVERSIONS BY USER!
def asnumpy(arr):
    """Explicitly convert array to CPU numpy array - use when you need CPU data"""
    return to_cpu(arr)

def asarray(arr):
    """Convert to appropriate array type (GPU if available, CPU otherwise)"""
    if GPU_AVAILABLE:
        return to_gpu(arr)
    return _cpu_np.asarray(arr)

def is_gpu_array(arr):
    """Check if array is on GPU"""
    if not GPU_AVAILABLE:
        return False
    return hasattr(_cupy, 'ndarray') and isinstance(arr, _cupy.ndarray)



# Expose key functions at module level
__all__ = ['np', 'numpy', 'fft', 'signal', 'interpolate', 'stats', 
           'to_cpu', 'to_gpu', 'asnumpy', 'asarray', 'is_gpu_array',
           'cpu_safe_import', 'ensure_cpu', 'GPU_AVAILABLE', 'USE_GPU']



''' OLD SCRIPT

if USE_CUDA:
    try:

        if sys.platform.startswith("darwin"):
            raise ImportError("CUDA libraries are not supported on macOS, use Linux machine for GPU support.")
        import cupy as np
        import cupyx.scipy.fft as fft
        import cupyx.scipy.signal as signal
        # RAPIDs libraries (WORK ONLY ON LINUX WITH NVIDIA GPUS)
        if sys.platform.startswith("linux"):
            import cudf.pandas
            cudf.pandas.install()
            import cuml.accel
            cuml.accel.install()
        elif sys.platform.startswith("win"):
            print("RAPIDs libraries are not supported on WIN, using cupy only.")
    except ImportError as e:
        print(f"Error importing GPU libraries: {e}")
        print("Falling back to CPU usage.")
        USE_CUDA = False
        import numpy as np
        import scipy.fft as fft
        import scipy.signal as signal
else:
    import numpy as np
    import scipy.fft as fft
    import scipy.signal as signal

# Recursively convert CuPy arrays to NumPy arrays
def to_cpu(obj):
    """Recursively convert CuPy arrays to NumPy arrays."""
    if USE_CUDA:
        import cupy
        if isinstance(obj, cupy.ndarray):
            return obj.get()
        elif isinstance(obj, dict):
            return {k: to_cpu(v) for k, v in obj.items()}
        elif isinstance(obj, (list, tuple, set)):
            return type(obj)(to_cpu(v) for v in obj)
        elif hasattr(obj, '__array__') and not isinstance(obj, _np.ndarray):
            return _np.asarray(obj)
    return obj

"""Recursively convert CuPy arrays to NumPy arrays."""
def to_cpu(obj):
    if USE_CUDA:
        import cupy
        if isinstance(obj, cupy.ndarray):
            return obj.get()  # returns a NumPy array
        elif isinstance(obj, dict):
            return {k: to_cpu(v) for k, v in obj.items()}
        elif isinstance(obj, (list, tuple, set)):
            return type(obj)(to_cpu(v) for v in obj)
        elif hasattr(obj, '__array__') and not isinstance(obj, _np.ndarray):
            return _np.asarray(obj)  # ensure it's a NumPy array
    return obj
'''
