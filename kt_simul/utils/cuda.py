import numpy as np

cuda_linalg_solver = False

try:
    import pycuda.autoinit
    import pycuda.gpuarray as gpuarray
    from scikits.cuda import linalg
    linalg.init()
    cuda_linalg_solver = True
except:
    pass


def linalg_solve_cuda(a, b):
    """Note that Cholesky factorization requires a matrix that is symmetric and positive definite.
    """

    if not cuda_linalg_solver:
        return np.linalg.solve(a, b)

    a_gpu = gpuarray.to_gpu(a)
    b_gpu = gpuarray.to_gpu(b)

    linalg.cho_solve(a_gpu, b_gpu)
    return b_gpu.get()
