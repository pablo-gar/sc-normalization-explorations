import numpy as np
from scipy import sparse
import numba
from math import inf 


def _get_mean_var(X, *, axis=0, RSM=False):
    if sparse.issparse(X):
        mean, var = sparse_mean_variance_axis(X, axis=axis, RSM=RSM)
    else:
        mean = np.mean(X, axis=axis, dtype=np.float64)
        
        if RSM:
            mean_of_squares = np.mean(np.square(X), axis=axis, dtype=np.float64)
            var = np.sqrt(mean_of_squares)
        else:
            mean_sq = np.multiply(X, X).mean(axis=axis, dtype=np.float64)
            var = mean_sq - mean ** 2
            
    # enforce R convention (unbiased estimator) for variance
    if not RSM:
        var *= X.shape[axis] / (X.shape[axis] - 1)
        
    return mean, var


def sparse_mean_variance_axis(mtx: sparse.spmatrix, axis: int, RSM: bool):
    """
    This code and internal functions are based on sklearns
    `sparsefuncs.mean_variance_axis`.

    Modifications:
    * allow deciding on the output type, which can increase accuracy when calculating the mean and variance of 32bit floats.
    * This doesn't currently implement support for null values, but could.
    * Uses numba not cython
    """
    assert axis in (0, 1)
    if isinstance(mtx, sparse.csr_matrix):
        ax_minor = 1
        shape = mtx.shape
    elif isinstance(mtx, sparse.csc_matrix):
        ax_minor = 0
        shape = mtx.shape[::-1]
    else:
        raise ValueError("This function only works on sparse csr and csc matrices")
    if axis == ax_minor:
        return sparse_mean_var_major_axis(
            mtx.data, mtx.indices, mtx.indptr, *shape, np.float64, RSM
        )
    else:
        return sparse_mean_var_minor_axis(mtx.data, mtx.indices, *shape, np.float64, RSM)


@numba.njit(cache=True)
def sparse_mean_var_minor_axis(data, indices, major_len, minor_len, dtype, RSM):
    """
    Computes mean and variance for a sparse matrix for the minor axis.

    Given arrays for a csr matrix, returns the means and variances for each
    column back.
    """
    non_zero = indices.shape[0]

    means = np.zeros(minor_len, dtype=dtype)
    lengths = np.zeros_like(means, dtype=dtype)
    variances = np.zeros_like(means, dtype=dtype)
    if RSM:
        mean_of_squares = np.zeros_like(means, dtype=dtype)

    counts = np.zeros(minor_len, dtype=np.int64)

    for i in range(non_zero):
        col_ind = indices[i]
        means[col_ind] += data[i]
        if data[i] > 0 :
            lengths[col_ind] += 1
        if RSM:
            mean_of_squares[col_ind] += data[i]*data[i]
            
    lengths[lengths==0] = 1

    for i in range(minor_len):
        #means[i] /= major_len
        means[i] /= lengths[i]
        if RSM:
            #variances[i] = np.sqrt(mean_of_squares[i] / major_len)
            variances[i] = np.sqrt(mean_of_squares[i] / lengths[i])
    
    # Calculate variance if not using RSM
    if not RSM:
        for i in range(non_zero):
            col_ind = indices[i]
            diff = data[i] - means[col_ind]
            variances[col_ind] += diff * diff
            counts[col_ind] += 1

        for i in range(minor_len):
            #variances[i] += (major_len - counts[i]) * means[i] ** 2
            #variances[i] /= major_len
            variances[i] += (lengths[i] - counts[i]) * means[i] ** 2
            variances[i] /= lengths[i]

    return means, variances


@numba.njit(cache=True)
def sparse_mean_var_major_axis(data, indices, indptr, major_len, minor_len, dtype, RSM):
    """
    Computes mean and variance for a sparse array for the major axis.

    Given arrays for a csr matrix, returns the means and variances for each
    row back.
    """
    means = np.zeros(major_len, dtype=dtype)
    lengths = np.zeros_like(means, dtype=dtype)
    variances = np.zeros_like(means, dtype=dtype)
    if RSM:
        mean_of_squares = np.zeros_like(means, dtype=dtype)

    for i in range(major_len):
        startptr = indptr[i]
        endptr = indptr[i + 1]
        counts = endptr - startptr

        for j in range(startptr, endptr):
            means[i] += data[j]
            if data[j] > 0:
                lengths[i] += 1
            if RSM:
                mean_of_squares[i] += data[j]*data[j]
                
        #means[i] /= minor_len
        lengths[lengths==0] = 1
        means[i] /= lengths[i]
        if RSM:
            #variances[i] = np.sqrt(mean_of_squares[i] / minor_len)
            variances[i] = np.sqrt(mean_of_squares[i] / lengths[i])
        else:
            for j in range(startptr, endptr):
                diff = data[j] - means[i]
                variances[i] += diff * diff

            #variances[i] += (minor_len - counts) * means[i] ** 2
            #variances[i] /= minor_len
            variances[i] += (lengths[i] - counts) * means[i] ** 2
            variances[i] /= lengths[i]

    return means, variances
