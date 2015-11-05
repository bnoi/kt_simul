import numpy as np


def transformations_matrix(center, vec):
    """Build transformation matrix:
    - translation : from (0, 0) to a point (center)
    - rotation : following angle between (1, 0) and vec

    Parameters
    ----------
    center : list or np.ndarray
    vec : list or np.ndarray

    Returns
    -------
    The transformation matrix, np.ndarray.
    """

    # Setup vectors
    origin_vec = np.array([1, 0])
    current_vec = vec / np.linalg.norm(vec)

    # Find the rotation angle
    a = origin_vec
    b = current_vec
    theta = np.arctan2(a[1], a[0]) + np.arctan2(b[1], b[0])

    # Build rotation matrix
    R = np.array([[np.cos(theta), -np.sin(theta), 0],
                  [np.sin(theta), np.cos(theta), 0],
                  [0, 0, 1]], dtype="float")

    # Build translation matrix
    T = np.array([[1, 0, -center[0]],
                  [0, 1, -center[1]],
                  [0, 0, 1]], dtype="float")

    # Make transformations from R and T in one
    A = np.dot(T.T, R)

    return A
