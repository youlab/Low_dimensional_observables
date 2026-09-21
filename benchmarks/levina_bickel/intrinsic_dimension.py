import numpy as np
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt

# Function to compute k-Nearest Neighbors
def kNN(x, n_neighbors, n_jobs):
    # Initialize and fit the NearestNeighbors model to the dataset
    neigh = NearestNeighbors(
        n_neighbors=n_neighbors,  # Number of neighbors to use
        n_jobs=n_jobs             # Number of parallel jobs to run (-1 for all processors)
    ).fit(x)

    # Find the k-nearest neighbors for each point in x
    dists, inds = neigh.kneighbors(x)

    # Return distances and indices of the nearest neighbors
    return dists, inds


# Function to estimate intrinsic dimensionality using the Levina-Bickel MLE method
def levina_bickel(x, dists, k):
    # Compute the log ratio of the distance to the k-th neighbor to the distances of 1st to (k-1)-th neighbors
    m = np.log(dists[:, k:k+1] / dists[:, 1:k])

    # Estimate local intrinsic dimensionality for each point
    m = (k - 1) / m.sum(axis=1)

    # Compute the mean intrinsic dimensionality across all points
    dim = np.mean(m)

    # Return the estimated intrinsic dimension
    return dim


# Main function to fit the model and estimate intrinsic dimensionality
def fit(x, k_list=np.array([20]), n_jobs=4):
    # Ensure k_list is a numpy array, raise error if not
    if not isinstance(k_list, np.ndarray):
        raise TypeError(f"k_list should be np.ndarray but is {type(k_list)}")

    # Determine the maximum number of neighbors needed (+2 to accommodate the algorithm)
    k_max = k_list.max() + 2

    # Compute nearest neighbors distances and indices
    dists, inds = kNN(
        x=x,               # Input dataset
        n_neighbors=k_max, # Max neighbors to compute
        n_jobs=n_jobs      # Number of parallel jobs
    )

    # Estimate intrinsic dimensionality for each k in k_list using Levina-Bickel method
    dimensions = np.array([
        levina_bickel(
            x=x,
            dists=dists,
            k=k
        ) for k in k_list
    ])

    # Return distances, indices, and the array of intrinsic dimension estimates
    return dists, inds, dimensions
