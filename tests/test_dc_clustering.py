from __future__ import print_function
import numpy as np
from pybilt.common import distance_cutoff_clustering as dc_cluster

def test_dc_clustering():
    #Euclidean distance
    vectors = np.array([[0.1], [0.9], [3.0], [3.7]])
    dist_func = dc_cluster.distance_euclidean
    cutoff = 1.0
    cluster_ref = [[0, 1], [2, 3]]
    clusters = dc_cluster.distance_cutoff_clustering(vectors, cutoff, dist_func)
    print("Euclidean distance cutoff clustering:")
    print((" output matches the reference: {}".format(cluster_ref==clusters)))
    #Euclidean distance with periodic boundaries
    box = np.array([4.0])
    center = np.array([0.0])
    vectors = np.array([[-1.8],[-0.1],[0.1],[1.8]])
    dist_func = dc_cluster.distance_euclidean_pbc
    cluster_ref = [[0, 3], [1, 2]]
    clusters = dc_cluster.distance_cutoff_clustering(vectors, cutoff, dist_func, 1, box, center=center)
    print("Euclidean distance (with periodic boundaries) cutoff clustering:")
    print((" output matches the reference: {}".format(cluster_ref==clusters)))

    return

if __name__ == '__main__':
    test_dc_clustering()
