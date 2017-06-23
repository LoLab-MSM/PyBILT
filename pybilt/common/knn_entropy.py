import numpy as np
from scipy.special import gamma,psi
from scipy.spatial.distance import cdist

# Get the k nearest neighbors
# between points in a random variable/vector.
# Uses Euclidean style distances.
def k_nearest_neighbors(X, k=1):
    #length
    nX = len(X)
    #initialize knn dict
    knn = {key: [] for key in xrange(nX)}
    #make sure X has the right shape for the cdist function
    X = np.reshape(X, (nX,-1))
    dists_arr = cdist(X, X)
    distances = [[i,j,dists_arr[i,j]] for i in xrange(nX-1) for j in xrange(i+1,nX)]
    #sort distances
    distances.sort(key=lambda x: x[2])
    #pick up the k nearest
    for d in distances:
        i = d[0]
        j = d[1]
        dist = d[2]
        if len(knn[i]) < k:
            knn[i].append([j, dist])
        if len(knn[j]) < k:
            knn[j].append([i, dist])
    return knn


# knn kth neighbor distances for entropy calcs.
def kth_nearest_neighbor_distances(X, k=1):

    #length
    nX = len(X)
    #make sure X has the right shape for the cdist function
    X = np.reshape(X, (nX,-1))
    dists_arr = cdist(X, X)
    #sorts each row
    dists_arr.sort()
    return [dists_arr[i][k] for i in xrange(nX)]

def shannon_entropy(X, k=1, kth_dists=None):
    # KL entropy estimator from
    #  https://arxiv.org/pdf/1506.06501v1.pdf
    # Also see:
    # https://www.cs.tut.fi/~timhome/tim/tim/core/differential_entropy_kl_details.htm
    # and
    # Kozachenko, L. F. & Leonenko, N. N. 1987 Sample estimate of entropy
    #    of a random vector. Probl. Inf. Transm. 23, 95-101.
    # the kth nearest neighbor distances
    r_k = kth_dists
    if kth_dists is None:
        r_k = kth_nearest_neighbor_distances(X, k=k)
    #length
    n = len(X)
    #dimension
    d = 1
    if len(X.shape) == 2:
        d = X.shape[1]
    # volume of unit ball in d^n
    v_unit_ball = np.pi**(0.5*d)/gamma(0.5*d + 1.0)
    # log distances
    lr_k = np.log(r_k)
    # Shannon entropy estimate  
    H = psi(n) - psi(k) + np.log(v_unit_ball) + (np.float(d)/np.float(n))*( lr_k.sum())

    return H

# Not entirely sure this one is correct.
def shannon_entropy_pc(X, k=1, kth_dists=None):
    # entropy estimator from
    # F. Perez-Cruz, (2008). Estimation of Information Theoretic Measures
    #     for Continuous Random Variables. Advances in Neural Information
    #     Processing Systems 21 (NIPS). Vancouver (Canada), December.
    #    https://papers.nips.cc/paper/3417-estimation-of-information-theoretic-measures-for-continuous-random-variables.pdf
    # 
    #
    r_k = kth_dists
    if kth_dists is None:
        r_k = np.array(kth_nearest_neighbor_distances(X, k=k))
    n = len(X)
    d = 1
    if len(X.shape) == 2:
        d = X.shape[1]
    #volume of the unit ball
    v_unit_ball = np.pi**(0.5*d)/gamma(0.5*d + 1.0)
    # probability estimator using knn distances 
    p_k_hat = (k / (n -1.0)) * (1.0/v_unit_ball) * (1.0/r_k**d)
    # log probability
    log_p_k_hat = np.log(p_k_hat)    
    #entropy estimator
    h_k_hat = log_p_k_hat.sum() / (-1.0*n)
    return h_k_hat
 
#mutual information
def mutual_information(var_tuple, k=2):
    nvar = len(var_tuple)
    #make sure the input arrays are properly shaped for hstacking
    var_tuple = tuple( var_tuple[i].reshape(len(var_tuple[i]),-1) for i in xrange(nvar) )
    #compute the individual entropies of each variable  
    Hx = [shannon_entropy(var_tuple[i],k=k) for i in xrange(nvar)]
    Hx = np.array(Hx)
    # and get the sum
    Hxtot = Hx.sum()
    #now get the entropy of the joint distribution
    joint = np.hstack(var_tuple)
    Hjoint = shannon_entropy(joint, k=k)
    #get the mutual information
    MI = Hxtot - Hjoint
    #set to zero if value is negative
    if MI < 0.0: MI = 0.0
    #return
    return MI 
  	
# conditional mutual information
def conditional_mutual_information(var_tuple, cond_tuple, k=2):
    nvar = len(var_tuple)
    ncon = len(cond_tuple)
    #make sure the input arrays are properly shaped for hstacking
    var_tuple = tuple( var_tuple[i].reshape(len(var_tuple[i]),-1) for i in xrange(nvar) )
    cond_tuple = tuple( cond_tuple[i].reshape(len(cond_tuple[i]),-1) for i in xrange(ncon) )
    # compute pair joint entropies
    Hxz = [shannon_entropy(np.hstack(var_tuple[i]+cond_tuple), k=k) for i in xrange(nvar)]
    Hxz = np.array(Hxz)
    jtup = var_tuple + cond_tuple
    joint = np.hstack( jtup )
    Hj = shannon_entropy(joint, k=k)
    Hz = 0.0
    if len(cond_tuple) > 1:
        joint = np.hstack(cond_tuple)
        Hz = shannon_entropy(joint, k=k)
    else:
        Hz = shannon_entropy(cond_tuple[0], k=k)
    Hxzsum = Hxz.sum()
   # print "Hxzsum: ",Hxzsum, " Hj: ",Hj, " Hz: ",Hz
    MIc = Hxzsum - Hj - Hz
   # print "MIc: ",MIc
    if MIc < 0.0: MIc = 0.0
    #print MIc
    return MIc
