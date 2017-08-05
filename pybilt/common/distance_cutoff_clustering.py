import numpy as np

def distance_euclidean(v_a, v_b):
    d_v = v_a - v_b
    return np.sqrt(np.dot(d_v, d_v))

def distance_euclidean_pbc(v_a, v_b, box_lengths, center='zero'):
    if isinstance(center, str):
        if center == 'zero':
            center = np.zeros(len(v_a))
        elif center == 'box_half':
            center = box_lengths/2.0
    #shift center to zero for minimum image
    v_a = v_a - center
    v_b = v_b - center
    #difference
    d_v = v_a - v_b
    d_v_a = np.absolute(d_v)
    dim = len(v_a)
    #check for minimum image
    for i in range(dim):
        v_i = d_v_a[i]
        box_i = box_lengths[i]
        box_i_h = box_i/2.0
        if v_i > box_i_h:
            d_v[i] = box_i - np.absolute(v_a[i]) - np.absolute(v_b[i])
    r = np.sqrt(np.dot(d_v, d_v))
    return r



def distance_cutoff_clustering(vectors, cutoff, dist_func, min_size=1, *df_args, **df_kwargs):
    #compute the boolean distance-cutoff matrix
    #print(df_args)
    #print(df_kwargs)
    nvecs = len(vectors)
    dist_bool = np.zeros((nvecs,nvecs), dtype=np.bool)
    for i in range(nvecs-1):
        vec_a = vectors[i]
        for j in range(i+1, nvecs):
            vec_b = vectors[j]
            dist = dist_func(vec_a, vec_b, *df_args, **df_kwargs)
            if dist <= cutoff:
                dist_bool[i][j] = True
                dist_bool[j][i] = True
            else:
                dist_bool[i][j] = False
                dist_bool[j][i] = False
    clusters = []
    master = []

    for i in range(nvecs):
        master.append([i, False])
    #print(master)
    #clustind = 0
    neighbors = []

    while len(master)>0:
        start = master[0][0]
        master[0][1] = True
        #print(master)
        #rest the neigbor list
        neighbors = []
        #seed neighborlist with start
        neighbors.append(start)
        #now loop over the neigbor list and build neigbors and neighbors of neighbors
        i = 0
        while i < len(neighbors):
            #print(neighbors)
            a = neighbors[i]
            #vec_a = vectors[a]
            for j in range(len(master)):
                b = master[j][0]
                if not master[j][1]:
                    #print "a ",a," b ",b
                    #vec_b = vectors[b]
                    if dist_bool[a][b]:
                        neighbors.append(b)
                        master[j][1] = True
                        #print(neighbors)
            #print(master)
            i+=1
        master = list([v for v in master if not v[1]])
        #print(master)
        if len(neighbors) > min_size:
            clusters.append(list(neighbors))
            #clusters[clustind] = list(neighbors)
    return clusters
