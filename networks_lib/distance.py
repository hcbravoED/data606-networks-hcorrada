import numpy as np
from collections import deque

## TODO: Implement this function
##
## input:
##   mat (np.array): adjacency matrix for graph
## 
## returns:
##   (np.array): distance matrix
##
## Note: You can assume input matrix is binary, square and symmetric 
##       Your output should be square and symmetric
def bfs_distance(mat):
    num_vertices = mat.shape[0]    
    res = np.full((num_vertices, num_vertices), np.inf)
    
    # Finish this loop
    for i in range(num_vertices):
        # keep track of visited vertices in this array
        visited = np.full((num_vertices), False)
    
        # mark the vertex we are starting from as visited
        visited[i] = True
        
        # add this node to the queue
        q = deque([(i,0)])

        # while there are nodes in the queue
        while len(q) > 0:
            # pop vertex and distance from queue
            vertex, d = q.popleft()
            
            # add the distance for this vertex to the result matrix
            res[i, vertex] = d

            # get the list of neighbors for this vertex and process them
            new_vertices = np.nonzero(mat[vertex, :])[0]
            for v in new_vertices:
                if not visited[v]:
                    # if neighbor has not been visited, add to the queue (increasing distance by 1)
                    visited[v] = True
                    q.append((v, d+1))
    return res

## TODO: Implement this function
##
## input:
##   mat (np.array): adjacency matrix for graph
## 
## returns:
##   (list of np.array): list of components
##
## Note: You can assume input matrix is binary, square and symmetric 
##       Your output should be square and symmetric
def get_components(mat):
    dist_mat = bfs_distance(mat)
    
    num_vertices = mat.shape[0]
    available = [True for _ in range(num_vertices)]

    components = []
    
    # finish this loop
    while any(available):
        # grab an available vertex
        cur_vertex = np.arange(num_vertices)[available][0]
        
        # find reachable vertices
        reachables = np.where(np.isfinite(dist_mat[cur_vertex,:]))[0]
        
        # add to the list of components
        components.append(reachables)

        # mark all the nodes reached as not available
        for reachable in reachables:
            available[reachable] = False
        
    return components
