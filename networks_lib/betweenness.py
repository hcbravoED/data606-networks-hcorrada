import numpy as np
from collections import deque

# first pass of the newman-girvan algorithm 
# path counting algorithm starting from vertex 
def _eb_down(starting_vertex, mat):
    num_vertices = mat.shape[0]

    # store the number of paths to each vertex
    num_paths = np.full((num_vertices), -1, dtype=np.int)
    
    # store the search level for each vector
    # initialized to infinity (we use this as 'not visited')
    search_level = np.full((num_vertices), np.inf)

    # the search level for starting vertex is 0
    search_level[starting_vertex] = 0
    
    # initialize the bfs queue
    # with vertex id and level
    q = deque([(starting_vertex, 0)])

    # while the queue is not empty
    while len(q) > 0:
        # pop vertex and search level from queue
        vertex, level = q.popleft()
        
        # get the neighbors (both parents and children)
        neighbors = np.nonzero(mat[vertex, :])[0]

        # figure out which are parents (they are in the previous level)
        parents = [neighbor for neighbor in neighbors if search_level[neighbor] == level - 1]
        
        # sum the number of paths for the parents
        num_paths[vertex] = np.sum(num_paths[parents]) if len(parents) > 0 else 1

        # figure out which children to add to the queue (hasn't been visited yet)
        children = [neighbor for neighbor in neighbors if np.isinf(search_level[neighbor])]
        for child in children:
            # add to the queue (setting the level marks as visited)
            search_level[child] = level + 1
            q.append((child, level+1))

    # return the number of paths calculation and the search level for each vertex
    return num_paths, search_level

# second pass of the newman-girvan algorithm
# path counting algorithm starting from a given vertex
def _eb_up(starting_vertex, mat, num_paths, search_level):
    num_vertices = mat.shape[0]
    
    # store vertex credits here
    vertex_credit = np.full((num_vertices), np.inf)
    
    # store edge credits (final result) here
    edge_credit = np.zeros(mat.shape)

    # start at the lowest level
    cur_level = np.amax(search_level)
    
    # go through each level
    while cur_level >= 0:
        # get the vertices in the current level
        cur_vertices = np.where(search_level == cur_level)[0]

        for vertex in cur_vertices:
            # get its neighbors
            neighbors = np.nonzero(mat[vertex, :])[0]

            # deal with children first
            children = [neighbor for neighbor in neighbors if search_level[neighbor] == cur_level + 1]
            
            # set vertex credit to 1 + the sum of the edge credits incoming from children
            vertex_credit[vertex] = 1 + (np.sum(edge_credit[vertex, children]) if len(children) > 0 else 0)

            # figure out which neighbors are parents
            parents = [neighbor for neighbor in neighbors if search_level[neighbor] == cur_level - 1]
            if len(parents) == 0:
                continue

            # this is the total credit to parcel out
            credit_factor = 1.0 * vertex_credit[vertex] / np.sum(num_paths[parents])
            
            # divide among edges outgoing to each parent
            for parent in parents:
                # set the credit of outgoing edge
                edge_credit[vertex,parent] = credit_factor * num_paths[parent]
                edge_credit[parent,vertex] = edge_credit[vertex,parent]
                
        # go to the next level
        cur_level = cur_level - 1
        
    return edge_credit

## TODO: Implement this function
##
## Implements the breadth-first algorithm of Girvan-Newman to compute
##   number (fractional) of shortest paths starting from a given vertex
##   that go through each edge of the graph
##
## Input:
##   - vertex (int): index of vertex paths start from
##   - mat (np.array): n-by-n adjacency matrix
##
## Output:
##   (np.array): n-by-n edge count matrix
##
## Note: assume input adjacency matrix is binary and symmetric
def edge_counts(vertex, mat):
    num_vertices = mat.shape[0]
    num_paths, search_level = _eb_down(vertex, mat)
    return _eb_up(vertex, mat, num_paths, search_level)

## Compute edge betweeness for a graph
## 
## Input: 
##   - mat (np.array): n-by-n adjacency matrix. 
##
## Output:
##   (np.array): n-by-n matrix of edge betweenness
##
## Notes: Input matrix is assumed binary and symmetric
def edge_betweenness(mat):
    res = np.zeros(mat.shape)
    num_vertices = mat.shape[0]
    for i in range(num_vertices):
        res += edge_counts(i, mat)
    return res / 2.