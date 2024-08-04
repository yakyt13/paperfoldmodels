import numpy as np
import openmesh as om
import networkx as nx
from tabuunfold import *
from usefullfunctions import *
import multiprocessing
from functools import partial


def find_best_spanning_tree(mesh, dual_graph, num_iterations=100):
    best_tree = None
    best_unfolding = None
    min_overlaps = float('inf')

    for _ in range(num_iterations):
        # Generate a random spanning tree
        spanning_tree = randomSpanningTree(dual_graph)

        # Unfold the mesh using the spanning tree
        unfolded_mesh, _, _, _, _ = unfoldSpanningTree(mesh, spanning_tree)

        # Count the number of overlaps in the unfolded mesh
        num_overlaps, _ = detectOverlaps(unfolded_mesh, spanning_tree)

        # Update the best tree if the current one has fewer overlaps
        if num_overlaps < min_overlaps:
            min_overlaps = num_overlaps
            best_tree = spanning_tree
            best_unfolding = unfolded_mesh

    return best_tree, best_unfolding, min_overlaps


def unfold(mesh, count=1):
    # Calculate the number of faces, edges and corners, as well as the lengths of the longest shortest edge
    # To be deleted
    numEdges = mesh.n_edges()
    numVertices = mesh.n_vertices()
    numFaces = mesh.n_faces()

    # Create the dual graph of the network and calculate the weights
    dualGraph = nx.Graph()

    # For the weights: Calculate the longest and shortest edges of the triangle
    minLength = 1000
    maxLength = 0
    for edge in mesh.edges():
        edgelength = mesh.calc_edge_length(edge)
        if edgelength < minLength:
            minLength = edgelength
        if edgelength > maxLength:
            maxLength = edgelength
    # minLength is length of shortest edge of model
    # maxLength is length of longest edge of model

    edges_list = list()

    # All edges in the network
    for edge in mesh.edges():
        # The two sides that border the edge
        # In a triangular mesh, each edge consists of two halfedges.
        # halfedge associated with one of the two faces that share the edge
        # These halfedges have opposite directions.
        face1 = mesh.face_handle(mesh.halfedge_handle(edge, 0))     # Get the first halfedge
        face2 = mesh.face_handle(mesh.halfedge_handle(edge, 1))     # Get the second halfedge

        # The weight
        # We normalize edge length and then invert the normalized value (1.0 - edge normalized weight)
        # So the shortest edges gets a weight close to 1 and the longest edge gets a weight close to 0
        # This type of normalization can be useful in algorithms where edge weights play a role in decision-making
        edgeweight = 1.0 - (mesh.calc_edge_length(edge) - minLength) / (maxLength - minLength)

        # Calculate the center points of the sides (only necessary for visualization)
        center1 = (0, 0)
        for vertex in mesh.fv(face1):
            center1 = center1 + 0.3333333333333333 * np.array([mesh.point(vertex)[0], mesh.point(vertex)[2]])
        center2 = (0, 0)
        for vertex in mesh.fv(face2):
            center2 = center2 + 0.3333333333333333 * np.array([mesh.point(vertex)[0], mesh.point(vertex)[2]])

        # Add the new nodes and edges to the dual graph
        # Node in dual graph represents edge in the original model
        # edge in dual graph represents edge between two neighboring faces of original model
        dualGraph.add_node(face1.idx(), pos=center1)        # First face
        dualGraph.add_node(face2.idx(), pos=center2)        # Second face
        dualGraph.add_edge(face1.idx(), face2.idx(), idx=edge.idx(), weight=edgeweight)     # Edge between two faces
        e = list()
        e.append(face1.idx())
        e.append(face2.idx())
        e.append(edge.idx())
        edges_list.append(e)

    # print(f'processing tree #{count}')
    # spanningTree, overlaps = tabuUnfolding(mesh, randomSpanningTree(dualGraph), edges_list)
    # while overlaps:
    #     print(f'No solution found. Remaining overlaps: {overlaps}')
    #     print('')
    #     count += 1
    #     print(f'processing tree #{count}')
    #     spanningTree, overlaps = tabuUnfolding(mesh, randomSpanningTree(dualGraph), edges_list)

    spanningTree, _ = tabuUnfolding(mesh, randomSpanningTree(dualGraph), edges_list)


    # num_parallel_runs = 4
    # with multiprocessing.Pool(num_parallel_runs) as pool:
    #     result = pool.map(partial(tabuUnfolding, mesh, randomSpanningTree(dualGraph)), edges_list)
    # for overlaps, tree in result:
    #     if overlaps == 0:
    #         spanningTree = tree
    #         break

    # spanningTree = nx.minimum_spanning_tree(dualGraph)

    fullUnfolding = unfoldSpanningTree(mesh, spanningTree)
    [unfoldedMesh, isFoldingEdge, connections, glueNumber, foldingDirection] = fullUnfolding

    # # Resolve the intersections
    # # Find all intersections
    # epsilon = 1E-12  # accuracy
    # faceIntersections = []
    # for face1 in unfoldedMesh.faces():
    #     for face2 in unfoldedMesh.faces():
    #         if face2.idx() < face1.idx():  # So that we don't go through the pairs twice
    #             # Get the triangle faces
    #             triangle1 = []
    #             triangle2 = []
    #             for halfedge in unfoldedMesh.fh(face1):
    #                 triangle1.append(unfoldedMesh.point(unfoldedMesh.from_vertex_handle(halfedge)))
    #             for halfedge in unfoldedMesh.fh(face2):
    #                 triangle2.append(unfoldedMesh.point(unfoldedMesh.from_vertex_handle(halfedge)))
    #             if triangleIntersection(triangle1, triangle2, epsilon):
    #                 faceIntersections.append([connections[face1.idx()], connections[face2.idx()]])
    #
    #
    # # Find the paths
    # # We find the minimum number of cuts to resolve each self-intersection
    #
    # # Find all paths between intersecting triangles
    # paths = []
    # for intersection in faceIntersections:
    #     paths.append(
    #         nx.algorithms.shortest_paths.shortest_path(spanningTree, source=intersection[0], target=intersection[1]))
    #
    # # Find all edges in all paths
    # edgepaths = []
    # for path in paths:
    #     edgepath = []
    #     for i in range(len(path) - 1):
    #         edgepath.append((path[i], path[i + 1]))
    #     edgepaths.append(edgepath)
    #
    # # List of all edges in all paths
    # allEdgesInPaths = list(set().union(*edgepaths))
    #
    # # Count how often each edge appears
    # numEdgesInPaths = []
    # for edge in allEdgesInPaths:
    #     num = 0
    #     for path in edgepaths:
    #         if edge in path:
    #             num = num + 1
    #     numEdgesInPaths.append(num)
    #
    # S = []
    # C = []
    #
    # while len(C) != len(paths):
    #     # Calculate the weights to decide which edge we cut
    #     cutWeights = np.empty(len(allEdgesInPaths))
    #     for i in range(len(allEdgesInPaths)):
    #         currentEdge = allEdgesInPaths[i]
    #
    #         # Count how many of the paths in which the edge occurs have already been cut
    #         numInC = 0
    #         for path in C:
    #             if currentEdge in path:
    #                 numInC = numInC + 1
    #
    #         # Determine the weight
    #         if (numEdgesInPaths[i] - numInC) > 0:
    #             cutWeights[i] = 1 / (numEdgesInPaths[i] - numInC)
    #         else:
    #             cutWeights[i] = 1000  # 1000 = infinite
    #     # Find the edge with the smallest weight
    #     minimalIndex = np.argmin(cutWeights)
    #     S.append(allEdgesInPaths[minimalIndex])
    #     # Find all paths where the edge occurs and add them to C
    #     for path in edgepaths:
    #         if allEdgesInPaths[minimalIndex] in path and not path in C:
    #             C.append(path)
    #
    # # Now we remove the cut edges from the minimal spanning tree
    # spanningTree.remove_edges_from(S)

    # Find the related components
    # connectedComponents = nx.algorithms.components.connected_components(new_tree)
    connectedComponents = nx.algorithms.components.connected_components(spanningTree)
    connectedComponentList = list(connectedComponents)


    # Processing of the components
    unfoldings = []
    for component in connectedComponentList:
        unfoldings.append(unfoldSpanningTree(mesh, spanningTree.subgraph(component)))

    return fullUnfolding, unfoldings

