import numpy as np
import openmesh as om
import networkx as nx
from tabuunfold import *
from usefullfunctions import *
from usefullfunctions import detectOverlaps


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


# Function that unwinds a spanning tree
def unfoldSpanningTree(mesh, spanningTree):
    unfoldedMesh = om.TriMesh()  # The unwound network

    numFaces = mesh.n_faces()
    sizeTree = spanningTree.number_of_edges()
    numUnfoldedEdges = 3 * numFaces - sizeTree

    isFoldingEdge = np.zeros(numUnfoldedEdges, dtype=bool)  # Specifies whether an edge is folded or cut
    glueNumber = np.empty(numUnfoldedEdges, dtype=int)  # Saves which edge is glued together
    foldingDirection = np.empty(numUnfoldedEdges, dtype=int)  # Valley folding or mountain folding

    connections = np.empty(numFaces, dtype=int)  # Stores which original triangle belongs to the unfolded one

    # Choose the first triangle arbitrarily
    startingNode = list(spanningTree.nodes())[0]
    startingTriangle = mesh.face_handle(startingNode)

    # We unwind the first triangle

    # All half edges of the first triangle
    firstHalfEdge = mesh.halfedge_handle(startingTriangle)
    secondHalfEdge = mesh.next_halfedge_handle(firstHalfEdge)
    thirdHalfEdge = mesh.next_halfedge_handle(secondHalfEdge)
    originalHalfEdges = [firstHalfEdge, secondHalfEdge, thirdHalfEdge]

    # Calculate the lengths of the edges, this determines the shape of the triangle (congruence)
    edgelengths = [mesh.calc_edge_length(firstHalfEdge), mesh.calc_edge_length(secondHalfEdge),
                   mesh.calc_edge_length(thirdHalfEdge)]

    # The first two points
    firstUnfoldedPoint = np.array([0, 0, 0])
    secondUnfoldedPoint = np.array([edgelengths[0], 0, 0])

    # We calculate the third point of the triangle from the first two. There are two possibilities
    [thirdUnfolded0, thirdUnfolded1] = getThirdPoint(firstUnfoldedPoint, secondUnfoldedPoint, edgelengths[0],
                                                     edgelengths[1],
                                                     edgelengths[2])
    if thirdUnfolded0[1] > 0:
        thirdUnfoldedPoint = thirdUnfolded0
    else:
        thirdUnfoldedPoint = thirdUnfolded1

    # Add the new corners to the unfolded mesh
    # firstUnfoldedVertex = unfoldedMesh.add_vertex(secondUnfoldedPoint)
    # secondUnfoldedVertex = unfoldedMesh.add_vertex(thirdUnfoldedPoint)
    # thirdUnfoldedVertex = unfoldedMesh.add_vertex(firstUnfoldedPoint)

    firstUnfoldedVertex = unfoldedMesh.add_vertex(firstUnfoldedPoint)
    secondUnfoldedVertex = unfoldedMesh.add_vertex(secondUnfoldedPoint)
    thirdUnfoldedVertex = unfoldedMesh.add_vertex(thirdUnfoldedPoint)

    # Create the face
    unfoldedFace = unfoldedMesh.add_face(firstUnfoldedVertex, secondUnfoldedVertex, thirdUnfoldedVertex)

    # Save properties of the face and edges
    # The half edges in the unrolled grid
    firstUnfoldedHalfEdge = unfoldedMesh.next_halfedge_handle(unfoldedMesh.opposite_halfedge_handle(unfoldedMesh.halfedge_handle(firstUnfoldedVertex)))
    secondUnfoldedHalfEdge = unfoldedMesh.next_halfedge_handle(firstUnfoldedHalfEdge)
    thirdUnfoldedHalfEdge = unfoldedMesh.next_halfedge_handle(secondUnfoldedHalfEdge)

    unfoldedHalfEdges = [firstUnfoldedHalfEdge, secondUnfoldedHalfEdge, thirdUnfoldedHalfEdge]

    # Associated triangle in the 3D mesh
    connections[unfoldedFace.idx()] = startingTriangle.idx()
    # Folding direction and glue number
    addVisualisationData(mesh, unfoldedMesh, originalHalfEdges, unfoldedHalfEdges, glueNumber, foldingDirection)

    halfEdgeConnections = {firstHalfEdge.idx(): firstUnfoldedHalfEdge.idx(),
                           secondHalfEdge.idx(): secondUnfoldedHalfEdge.idx(),
                           thirdHalfEdge.idx(): thirdUnfoldedHalfEdge.idx()}

    # We go through the tree
    for dualEdge in nx.dfs_edges(spanningTree, source=startingNode):
        foldingEdge = mesh.edge_handle(spanningTree[dualEdge[0]][dualEdge[1]]['idx'])
        # Find the corresponding half-edge in the starting triangle
        foldingHalfEdge = mesh.halfedge_handle(foldingEdge, 0)
        if not (mesh.face_handle(foldingHalfEdge).idx() == dualEdge[0]):
            foldingHalfEdge = mesh.halfedge_handle(foldingEdge, 1)

        # Find the corresponding unwrapped half edge
        unfoldedLastHalfEdge = unfoldedMesh.halfedge_handle(halfEdgeConnections[foldingHalfEdge.idx()])

        # Find the point in the unwrapped triangle that is not on the fold edge
        oppositeUnfoldedVertex = unfoldedMesh.to_vertex_handle(unfoldedMesh.next_halfedge_handle(unfoldedLastHalfEdge))

        # We turn the half edges over to lie in the new triangle
        foldingHalfEdge = mesh.opposite_halfedge_handle(foldingHalfEdge)
        unfoldedLastHalfEdge = unfoldedMesh.opposite_halfedge_handle(unfoldedLastHalfEdge)

        # The two corners of the folded edge
        unfoldedFromVertex = unfoldedMesh.from_vertex_handle(unfoldedLastHalfEdge)
        unfoldedToVertex = unfoldedMesh.to_vertex_handle(unfoldedLastHalfEdge)

        # Calculate the edge lengths in the new triangle
        secondHalfEdgeInFace = mesh.next_halfedge_handle(foldingHalfEdge)
        thirdHalfEdgeInFace = mesh.next_halfedge_handle(secondHalfEdgeInFace)

        originalHalfEdges = [foldingHalfEdge, secondHalfEdgeInFace, thirdHalfEdgeInFace]

        edgelengths = [mesh.calc_edge_length(foldingHalfEdge), mesh.calc_edge_length(secondHalfEdgeInFace),
                       mesh.calc_edge_length(thirdHalfEdgeInFace)]

        # We calculate the two possibilities for the third point in the triangle
        [newUnfoldedVertex0, newUnfoldedVertex1] = getThirdPoint(unfoldedMesh.point(unfoldedFromVertex),
                                                                 unfoldedMesh.point(unfoldedToVertex), edgelengths[0],
                                                                 edgelengths[1], edgelengths[2])


        newUnfoldedVertex = unfoldedMesh.add_vertex(newUnfoldedVertex0)

        # Make the face
        newface = unfoldedMesh.add_face(unfoldedFromVertex, unfoldedToVertex, newUnfoldedVertex)

        secondUnfoldedHalfEdge = unfoldedMesh.next_halfedge_handle(unfoldedLastHalfEdge)
        thirdUnfoldedHalfEdge = unfoldedMesh.next_halfedge_handle(secondUnfoldedHalfEdge)
        unfoldedHalfEdges = [unfoldedLastHalfEdge, secondUnfoldedHalfEdge, thirdUnfoldedHalfEdge]

        # Store edge and side information
        # Dashed line in the output
        unfoldedLastEdge = unfoldedMesh.edge_handle(unfoldedLastHalfEdge)
        isFoldingEdge[unfoldedLastEdge.idx()] = True

        # Glue number and fold direction
        addVisualisationData(mesh, unfoldedMesh, originalHalfEdges, unfoldedHalfEdges, glueNumber, foldingDirection)

        # Related page
        connections[newface.idx()] = dualEdge[1]

        # Identify the half edges
        for i in range(3):
            halfEdgeConnections[originalHalfEdges[i].idx()] = unfoldedHalfEdges[i].idx()

    return [unfoldedMesh, isFoldingEdge, connections, glueNumber, foldingDirection]


def unfold(mesh):
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

    # TODO: add tabu search to resolve intersections
    # TODO: tabu search search for best spanning tree

    # Calculate the initial spanning tree
    spanningTree = randomSpanningTree(dualGraph)
    # spanningTree = nx.minimum_spanning_tree(dualGraph)

    # Unfold the tree
    fullUnfolding = unfoldSpanningTree(mesh, spanningTree)
    [unfoldedMesh, isFoldingEdge, connections, glueNumber, foldingDirection] = fullUnfolding

    # While spanning tree has overlaps we try to change tree by moving sub-trees/nodes to another parents
    # If after all possible changes spanning tree still has overlaps we take other random spanning tree
    overlaps, _ = detectOverlaps(unfoldedMesh, spanningTree)
    while overlaps:
        [spanningTree, overlaps] = tabuUnfolding(unfoldedMesh, spanningTree)
        if overlaps == 0:
            break

        # Case when overlaps > 0
        else:
            spanningTree = randomSpanningTree(dualGraph)
            fullUnfolding = unfoldSpanningTree(mesh, spanningTree)
            [unfoldedMesh, isFoldingEdge, connections, glueNumber, foldingDirection] = fullUnfolding
            overlaps, _ = detectOverlaps(unfoldedMesh, spanningTree)

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
    connectedComponents = nx.algorithms.components.connected_components(spanningTree)
    connectedComponentList = list(connectedComponents)

    # Processing of the components
    unfoldings = []
    for component in connectedComponentList:
        unfoldings.append(unfoldSpanningTree(mesh, spanningTree.subgraph(component)))

    return fullUnfolding, unfoldings

