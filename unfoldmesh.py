from tabuunfold import *


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
    # print(len(mesh.edges()))
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
        # dualGraph.add_node(face1.idx(), pos=center1)        # First face
        # dualGraph.add_node(face2.idx(), pos=center2)        # Second face
        dualGraph.add_node(face1.idx(), pos=center1, idx=face1.idx())       # First face
        dualGraph.add_node(face2.idx(), pos=center2, idx=face2.idx())       # Second face
        dualGraph.add_edge(face1.idx(), face2.idx(), idx=edge.idx(), weight=edgeweight)     # Edge between two faces
        e = list()
        e.append(face1.idx())
        e.append(face2.idx())
        e.append(edge.idx())
        edges_list.append(e)

    print(f'processing tree #{count}')
    spanningTree, overlaps = tabuUnfolding(mesh, randomSpanningTree(dualGraph), edges_list)
    while overlaps:
        print(f'No solution found. Remaining overlaps: {overlaps}')
        print('')
        count += 1
        print(f'processing tree #{count}')
        spanningTree, overlaps = tabuUnfolding(mesh, randomSpanningTree(dualGraph), edges_list)

    fullUnfolding = unfoldSpanningTree(mesh, spanningTree)
    [unfoldedMesh, isFoldingEdge, connections, glueNumber, foldingDirection] = fullUnfolding

    connectedComponents = nx.algorithms.components.connected_components(spanningTree)
    connectedComponentList = list(connectedComponents)

    # Processing of the components
    unfoldings = []
    for component in connectedComponentList:
        unfoldings.append(unfoldSpanningTree(mesh, spanningTree.subgraph(component)))

    return fullUnfolding, unfoldings

