import numpy as np
import openmesh as om
import networkx as nx
import random
import math
from collections import deque


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


class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, u):
        if self.parent[u] != u:
            self.parent[u] = self.find(self.parent[u])
        return self.parent[u]

    def union(self, u, v):
        root_u = self.find(u)
        root_v = self.find(v)
        if root_u != root_v:
            if self.rank[root_u] > self.rank[root_v]:
                self.parent[root_v] = root_u
            elif self.rank[root_u] < self.rank[root_v]:
                self.parent[root_u] = root_v
            else:
                self.parent[root_v] = root_u
                self.rank[root_u] += 1


# Compute the third point of a triangle when two points and all edge lengths are given
def getThirdPoint(v0, v1, l01, l12, l20):
    v2rotx = (l01 ** 2 + l20 ** 2 - l12 ** 2) / (2 * l01)
    v2roty0 = np.sqrt((l01 + l20 + l12) * (l01 + l20 - l12) * (l01 - l20 + l12) * (-l01 + l20 + l12)) / (2 * l01)

    v2roty1 = - v2roty0

    theta = np.arctan2(v1[1] - v0[1], v1[0] - v0[0])

    v2trans0 = np.array(
        [v2rotx * np.cos(theta) - v2roty0 * np.sin(theta), v2rotx * np.sin(theta) + v2roty0 * np.cos(theta), 0])
    v2trans1 = np.array(
        [v2rotx * np.cos(theta) - v2roty1 * np.sin(theta), v2rotx * np.sin(theta) + v2roty1 * np.cos(theta), 0])
    return [v2trans0 + v0, v2trans1 + v0]


# Check if two lines intersect
def lineIntersection(v1, v2, v3, v4, epsilon):
    d = (v4[1] - v3[1]) * (v2[0] - v1[0]) - (v4[0] - v3[0]) * (v2[1] - v1[1])
    u = (v4[0] - v3[0]) * (v1[1] - v3[1]) - (v4[1] - v3[1]) * (v1[0] - v3[0])
    v = (v2[0] - v1[0]) * (v1[1] - v3[1]) - (v2[1] - v1[1]) * (v1[0] - v3[0])
    if d < 0:
        u, v, d = -u, -v, -d
    return ((0 + epsilon) <= u <= (d - epsilon)) and ((0 + epsilon) <= v <= (d - epsilon))


# Check if a point lies inside a triangle
def pointInTriangle(A, B, C, P, epsilon):
    v0 = [C[0] - A[0], C[1] - A[1]]
    v1 = [B[0] - A[0], B[1] - A[1]]
    v2 = [P[0] - A[0], P[1] - A[1]]
    cross = lambda u, v: u[0] * v[1] - u[1] * v[0]
    u = cross(v2, v0)
    v = cross(v1, v2)
    d = cross(v1, v0)
    if d < 0:
        u, v, d = -u, -v, -d
    return u >= (0 + epsilon) and v >= (0 + epsilon) and (u + v) <= (d - epsilon)


# Check if two triangles intersect
def triangleIntersection(t1, t2, epsilon):
    if lineIntersection(t1[0], t1[1], t2[0], t2[1], epsilon): return True
    if lineIntersection(t1[0], t1[1], t2[0], t2[2], epsilon): return True
    if lineIntersection(t1[0], t1[1], t2[1], t2[2], epsilon): return True
    if lineIntersection(t1[0], t1[2], t2[0], t2[1], epsilon): return True
    if lineIntersection(t1[0], t1[2], t2[0], t2[2], epsilon): return True
    if lineIntersection(t1[0], t1[2], t2[1], t2[2], epsilon): return True
    if lineIntersection(t1[1], t1[2], t2[0], t2[1], epsilon): return True
    if lineIntersection(t1[1], t1[2], t2[0], t2[2], epsilon): return True
    if lineIntersection(t1[1], t1[2], t2[1], t2[2], epsilon): return True
    inTri = True
    inTri = inTri and pointInTriangle(t1[0], t1[1], t1[2], t2[0], epsilon)
    inTri = inTri and pointInTriangle(t1[0], t1[1], t1[2], t2[1], epsilon)
    inTri = inTri and pointInTriangle(t1[0], t1[1], t1[2], t2[2], epsilon)
    if inTri == True: return True
    inTri = True
    inTri = inTri and pointInTriangle(t2[0], t2[1], t2[2], t1[0], epsilon)
    inTri = inTri and pointInTriangle(t2[0], t2[1], t2[2], t1[1], epsilon)
    inTri = inTri and pointInTriangle(t2[0], t2[1], t2[2], t1[2], epsilon)
    if inTri == True: return True
    return False


# Functions for visualization and output
def addVisualisationData(mesh, unfoldedMesh, originalHalfedges, unfoldedHalfedges, glueNumber, foldingDirection):
    for i in range(3):
        # Folding direction
        if mesh.calc_dihedral_angle(originalHalfedges[i]) < 0:
            foldingDirection[unfoldedMesh.edge_handle(unfoldedHalfedges[i]).idx()] = -1
        else:
            foldingDirection[unfoldedMesh.edge_handle(unfoldedHalfedges[i]).idx()] = 1

        # Information, which edges belong together
        glueNumber[unfoldedMesh.edge_handle(unfoldedHalfedges[i]).idx()] = mesh.edge_handle(originalHalfedges[i]).idx()


# Find bounding box of face
def findBoundingBox(mesh):
    firstpoint = mesh.point(mesh.vertex_handle(0))
    xmin = firstpoint[0]
    xmax = firstpoint[0]
    ymin = firstpoint[1]
    ymax = firstpoint[1]
    for vertex in mesh.vertices():
        coordinates = mesh.point(vertex)
        if (coordinates[0] < xmin):
            xmin = coordinates[0]
        if (coordinates[0] > xmax):
            xmax = coordinates[0]
        if (coordinates[1] < ymin):
            ymin = coordinates[1]
        if (coordinates[1] > ymax):
            ymax = coordinates[1]
    boxSize = np.maximum(np.abs(xmax - xmin), np.abs(ymax - ymin))

    return [xmin, ymin, boxSize]


def writeSVG(filename, unfolding, size, printNumbers):
    mesh = unfolding[0]     # unfolding[0] = unfoldedMesh

    # print("write svg start")
    # face_handle = mesh.face_handle(20)
    # edges = []
    # for he in mesh.fh(face_handle):
    #     edge_handle = mesh.edge_handle(he)
    #     edges.append(edge_handle)
    # for e in edges:
    #     print(e.idx())

    isFoldingEdge = unfolding[1]
    glueNumber = unfolding[3]
    foldingDirection = unfolding[4]
    # Calculate the bounding box
    [xmin, ymin, boxSize] = findBoundingBox(unfolding[0])

    if size > 0:
        boxSize = size

    strokewidth = 0.002 * boxSize
    dashLength = 0.008 * boxSize
    spaceLength = 0.02 * boxSize

    textDistance = 0.02 * boxSize
    textStrokewidth = 0.05 * strokewidth
    textLength = 0.001 * boxSize
    fontsize = 0.015 * boxSize

    frame = 0.05 * boxSize

    # Open file in write mode (write)
    file = open(filename, 'w')

    # Write xml header
    file.write("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n")

    # The paper size and the scaling
    file.write("<svg width=\"30.5cm\" height=\"30.5cm\" viewBox = \"" + str(xmin - frame) + " " + str(
        ymin - frame) + " " + str(boxSize + 2 * frame) + " " + str(
        boxSize + 2 * frame) + "\" version = \"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n")

    # Go over all edges of the grid
    for edge in mesh.edges():
        # The two end points
        # if edge.idx() == 40:
        #     continue
        he = mesh.halfedge_handle(edge, 0)
        vertex0 = mesh.point(mesh.from_vertex_handle(he))
        vertex1 = mesh.point(mesh.to_vertex_handle(he))

        # Write a line between the two corners
        file.write("<path d =\"M " + str(vertex0[0]) + "," + str(vertex0[1]) + " " + str(vertex1[0]) + "," + str(
            vertex1[1]) + "\" style=\"fill:none;stroke:")

        # Color depending on fold direction
        if foldingDirection[edge.idx()] > 0:
            file.write("#ff0000")
        elif foldingDirection[edge.idx()] < 0:
            file.write("#0066ff")

        file.write(";stroke-width:" + str(
            strokewidth) + ";stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:4;stroke-dasharray:")

        # Dashed lines for fold edges
        if isFoldingEdge[edge.idx()]:
            file.write((str(dashLength) + ", " + str(spaceLength)))
        else:
            file.write("none")

        file.write(";stroke-dashoffset:0;stroke-opacity:1")
        file.write("\" />\n")

        # The number of the edge to be glued together
        if not isFoldingEdge[edge.idx()]:
            # Find halfedge in the face
            halfEdge = mesh.halfedge_handle(edge, 0)
            if mesh.face_handle(halfEdge).idx() == -1:
                halfEdge = mesh.opposite_halfedge_handle(halfEdge)
            vector = mesh.calc_edge_vector(halfEdge)
            # normalize
            vector = vector / np.linalg.norm(vector)
            midPoint = 0.5 * (
                    mesh.point(mesh.from_vertex_handle(halfEdge)) + mesh.point(mesh.to_vertex_handle(halfEdge)))
            rotatedVector = np.array([-vector[1], vector[0], 0])
            angle = np.arctan2(vector[1], vector[0])
            position = midPoint + textDistance * rotatedVector
            rotation = 180 / np.pi * angle

            if (printNumbers):
                file.write("<text x=\"" + str(position[0]) + "\" y=\"" + str(position[1]) + "\" font-size=\"" + str(
                    fontsize) + "\" stroke-width=\"" + str(textStrokewidth) + "\" transform=\"rotate(" + str(
                    rotation) + "," + str(position[0]) + "," + str(position[1]) + ")\">" + str(
                    glueNumber[edge.idx()]) + "</text>\n")

    file.write("</svg>")
    file.close()


# Function to generate all possible panning trees of graph
# Proccess_function is function that we run on each face
# In this project process_function = tabuunfold
def all_spanning_trees(mesh, graph, process_function, edges_list, t=1):
    def dfs(subgraph, edges, visited, t):
        if len(subgraph.nodes) == len(graph.nodes):
            print(f'processing Tree #{t}')
            result = process_function(mesh, subgraph.copy(), edges_list)
            if not result:  # If found unfold without overlaps
                return subgraph.copy()
            return False
        for u, v, data in edges:
            if u not in visited or v not in visited or not nx.has_path(subgraph, u, v):
                new_visited = visited | {u, v}
                subgraph.add_edge(u, v, **data)
                remaining_edges = [e for e in edges if e != (u, v, data)]
                found_tree = dfs(subgraph, remaining_edges, new_visited, t)
                if found_tree:
                    return found_tree
                subgraph.remove_edge(u, v)
        return False

    nodes = list(graph.nodes)
    if not nodes:
        return None

    initial_edge_list = list(graph.edges(data=True))
    initial_visited = {nodes[0]}
    initial_subgraph = nx.Graph()
    initial_subgraph.add_node(nodes[0])

    return dfs(initial_subgraph, initial_edge_list, initial_visited, t)


# Function to generate random panning tree of graph
def randomSpanningTree(dual_graph):
    edges = list(dual_graph.edges(data=True))
    random.shuffle(edges)
    uf = UnionFind(len(dual_graph.nodes))
    mst = dual_graph.__class__()
    mst.graph.update(dual_graph.graph)
    mst.add_nodes_from(dual_graph.nodes.items())
    for u, v, data in edges:
        if uf.find(u) != uf.find(v):
            mst.add_edge(u, v, **data)
            uf.union(u, v)
            if mst.number_of_edges() == len(dual_graph.nodes) - 1:
                break

    return mst


# Function to change root of a tree
# def changeRoot(tree, newRoot):
#     def dfs(cur_node, parent):
#         children = []
#         for n in tree.neighbors(cur_node):
#             if n != parent:
#                 children.append(n)
#                 dfs(n, cur_node)
#         for c in children:
#             tree.add_edge(cur_node, c)
#     new_tree = nx.Graph()
#     new_tree.add_nodes_from(tree.nodes())
#     dfs(newRoot, None)
#     return tree


# Returns root of current tree
def getRoot(tree):
    # Get the root of the tree
    return tree.graph['root']


# Set root of current tree to specific node
def setRoot(tree, new_root):
    tree.graph['root'] = new_root
    tree.nodes[new_root]['parent'] = None
    return tree


# Function to change root of a tree
def changeTreeRoot(tree, new_root):
    if not tree.has_node(new_root):
        print(new_root)
        raise ValueError("The specified new root is not in the tree.")

        # Create an empty graph for the new tree
    new_tree = nx.Graph()

    # Perform DFS to set the parent-child relationships
    def dfs(current_node, parent):
        # Set the parent attribute for the current node
        if parent is not None:
            tree.nodes[current_node]['parent'] = parent
        else:
            tree.nodes[current_node]['parent'] = None  # Root node has no parent

        for neighbor in tree.neighbors(current_node):
            if neighbor != parent:
                # Add the edge to the new tree
                edge_data = tree.get_edge_data(current_node, neighbor)
                new_tree.add_edge(current_node, neighbor, **edge_data)
                # Continue DFS traversal
                dfs(neighbor, current_node)

    # Start DFS from the new root
    dfs(new_root, None)

    # Copy the node attributes (including 'parent') to the new tree
    for node in new_tree.nodes:
        new_tree.nodes[node].update(tree.nodes[node])
    setRoot(new_tree, new_root)
    return new_tree


# Return a subtree of given node
def getSubTree(tree, root):
    # Use DFS to traverse the subtree from the root node
    def dfs(current_node, parent, subtree_nodes):
        subtree_nodes.add(current_node)
        for neighbor in tree.neighbors(current_node):
            if neighbor != parent:
                dfs(neighbor, current_node, subtree_nodes)

    # Set to keep track of nodes in the subtree
    subtree_nodes = set()

    # Start DFS from the root node
    parent = tree.nodes[root]['parent']
    dfs(root, parent, subtree_nodes)
    subtree_nodes.remove(root)
    subtree_nodes.add(parent)
    return subtree_nodes


def calculateBoundingBox(triangle):
    min_corner = np.min(triangle, axis=0)
    max_corner = np.max(triangle, axis=0)
    return min_corner, max_corner


def hashFunction(point, cell_size):
    return tuple((point // cell_size).astype(int))


def insertIntoGrid(grid, face, bounding_box, cell_size):
    min_hash = hashFunction(bounding_box[0], cell_size)
    max_hash = hashFunction(bounding_box[1], cell_size)
    for x in range(min_hash[0], max_hash[0] + 1):
        for y in range(min_hash[1], max_hash[1] + 1):
            for z in range(min_hash[2], max_hash[2] + 1):
                cell = (x, y, z)
                if cell not in grid:
                    grid[cell] = []
                grid[cell].append(face)


def detectOverlaps(unfoldedMesh, spanningTree, cell_size=1.0):
    epsilon = 1E-12  # accuracy
    collidingFaces = set()

    faces = list(unfoldedMesh.faces())
    num_faces = len(faces)

    grid = {}
    bounding_boxes = {}

    # Calculate bounding boxes and insert faces into the grid
    for face in faces:
        triangle = [np.array(unfoldedMesh.point(unfoldedMesh.from_vertex_handle(halfedge))) for halfedge in unfoldedMesh.fh(face)]
        bounding_box = calculateBoundingBox(triangle)
        bounding_boxes[face.idx()] = bounding_box
        insertIntoGrid(grid, face, bounding_box, cell_size)

    # Detect overlaps
    for face1 in faces:
        triangle1 = [np.array(unfoldedMesh.point(unfoldedMesh.from_vertex_handle(halfedge))) for halfedge in unfoldedMesh.fh(face1)]
        bounding_box1 = bounding_boxes[face1.idx()]
        min_hash1 = hashFunction(bounding_box1[0], cell_size)
        max_hash1 = hashFunction(bounding_box1[1], cell_size)

        checked_faces = set()
        for x in range(min_hash1[0], max_hash1[0] + 1):
            for y in range(min_hash1[1], max_hash1[1] + 1):
                for z in range(min_hash1[2], max_hash1[2] + 1):
                    cell = (x, y, z)
                    if cell in grid:
                        for face2 in grid[cell]:
                            if face2.idx() <= face1.idx() or face2.idx() in checked_faces:
                                continue
                            checked_faces.add(face2.idx())
                            triangle2 = [np.array(unfoldedMesh.point(unfoldedMesh.from_vertex_handle(halfedge)))
                                         for halfedge in unfoldedMesh.fh(face2)]
                            if triangleIntersection(triangle1, triangle2, epsilon):
                                collidingFaces.add(face1.idx())
                                collidingFaces.add(face2.idx())

    return len(collidingFaces), collidingFaces


def print_loading_bar(current_count, max_count, overlaps):
    # Calculate the percentage
    percentage = (current_count / max_count) * 100
    # Calculate the number of '#' to print
    num_hashes = int(percentage * 0.5)
    # Calculate the number of '_' to print
    num_underscores = 50 - num_hashes

    if overlaps == 0:
        num_hashes = 50
        num_underscores = 0
        percentage = 100

    # Print the loading bar
    print(f"\rProcessing [{'#' * num_hashes}{'_' * num_underscores} {int(percentage)}%] Overlaps: {overlaps}", end='')

