import numpy as np
import openmesh as om
import networkx as nx
import random
from unfoldmesh import *


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


def faceBoundingBox(mesh, face):
    vertices = [v for v in mesh.fv(face)]
    firstpoint = mesh.point(vertices[0])
    xmin = firstpoint[0]
    xmax = firstpoint[0]
    ymin = firstpoint[1]
    ymax = firstpoint[1]
    for vertex in vertices:
        coordinates = mesh.point(vertex)
        if coordinates[0] < xmin:
            xmin = coordinates[0]
        if coordinates[0] > xmax:
            xmax = coordinates[0]
        if coordinates[1] < ymin:
            ymin = coordinates[1]
        if coordinates[1] > ymax:
            ymax = coordinates[1]
    # xmin - Lower Bounding Box coordinate. found at index 0
    # xmax - Upper Bounding Box coordinate found at index 2
    # same for ymin and ymax if necessary. indexes 1 and 3
    return (xmin, ymin, xmax, ymax)


def writeSVG(filename, unfolding, size, printNumbers):
    mesh = unfolding[0]
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


def all_spanning_trees(graph):
    def dfs(subgraph, edges, visited):
        if len(subgraph.nodes) == len(graph.nodes):
            yield subgraph.copy()
            return
        for u, v, data in edges:
            if u not in visited or v not in visited:
                new_visited = visited.copy()
                new_visited.add(u)
                new_visited.add(v)
                subgraph.add_edge(u, v, **data)
                new_edges = [e for e in edges if e != (u, v, data)]
                yield from dfs(subgraph, new_edges, new_visited)
                subgraph.remove_edge(u, v)

    nodes = list(graph.nodes)
    if not nodes:
        return

    initial_edge_list = list(graph.edges(data=True))
    initial_visited = {nodes[0]}
    initial_subgraph = nx.Graph()
    initial_subgraph.add_node(nodes[0])

    yield from dfs(initial_subgraph, initial_edge_list, initial_visited)


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


def detectOverlaps(unfoldedMesh, spanningTree):
    overlaps = 0
    epsilon = 1E-12  # accuracy
    faceBoxes = []
    collidingFaces = set()
    for face in unfoldedMesh.faces():
        bbox = faceBoundingBox(unfoldedMesh, face)
        faceBoxes.append((face, bbox))
    faceBoxes.sort(key=lambda fb: fb[1][0])
    n = len(faceBoxes)
    for i in range(n):
        fi, boxi = faceBoxes[i]
        for j in range(i+1, n):
            fj, boxj = faceBoxes[j]
            if boxj[0] > boxi[2]:
                break
            # Check if fi and fj are connected in the spanning tree
            if nx.has_path(spanningTree, fi.idx(), fj.idx()):
                triangle1 = []
                triangle2 = []
                for halfedge in unfoldedMesh.fh(fi):
                    triangle1.append(unfoldedMesh.point(unfoldedMesh.from_vertex_handle(halfedge)))
                for halfedge in unfoldedMesh.fh(fj):
                    triangle2.append(unfoldedMesh.point(unfoldedMesh.from_vertex_handle(halfedge)))
                if triangleIntersection(triangle1, triangle2, epsilon):
                    overlaps += 1
                    collidingFaces.add(fi.idx())
                    collidingFaces.add(fj.idx())
    return overlaps, collidingFaces

