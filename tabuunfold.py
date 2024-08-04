import math
import random
from collections import deque
import openmesh as om
from usefullfunctions import *


# Tabu search steps:
# 1) Input: mesh model. Each face must be triangular and planar
# 2) The initial unfolder: unfolded model. If it has overlaps it is acceptable
# 3) Selecting the best step: unfolded model with least number of overlaps
# 4) Local Minima: Tabu search store m past steps (tabu-list). Any step in the tabu list cannot be undone.
#                 need to choose value of m so that algorithm won't fall back to local minima
#                 and won't block possible moves.
#                 we choose m = val * log_val(|F|). val is the average valence in dual-graph.
# 5) Switching root node if needed
# 6) Efficient overlaps detection: After each change in unfold-tree it is necessary to determine
#                                  the resulting number of overlaps.
#   This part is most time-consuming, so we need to implement efficient algorithm to resolve this
#   To do this we will implement Sweep Line Overlap Detection algorithm which improve complexity to O(nlogn)

def parentChildDependency(tree):
    nodes = list(tree.nodes)
    root = nodes[0]
    nx.set_node_attributes(tree, None, 'parent')
    queue = deque([root])
    tree.nodes[root]['parent'] = None

    while queue:
        current = queue.popleft()
        for n in tree.neighbors(current):
            if tree.nodes[n]['parent'] is None and n != root:
                tree.nodes[n]['parent'] = current
                queue.append(n)
    tree.graph['root'] = root
    return tree


def calculateAverageValence(spanningTree):
    return sum(spanningTree.degree(face) for face in spanningTree.nodes()) / len(spanningTree.nodes())


def stuckInRootLock(spanningTree, edge_list, collidingFaces, root=None):
    if root is None:
        root = getRoot(spanningTree)
    if root not in collidingFaces:
        return False
    neighbors = getNeighbors(edge_list, root)
    neighbors = set(neighbors)
    sub_tree = getSubTree(spanningTree, root)
    return neighbors.issubset(sub_tree)


def stuckInMemoryLock(tabu_list, max_iterations):
    return len(tabu_list) > max_iterations


def rerootIntoStuckSubTree(spanningTree):
    root = getRoot(spanningTree)
    sub_tree = getSubTree(spanningTree, root)
    if None in sub_tree:
        sub_tree.remove(None)
    new_root = random.choice(list(sub_tree))
    if new_root is None:
        print(sub_tree)
    return changeTreeRoot(spanningTree, new_root)


def reRoot(tree, new_root):
    current = new_root
    parent = tree.nodes[current]['parent']

    while parent is not None:
        # Get the current node's parent and its parent
        grandparent = tree.nodes[parent]['parent']

        # Reverse the link
        tree.nodes[parent]['parent'] = current
        tree.nodes[current]['parent'] = grandparent

        # Move up the tree
        current = parent
        parent = grandparent

    # Update the new root node
    tree = setRoot(tree, new_root)
    return tree


def randomReroot(tree):
    nodes = list(tree.nodes)
    root = getRoot(tree)
    nodes.remove(root)
    root = random.choice(nodes)
    return changeTreeRoot(tree, root)


def printTreeDFS(tree):
    visited = set()
    root = getRoot(tree)

    def dfs(node):
        if node in visited:
            return
        visited.add(node)
        parent = tree.nodes[node]['parent']
        print(f'node: {node}, parent: {parent}')
        for neighbor in tree.neighbors(node):
            if neighbor not in visited:
                dfs(neighbor)

    dfs(root)


def initHistoryValues(tree, overlaps, collidingFaces):
    # Initialize history values
    bH, xH = None, None
    treeH = tree
    oH = float('inf')
    collH = collidingFaces
    return bH, xH, oH, treeH, collH


def selectRandomCollidingFace(collidingFaces):
    if not list(collidingFaces):
        return None
    return random.choice(list(collidingFaces))


def getNeighbors(edges_list, face_idx):
    neighbors = list()
    for edge in edges_list:
        edge = list(edge)
        if edge[0] == face_idx:
            neighbors.append(edge[1])
        if edge[1] == face_idx:
            neighbors.append(edge[0])
    return neighbors


def bestNeighbor(spanningTree, face,  tabu_list, edges_list):
    # Find the best neighbor filtered by the tabu list
    neighbors = getNeighbors(edges_list, face)
    neighbors = set(neighbors)
    sub_tree = getSubTree(spanningTree, face)
    for neighbor in neighbors:
        if neighbor not in sub_tree and {face, neighbor} not in tabu_list:
            return neighbor
    return None


def move(spanningTree, face, neighbor, edges_list):
    root = getRoot(spanningTree)
    predecessors = nx.predecessor(spanningTree, root)
    parent = predecessors.get(face, [])
    idx = findEdgeIndex(neighbor, face, edges_list)
    if parent and parent[0] is not None:
        spanningTree.remove_edge(parent[0], face)
    spanningTree.add_edge(neighbor, face, idx=idx)
    spanningTree.nodes[face]['parent'] = neighbor
    return spanningTree


def memorize(face, neighbor, tabu_list, max_iterations):
    # Memorize the move
    if face is not None and {face, neighbor} not in tabu_list:
        tabu_list.append({face, neighbor})
    if len(tabu_list) > max_iterations:
        tabu_list.pop(0)


def findParentOfFace(spaningTree, face):
    parent = spaningTree.nodes[face]['parent']
    return parent


def findEdgeIndex(p1, p2, edges_list):
    for edge in edges_list:
        edge = list(edge)
        if (edge[0] == p1 and edge[1] == p2) or (edge[0] == p2 and edge[1] == p1):
            return edge[2]
    return 0


def tabuUnfolding(mesh, spanningTree, edges_list, count=0, max_count=None):
    # Main loop of the Tabu Search
    spanningTree = parentChildDependency(spanningTree)
    # initial unfolding and overlaps
    fullUnfolding = unfoldSpanningTree(mesh, spanningTree)
    [unfoldedMesh, isFoldingEdge, connections, glueNumber, foldingDirection] = fullUnfolding
    overlaps, collidingFaces = detectOverlaps(unfoldedMesh, spanningTree)

    if max_count is None:
        max_count = len(spanningTree.nodes) * 2
    print(f'max iterations for current tree: {max_count}')

    tabu_list = []
    tabu_len = 3 * math.log(len(spanningTree.nodes), 3) + 1

    while overlaps > 0 and count < max_count:
    # while overlaps > 0:
        if overlaps > len(spanningTree.nodes):
            print(collidingFaces)
        if stuckInRootLock(spanningTree, edges_list, collidingFaces):
            spanningTree = rerootIntoStuckSubTree(spanningTree)
        else:
            spanningTree = randomReroot(spanningTree)

        face = selectRandomCollidingFace(collidingFaces)
        bH, xH, oH, treeH, collH = initHistoryValues(spanningTree, overlaps, collidingFaces)
        temp_overlaps = 0

        if face is None:
            print("No overlaps found")
            break

        # if stuckInMemoryLock
        # if not bestNeighbor(spanningTree, face, tabu_list, edges_list):
        #     tabu_list = []

        while face != getRoot(spanningTree):
            b = bestNeighbor(spanningTree, face, tabu_list, edges_list)  # Filtered by tabu list
            parent = findParentOfFace(spanningTree, face)
            if b is None:
                face = parent
                continue
            temp_tree = move(spanningTree, face, b, edges_list)
            fullUnfolding = unfoldSpanningTree(mesh, temp_tree)
            [unfoldedMesh, _, _, _, _] = fullUnfolding
            temp_overlaps, temp_colliding_faces = detectOverlaps(unfoldedMesh, temp_tree)

            if temp_overlaps < overlaps:
                spanningTree = temp_tree    # move(face, neighbor)
                overlaps = temp_overlaps
                collidingFaces = temp_colliding_faces
                memorize(face, b, tabu_list, tabu_len)
                break

            if temp_overlaps < oH:
                bH = b
                xH = face
                oH = temp_overlaps
                treeH = temp_tree
                collH = temp_colliding_faces

            face = parent

        # No improving solution found within the loop
        if temp_overlaps >= oH:
            spanningTree = treeH
            overlaps = oH
            collidingFaces = collH
            memorize(xH, bH, tabu_list, tabu_len)
        # if not count % 25:
        print(f"Performed {count} moves. Current overlaps: {overlaps}")
        # print('#', end='')
        count += 1
    print('')
    return spanningTree, overlaps
    # return overlaps
