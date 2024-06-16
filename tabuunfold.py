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

def setInitialRoot(tree):
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

def stuckInRootLock():
    # Implement a condition to check if stuck in root lock
    pass


def stuckInMemoryLock():
    # Implement a condition to check if stuck in memory lock
    pass


def rerootIntoStuckSubTree(x):
    # Implement rerooting into a stuck subtree
    pass


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
    # Select a random node that is not the current root
    nodes = list(tree.nodes)
    nodes.remove(tree.graph['root'])
    new_root = random.choice(nodes)

    # Reroot the tree
    tree = reRoot(tree, new_root)
    tree.graph['root'] = new_root
    return tree


def initHistoryValues():
    # Initialize history values
    bH, xH = None, None
    return bH, xH


def selectRandomCollidingFace(collidingFaces):
    return random.choice(list(collidingFaces))


def getNeighbors(unfoldedMesh, face_idx):

    neighbors = set()
    face_handle = unfoldedMesh.face_handle(face_idx)
    for halfedge_handle in unfoldedMesh.fh(face_handle):
        opp_halfedge_handle = unfoldedMesh.opposite_halfedge_handle(halfedge_handle)
        adj_face_handle = unfoldedMesh.face_handle(opp_halfedge_handle)
        if adj_face_handle.is_valid() and adj_face_handle != face_handle:
            neighbors.add(adj_face_handle.idx())

    return list(neighbors)


def bestNeighbor(neighbors, face,  tabu_list):
    # Find the best neighbor filtered by the tabu list
    for neighbor in neighbors:
        if (face, neighbor) not in tabu_list:
            return neighbor
    return None


def move(spaningTree, face, neighbor):
    root = getRoot(spaningTree)
    predecessors = nx.predecessor(spaningTree, root)
    parent = predecessors.get(face, [])
    if parent[0] is not None:
        spaningTree.remove_edge(parent[0], face)
    spaningTree.add_edge(neighbor, face)
    return spaningTree


def memorize(face, neighbor, tabu_list, max_iterations):
    # Memorize the move
    tabu_list.append((face, neighbor))
    if len(tabu_list) > max_iterations:
        tabu_list.pop(0)


def getRoot(tree):
    # Get the root of the tree
    return tree.graph['root']


def setRoot(tree, new_root):
    tree.graph['root'] = new_root
    tree.nodes[new_root]['parent'] = None
    return tree


def findParentOfFace(spaningTree, face):
    parent = spaningTree.nodes[face]['parent']
    return parent


def tabuUnfolding(unfoldedMesh, spaningTree):
    tabu_list = []
    max_iterations = 3 * math.log(len(unfoldedMesh.faces()), 3)
    history_values = {}
    overlaps, collidingFaces = detectOverlaps(unfoldedMesh, spaningTree)
    spaningTree = setInitialRoot(spaningTree)

    # Main loop of the Tabu Search
    # while overlaps > 0:
    #     if stuckInRootLock():
    #         rerootIntoStuckSubTree(spaningTree)
    #     else:
    #         randomReroot(spaningTree)
    #
    #     if stuckInMemoryLock():
    #         tabu_list.clear()
    #
    #     bH, xH = initHistoryValues()
    #     face = selectRandomCollidingFace(collidingFaces)
    #
    #     while face != getRoot(spaningTree):
    #         neighbor = bestNeighbor(face)
    #         if move(spaningTree, face, neighbor) < overlaps:
    #             spaningTree = move(spaningTree, face, neighbor)
    #             memorize(face, neighbor)
    #             break
    #
    #         if move(spaningTree, face, neighbor) < move(spaningTree, xH, bH):
    #             bH = neighbor
    #             xH = face
    #
    #         face = face.parent
    #
    #     # No improving solution found within the loop
    #     spaningTree = move(spaningTree, xH, bH)
    #     memorize(xH, bH, tabu_list, max_iterations)

    face = random.choice(list(collidingFaces))
    neighbors = getNeighbors(unfoldedMesh, face)
    best_neighbor = bestNeighbor(neighbors, face, tabu_list)
    new_spaning_tree = move(spaningTree, face, best_neighbor)

    print(f"All neighbors: {neighbors}")
    print(f"Best neighbor: {best_neighbor}")
    print(f"Face parent: {findParentOfFace(spaningTree, face)}")
    print('Root: ', getRoot(spaningTree))
    print(f"Face index: {selectRandomCollidingFace(collidingFaces)}")
    print(f"Average valence: {calculateAverageValence(spaningTree)}")
    print(f"Max iterations: {max_iterations}")
    return spaningTree, 0

