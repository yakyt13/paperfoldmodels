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

# Di-graph doesn't have parent-child dependency
# This function add 'parent' field to each node
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

# Used to determine tabu-list length
def calculateAverageValence(spanningTree):
    return sum(spanningTree.degree(face) for face in spanningTree.nodes()) / len(spanningTree.nodes())


# Check stuck situation when tree root overlapping
def stuckInRootLock(spanningTree, edge_list, collidingFaces, root=None):
    if root is None:
        root = getRoot(spanningTree)
    if root not in collidingFaces:
        return False
    neighbors = getNeighbors(edge_list, root)
    neighbors = set(neighbors)
    sub_tree = getSubTree(spanningTree, root)
    # sub_tree = dfs_tree(spanningTree, root)
    return neighbors.issubset(sub_tree)


# Check if there is no available moves
def stuckInMemoryLock(spanningTree, tabu_list, collidingFaces, edges_list):
    for face in collidingFaces:
        if availableNeighbors(spanningTree, face, tabu_list, edges_list):
            return True
    return False


# Clear tabu-list
def clearMemory(tabu_list):
    tabu_list.clear()
    return tabu_list


# Change root of stuck subtree(whole tree in our case)
def rerootIntoStuckSubTree(spanningTree):
    root = getRoot(spanningTree)
    sub_tree = getSubTree(spanningTree, root)
    # sub_tree = dfs_tree(spanningTree, root)
    if None in sub_tree:
        sub_tree.remove(None)
    new_root = random.choice(list(sub_tree))
    return changeTreeRoot(spanningTree, new_root)


# Change tree root to given node
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


# Perform random reroot of tree
def randomReroot(tree):
    nodes = list(tree.nodes)
    root = getRoot(tree)
    nodes.remove(root)
    root = random.choice(nodes)
    return changeTreeRoot(tree, root)


# Print tree in DFS order
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


# Initialize history values. used in Tabu unfold algorithm
def initHistoryValues(tree, collidingFaces):
    # Initialize history values
    bH, xH = None, None
    treeH = tree
    oH = float('inf')
    collH = collidingFaces
    return bH, xH, oH, treeH, collH


# Pick a random colliding  face
def selectRandomCollidingFace(collidingFaces):
    if not list(collidingFaces):
        return None
    return random.choice(list(collidingFaces))


# Return list of neighbors of given face
def getNeighbors(edges_list, face_idx):
    neighbors = list()
    for edge in edges_list:
        edge = list(edge)
        if edge[0] == face_idx:
            neighbors.append(edge[1])
        if edge[1] == face_idx:
            neighbors.append(edge[0])
    return neighbors


# Return if there is available neighbors of given face
def availableNeighbors(spanningTree, face, tabu_list, edges_list):
    neighbors = getNeighbors(edges_list, face)
    neighbors = set(neighbors)
    sub_tree = getSubTree(spanningTree, face)
    # sub_tree = dfs_tree(spanningTree, face)
    for neighbor in neighbors:
        if neighbor not in sub_tree and {face, neighbor} not in tabu_list:
            return True
    return False


# Pick a random neighbor of fiven face
def pickNeighbor(spanningTree, face, tabu_list, edges_list):
    neighbors = getNeighbors(edges_list, face)
    sub_tree = getSubTree(spanningTree, face)
    # sub_tree = dfs_tree(spanningTree, face)
    result = [item for item in neighbors if item not in sub_tree]
    if not result:
        return None
    return random.choice(result)


# Return best neighbor of given face
def bestNeighbor(mesh, spanningTree, face, tabu_list, edges_list):
    neighbors = getNeighbors(edges_list, face)
    neighbors = set(neighbors)
    sub_tree = getSubTree(spanningTree, face)
    # sub_tree = dfs_tree(spanningTree, face)

    best_neighbor = None
    best_score = float('inf')
    best_tree = None
    best_unfold = None
    best_col_faces = None

    for neighbor in neighbors:
        if neighbor not in sub_tree and {face, neighbor} not in tabu_list:
            temp_tree = move(spanningTree, face, neighbor, edges_list)
            fullUnfolding = unfoldSpanningTree(mesh, temp_tree)
            unfoldedMesh = fullUnfolding[0]
            temp_overlaps, temp_col_faces = detectOverlaps(unfoldedMesh, temp_tree)

            if temp_overlaps < best_score:
                best_score = temp_overlaps
                best_neighbor = neighbor
                best_tree = temp_tree
                best_unfold = fullUnfolding
                best_col_faces = temp_col_faces

    return best_tree, best_score, best_col_faces, best_neighbor, best_unfold


# Perform a move. Attach face to given neighbor
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


# Memorize performed move (face, neighbor)
def memorize(face, neighbor, tabu_list, max_iterations):
    # Memorize the move
    if face is not None and {face, neighbor} not in tabu_list:
        tabu_list.append({face, neighbor})
    if len(tabu_list) > max_iterations:
        tabu_list.pop(0)


# Return parent of given face
def findParentOfFace(spaningTree, face):
    parent = spaningTree.nodes[face]['parent']
    return parent


# Return index of given face. each face in mesh(before unfolding) has personal index
def findEdgeIndex(p1, p2, edges_list):
    for edge in edges_list:
        edge = list(edge)
        if (edge[0] == p1 and edge[1] == p2) or (edge[0] == p2 and edge[1] == p1):
            return edge[2]
    return 0


# Main function to unfold mesh without overlaps
def tabuUnfolding(mesh, spanningTree, edges_list, count=0, max_count=None):
    # Main loop of the Tabu Search
    spanningTree = parentChildDependency(spanningTree)
    # initial unfolding and overlaps
    fullUnfolding = unfoldSpanningTree(mesh, spanningTree)
    [unfoldedMesh, isFoldingEdge, connections, glueNumber, foldingDirection] = fullUnfolding
    overlaps, collidingFaces = detectOverlaps(unfoldedMesh, spanningTree)

    if max_count is None:
        max_count = len(spanningTree.nodes) * 2
        # max_count = len(spanningTree.nodes)
    print(f'max iterations for current tree: {max_count}')

    tabu_list = []
    tabu_len = 3 * math.log(len(spanningTree.nodes), 3) + 1

    while overlaps > 0 and count < max_count:
    # while overlaps > 0:
        if stuckInRootLock(spanningTree, edges_list, collidingFaces):
            spanningTree = rerootIntoStuckSubTree(spanningTree)
        else:
            spanningTree = randomReroot(spanningTree)

        face = selectRandomCollidingFace(collidingFaces)
        face = connections[face]
        bH, xH, oH, treeH, collH = initHistoryValues(spanningTree, collidingFaces)
        temp_overlaps = 0

        if face is None:
            print("No overlaps found")
            break

        if stuckInMemoryLock(spanningTree, tabu_list, collidingFaces, edges_list):
            tabu_list = clearMemory(tabu_list)

        while face != getRoot(spanningTree):
            parent = findParentOfFace(spanningTree, face)
            temp_tree, temp_overlaps, temp_colliding_faces, b, fullUnfolding = (
                bestNeighbor(mesh, spanningTree, face, tabu_list, edges_list))  # Filtered by tabu list
            if b:
                [unfoldedMesh, _, connections, _, _] = fullUnfolding

            if not b:
                face = parent
                continue

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
        count += 1
        print_loading_bar(count, max_count, overlaps)
        # print(f'\rmove # {count}, overlap: {overlaps}', end='')
    print('')

    return spanningTree, overlaps
