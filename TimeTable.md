Work process:

week 1 (21-27/4): First of all took me about weak and half to run this project.
    There was a bug with openmesh. Openmesh did not want to install. There was an error with cmake for some reason.
    I tried to run it on vs, vs code, pycharm. Didn't help. Tried to build cmake solution. But it works for c++ projects.
    After all the problem was due to python version. You can download only last version of python(3.12) from official website.
    Found proper version of python(3.8) in microsoft store.
    Note: this project works only with python 3.8. in python 3.9 they moved function from math to other library, so it will lead to errors.

week 2 (28/4- 4/5): Started to inspect code flow, check how it works. Added some model from internet to see the result.
    This project take mst of dual_graph. unfold mesh by this mst. and for all overlaps it divides model to sub-part
    with no overlaps(connected components)

week 3 (5-11/5): Removed code that divides model to sub-part. It started to output only one unfolded mesh(no sub-parts)
    with overlaps(visible in svg file).
    First task was to implement function to calculate spaning tree on dual-graph.
    I had few approaches: calculate randon spaning tree, calculate all possible spaning trees
    First approach: all spanning tree function(in usefullfunctions.py) return list of all spanning trees of graph.
    But run-time is too high, so i decided to go to next approach.
    Second approach: random spanning tree. The main idea is to pick every time random edge from edge list,
    check if after adding this edge to spanning tree, there is no cycle.
    This approach more simple and the chance that it will generate two identical spanning trees is negligible
    But maybe I will change all_spanning_trees to fetch my task.

week 4 (12-18/5): Read the "Unfolding Polyhedra via Tabu Search" report. The main goal of this report was to implement 
    algorithm that take spanning tree and try to re-attach overlaping faces to their neighbors.
    To do this need to implement change to tabu search algorithm to fit the problem.
    First I need to implement detectOverlaps function. The main idea was to sort faces by lower x coordinate and iterate
    over all faces and check if there is overlaps(by using triangleIntersection function)
    Tried to do simple overlaps check b using two loop( for in for) and check intersection for all faces,
    but it works about three time slower, than function that use sorted faces.
    After I had an idea to find best spanning tree with the least number of overlaps.
    For this I implemented find_best_spanning_tree function, but after all it is not guaranteed that this best tree
    will lead to solution

week 5 (19-25/5): First moved a lot of function from unfoldmesh.py to usefullfunctions.py to clean the code a bit(for visual convenience only)
    Started to implement tabu unfold algorithm.
    First added the main skeleton of algorithm with names of all function to be implemented.
    Name of functions taken from report.
    Started with simple function. Added setInitialRoot function which take first node of tree and set it as root.
    Currently, the root not so import cause we gonna change root in every iteration of main while loop
    Implemented randomReroot function which take random node and set it as root of the tree.
    Plus added getRoot and setRoot function to set and get root of the tree

week 6 (26/5 - 1/7): Added selectRandomCollidingFace function. For this changed detectOverlaps function to bring and
    return list of overlaping faces
    Added getNeighbors and bestNeighbor function. getNeighbors return list of all neighbors of given face,
    bestNeighbor return first neighbor in the list that is not in tabu list
    
week 7 (2-8/6): Added calculateAverageValence, memorize, initHistoryValues function needed to tabu-unfold
    Started to implement move function which take random overlaping face and try to re-attach it to best neighbor.

week 8 (9-15/7): Continued to work on move, the function work and re-attach node to it's neighbor, but currently there 
    is a bug when I add new edge, need to add **data as well. need to resolve it, and unfold mesh using new spanning tree.

TO-DO:
week 9 (16-22/6): finish with move, implement stuckInRootLock, stuckInMemoryLock,rerootIntoStuckSubTree functions
    check algorithm for bugs and how it works

week 10 (23-30/6): finish with algorithm and check how it works on different models



Work time:  I don't want to take into account first week, but it took me more than 10 hours to resolve the problem
            for week 2-8 i work about 2-3 days a week for 4-8  hours. I guess I will calculate average work time per week
            as 2.5 days and 6 hours per for total 15 hours. but let's take 13 hours. I total for weeks 2-8 it will be 91 hours.
            
            For weeks 9-10 plan to do 20-25 hours and finish with project

