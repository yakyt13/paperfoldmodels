from argparse import ArgumentParser
from unfoldmesh import *


def main():
    parser = ArgumentParser()
    parser.add_argument("-f", "--file", dest="filename", default="models/kndC.obj",
                        help="path to the model", metavar="FILE")
    parser.add_argument("-n", "--numbers",
                        action="store_true", dest="printNumbers", default=False,
                        help="Print numbers on the cut edges")
    args = parser.parse_args()


    printNumbers = args.printNumbers

    # Import the mode
    mesh = om.read_trimesh(args.filename)

    fullUnfolded, unfoldedComponents = unfold(mesh)


    #Compute maxSize of the components
    # All components must be scaled to the same size as the largest component
    maxSize = 0
    for unfolding in unfoldedComponents:
        [xmin, ymin, boxSize] = findBoundingBox(unfolding[0])
        if boxSize > maxSize:
            maxSize = boxSize

    #Write SVG
    if printNumbers:
        basefilename = "unfoldingNumbers"
    else:
        basefilename = "unfolding"

    for i in range(len(unfoldedComponents)):
            writeSVG(basefilename + str(i) + ".svg", unfoldedComponents[i], maxSize, printNumbers)
            # writeSVG(basefilename + str(i) + ".svg", unfoldedComponents[i], maxSize, True)

    print("We wrote " + str(len(unfoldedComponents)) + " connected components.")

if __name__=='__main__':
    main()


