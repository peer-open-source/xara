# -----------------------------------------------------------------------------
# The following Python code is adapted from the work of Professor Terje Haukaas
# at the University of British Columbia in Vancouver, Canada. It is made freely
# available online at terje.civil.ubc.ca together with notes, examples, and
# additional Python code. Please be cautious when using this code; it may
# contain bugs and comes without warranty of any form.
# -----------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

class Edge:
    def __init__(self, xi, xj, thickness):
        self.vertices = (xi, xj)
        self.thickness = thickness

        # Compute derived quantities
        y1, z1 = xi
        y2, z2 = xj


        # Midpoint coordinates
        self.centroid = [0.5 * (y1 + y2),
                         0.5 * (z1 + z2)]

        self.length = np.sqrt((y2 - y1)**2 + (z2 - z1)**2)

        self.orientation = np.arctan((z2 - z1)/(y2 - y1)) \
                           if abs(y2 - y1) > 1e-8 else 0.5*np.pi

        self.area = self.thickness*self.length


def _omega_lookup(diagram, node):
    for i in range(len(diagram)):
        index = 0
        value = 0.0
        if int(diagram[i][0]) == int(node):
            return int(node), diagram[i][1]

    return index, value


# ------------------------------------------------------------------------
# IDENTIFY CELLS AND FLANGES
# ------------------------------------------------------------------------
def find_cells(edges, vertices):
    # Initialize lists of paths and cells in the cross-section
    pathList = [[int(edges[0, 0]),
                 int(edges[0, 1])]]
    naive_cells = []

    # Keep walking around in the cross-section, until all paths are covered
    stillPathsToWalk = True
    while stillPathsToWalk:

        foundMatch = False

        # Loop over all paths that have been walked
        for i in range(len(pathList)):

            # Extract path
            currentPath = pathList[i].copy()

            # Leave cells alone
            is_cell = False
            for k in range(1, len(currentPath)+1):
                num = currentPath.count(k)
                if num > 1:
                    is_cell = True

            if not is_cell:
                # Initialize before looping over edges
                leftCandidate = []
                rightCandidate = []

                # Loop over all edges and try to append the current path
                for j in range(edges.shape[0]):

                    # # Length of the current path
                    # lengthOfCurrentPath = len(currentPath)

                    # Get the left-most items in the path
                    leftEnd = currentPath[0]
                    nextToLeftEnd = currentPath[1]

                    # Get the right-most items in the path
                    rightEnd       = currentPath[ -1]
                    nextToRightEnd = currentPath[ -2]

                    # Detect link on the left-hand side, from first column of edges
                    if int(edges[j, 0]) == leftEnd:

                        # Make sure it is not the same edge
                        if int(edges[j, 1]) != nextToLeftEnd:

                            leftCandidate.append(int(edges[j, 1]))
                            foundMatch = True

                            # Is it a cell?
                            if int(edges[j, 1]) in currentPath:

                                cellCandidate = currentPath.copy()
                                cellCandidate.insert(0, int(edges[j, 1]))
                                for k in range(len(cellCandidate)):
                                    if cellCandidate[0] != int(edges[j, 1]):
                                        cellCandidate.pop(0)
                                naive_cells.append(cellCandidate)

                    # Detect link on the left-hand side, from second column of edges
                    if int(edges[j, 1]) == leftEnd:
                        # Make sure it is not the same edge
                        if int(edges[j, 0]) != nextToLeftEnd:

                            leftCandidate.append(int(edges[j, 0]))
                            foundMatch = True

                            # Is it a cell?
                            if int(edges[j, 0]) in currentPath:
                                cellCandidate = currentPath.copy()
                                cellCandidate.insert(0, int(edges[j, 0]))
                                for k in range(len(cellCandidate)):
                                    if cellCandidate[0] != int(edges[j, 0]):
                                        cellCandidate.pop(0)
                                naive_cells.append(cellCandidate)

                    # Detect link on the right-hand side, from first column of edges
                    if int(edges[j, 0]) == rightEnd:
                        # Make sure it is not the same edge
                        if int(edges[j, 1]) != nextToRightEnd:
                            rightCandidate.append(int(edges[j, 1]))
                            foundMatch = True
                            # Is it a cell?
                            if int(edges[j, 1]) in currentPath:
                                cellCandidate = currentPath.copy()
                                cellCandidate.append(int(edges[j, 1]))
                                for k in range(len(cellCandidate)):
                                    if cellCandidate[0] != int(edges[j, 1]):
                                        cellCandidate.pop(0)
                                naive_cells.append(cellCandidate)

                    # Detect link on the right-hand side, from second column of edges
                    if int(edges[j, 1]) == rightEnd:

                        # Make sure it is not the same edge
                        if int(edges[j, 0]) != nextToRightEnd:

                            rightCandidate.append(int(edges[j, 0]))
                            foundMatch = True

                            # Is it a cell?
                            if int(edges[j, 0]) in currentPath:
                                cellCandidate = currentPath.copy()
                                cellCandidate.append(int(edges[j, 0]))
                                for k in range(len(cellCandidate)):
                                    if cellCandidate[0] != int(edges[j, 0]):
                                        cellCandidate.pop(0)
                                naive_cells.append(cellCandidate)

                # See if we have any candidates for appending/branching the current path
                numLeftCandidates  = len(leftCandidate)
                numRightCandidates = len(rightCandidate)

                if numLeftCandidates == 1:
                    currentPath.insert(0, leftCandidate[0])

                elif numLeftCandidates > 1:
                    for k in range(1, numLeftCandidates):

                        newPath = currentPath.copy()

                        newPath.insert(0, leftCandidate[0])

                    currentPath.insert(0, leftCandidate[0])

                    # Insert the new path into the path list
                    pathList.append(newPath)

                if numRightCandidates == 1:
                    currentPath.append(rightCandidate[0])

                elif numRightCandidates > 1:
                    for k in range(1, numRightCandidates):
                        newPath = currentPath.copy()

                        newPath.append(rightCandidate[k])

                    currentPath.append(rightCandidate[0])

                    # Insert the new path into the path list
                    pathList.append(newPath)

                # Put the current path back into the path list
                pathList[i] = currentPath

        # Stop if there were no links anywhere
        if not foundMatch:
            stillPathsToWalk = False

    # Identify all flanges, i.e., segments not part of cells
    else:
        naive_flanges = []
        for i in range(len(pathList)):
            currentList = pathList[i]
            for j in range(len(pathList[i])-1):
                if len(naive_cells) == 0:
                    naive_flanges.append([currentList[j], currentList[j + 1]])
                else:
                    segmentAppearsInCell = False
                    for k in range(len(naive_cells)):
                        if currentList[j] in naive_cells[k] and currentList[j+1] in naive_cells[k]:
                            segmentAppearsInCell = True
                    if not segmentAppearsInCell:
                        naive_flanges.append([currentList[j], currentList[j+1]])

    # Create list of duplicate flanges
    flanges_to_remove = []
    for i in range(len(naive_flanges)):
        list1 = sorted(naive_flanges[i])
        for j in range(i+1, len(naive_flanges)):
            list2 = sorted(naive_flanges[j])
            if list1 == list2:
                flanges_to_remove.append(j)

    # Remove the duplicate flanges
    flanges_to_remove = np.unique(flanges_to_remove)
    flanges = []
    for i in range(len(naive_flanges)):
        if i not in flanges_to_remove:
            flanges.append(naive_flanges[i])

    # Create list of duplicate cells
    cells_to_remove = []
    for i in range(len(naive_cells)):
        list1 = np.unique(naive_cells[i])
        for j in range(i+1, len(naive_cells)):
            list2 = np.unique(naive_cells[j])
            if np.array_equal(list1, list2):
                cells_to_remove.append(j)

    # Remove the duplicate cells
    cells_to_remove = np.unique(cells_to_remove)
    uniqueCells = []
    for i in range(len(naive_cells)):
        if i not in cells_to_remove:
            uniqueCells.append(naive_cells[i])

    # Determine max-coordinates of the cells
    yMaxPerCell = []
    zMaxPerCell = []
    yMinPerCell = []
    zMinPerCell = []
    for cell in uniqueCells:
        yValues = []
        zValues = []
        # for j in range(len(cell)):
        #     yValues.append(vertices[int(cell[j])-1, 0])
        #     zValues.append(vertices[int(cell[j])-1, 1])
        for i in map(int, cell):
            yValues.append(vertices[i-1, 0])
            zValues.append(vertices[i-1, 1])
        yMaxPerCell.append(np.max(yValues))
        zMaxPerCell.append(np.max(zValues))
        yMinPerCell.append(np.min(yValues))
        zMinPerCell.append(np.min(zValues))

    # Create list of cells that encompass other cells
    cells_to_remove = []
    for i in range(len(uniqueCells)):
        for j in range(i+1, len(uniqueCells)):
            if yMaxPerCell[j] >= yMaxPerCell[i] and \
               yMinPerCell[j] <= yMinPerCell[i] and \
               zMaxPerCell[j] >= zMaxPerCell[i] and \
               zMinPerCell[j] <= zMinPerCell[i]:
                cells_to_remove.append(j)

    # Remove the duplicate cells
    cells_to_remove = np.unique(cells_to_remove)
    finalCells = []
    for i in range(len(uniqueCells)):
        if i not in cells_to_remove:
            finalCells.append(uniqueCells[i])

    # Output
    print('\n'"Paths:", pathList)
    print('\n'"Cells:", finalCells)
    print('\n'"Flanges:", flanges)
    return pathList, finalCells, flanges


def geometry(edges, vertices):
    # ------------------------------------------------------------------------
    # AREA AND CENTROID
    # ------------------------------------------------------------------------

    A  = 0.0
    yA = 0.0
    zA = 0.0
    areas = []
    centroids = []
    thetas = []
    lengths = []
    thicknesses = []
    for i,edge in enumerate(edges):

        # Coordinates of end 1
        end1 = edges[i, 0]
        y1 = vertices[int(end1)-1, 0]
        z1 = vertices[int(end1)-1, 1]

        # Coordinates of end 2
        end2 = edges[i, 1]
        y2 = vertices[int(end2)-1, 0]
        z2 = vertices[int(end2)-1, 1]

        # Midpoint coordinates
        yMidpoint = 0.5 * (y1 + y2)
        zMidpoint = 0.5 * (z1 + z2)
        centroids.append([yMidpoint, zMidpoint])

        # Length of segment
        delta_y = y2 - y1
        delta_z = z2 - z1
        L = np.sqrt(delta_y**2 + delta_z**2)
        lengths.append(L)

        # Angle relative to y-axis
        if delta_y != 0:
            theta = np.arctan(delta_z / delta_y)
        else:
            theta = 0.5 * np.pi

        thetas.append(theta)

        # Area of segment
        thickness = edges[i, 2]
        area = L * thickness
        areas.append(area)
        thicknesses.append(thickness)

        # Accumulate total area
        A += area

        # Accumulate A * y and A * z
        yA += yMidpoint * area
        zA += zMidpoint * area

    # Centroid location for the cross-section
    yCentroid = yA/A
    zCentroid = zA/A

    return A, areas, thicknesses, lengths, thetas, centroids, (yCentroid, zCentroid)

def moi(edges, lengths, thicknesses, thetas, centroid, centroids, areas):
    yCentroid, zCentroid = centroid

    # ------------------------------------------------------------------------
    # MOMENTS OF INERTIA and PRINCIPAL AXES ORIENTATION
    # ------------------------------------------------------------------------

    Iy = 0
    Iz = 0
    Iyz = 0
    for i in range(len(edges)):
        # Local moments of inertia
        Iyo  = lengths[i]**3 * thicknesses[i] * np.sin(thetas[i])**2 / 12.0
        Izo  = lengths[i]**3 * thicknesses[i] * np.cos(thetas[i])**2 / 12.0
        Iyzo = lengths[i]**3 * thicknesses[i] * np.cos(thetas[i]) * np.sin(thetas[i]) / 12.0

        # Add contributions from the parallel axis theorem (Steiner's sats)
        Iy += Iyo + ((centroids[i])[1] - zCentroid)** 2 * areas[i]
        Iz += Izo + ((centroids[i])[0] - yCentroid)** 2 * areas[i]
        Iyz += Iyzo + ((centroids[i])[0] - yCentroid) * ((centroids[i])[1] - zCentroid) * areas[i]

    # Apply caution with the principal axes rotation: what if Iy=Iz
    thetaPrincipal = 0.0
    if np.abs(Iyz) < 1e-10:
        thetaPrincipal = 0.0

    elif np.abs(Iy-Iz) < 1e-10:
        if np.abs(Iyz) < 1e-10:
            thetaPrincipal = 0.0
        else:
            thetaPrincipal = np.sign(Iyz) * 0.25 * np.pi
    else:
        thetaPrincipal = 0.5 * np.arctan(2.0*Iyz/(Iy - Iz))


    # Let counter-clockwise rotation of the system be positive
    thetaPrincipal *= -1

    # Transform the moments of inertia
    IyRotated =  0.5*(Iy + Iz) + 0.5*np.sqrt((Iz-Iy)**2 + 4.0* Iyz**2)
    IzRotated =  0.5*(Iy + Iz) - 0.5*np.sqrt((Iz-Iy)**2 + 4.0* Iyz**2)

    return thetaPrincipal, (Iy, Iz), (IyRotated, IzRotated)


def _venant_open(lengths, thicknesses):
    return sum(L * t**3 / 3.0 for L, t in zip(lengths, thicknesses))



def _venant_single(cells, flanges, centroid):
    J = 0
    cellArea = 0
    yCentroid, zCentroid = centroid

    # Loop over edges to calculate the area
    circleIntegral = 0
    for i in range(len(cells[0])-1):

        # Get the two points that are connected to this edge
        vertex1 = cells[0][i]
        vertex2 = cells[0][i+1]

        # Calculate the length of this edge
        y1 = vertices[vertex1-1, 0]
        z1 = vertices[vertex1-1, 1]
        y2 = vertices[vertex2-1, 0]
        z2 = vertices[vertex2-1, 1]
        delta_y = y2 - y1
        delta_z = z2 - z1
        length = np.sqrt(delta_y**2 + delta_z**2)

        # Get the thickness of this edge
        for i in range(edges.shape[0]):
            node1 = int(edges[i, 0])
            node2 = int(edges[i, 1])
            if node1 == vertex1 and node2 == vertex2:
                thickness = edges[i, 2]
            elif node1 == vertex2 and node2 == vertex1:
                thickness = edges[i, 2]

        # Add to the denominator integral
        circleIntegral += length / thickness

        # Add to the cell area, using the cross-product
        yCoordVector1 = y1 - yCentroid
        zCoordVector1 = z1 - zCentroid
        yCoordVector2 = y2 - yCentroid
        zCoordVector2 = z2 - zCentroid
        cellArea += 0.5 * np.abs(yCoordVector1*zCoordVector2-yCoordVector2*zCoordVector1)

    J = 4.0 * cellArea**2 / circleIntegral
    Jcell = J

    # Check if there are free flanges too
    for i in range(len(flanges)):

        # Get the two points that are connected to this edge
        vertex1 = flanges[i][0]
        vertex2 = flanges[i][1]

        # Calculate the length of this edge
        y1 = vertices[vertex1-1, 0]
        z1 = vertices[vertex1-1, 1]
        y2 = vertices[vertex2-1, 0]
        z2 = vertices[vertex2-1, 1]
        delta_y = y2 - y1
        delta_z = z2 - z1
        length = np.sqrt(delta_y**2 + delta_z**2)

        # Get the thickness of this edge
        for i in range(edges.shape[0]):
            node1 = int(edges[i, 0])
            node2 = int(edges[i, 1])
            if node1 == vertex1 and node2 == vertex2:
                thickness = edges[i, 2]
            elif node1 == vertex2 and node2 == vertex1:
                thickness = edges[i, 2]

        J += length * thickness**3 / 3.0

    # Contributions to J expressed as fractions
    cellFractionOfJ = Jcell/J
    return J, cellFractionOfJ, cellArea

def _rotate(vertices, edges, thetaPrincipal):
    rotated = vertices.copy()
    # ------------------------------------------------------------------------
    # ROTATE AXES IF Iyz HAS VALUE
    # ------------------------------------------------------------------------

    # Shift and rotate axis system to align with centroid and principal axes
    if np.abs(thetaPrincipal) > 1e-10:
        for i in range(len(vertices)):
            y = vertices[i-1, 0] - yCentroid
            z = vertices[i-1, 1] - zCentroid
            rotated[i - 1, 0] =  y * np.cos(thetaPrincipal) + z * np.sin(thetaPrincipal)
            rotated[i - 1, 1] = -y * np.sin(thetaPrincipal) + z * np.cos(thetaPrincipal)
        storeYcentroid = yCentroid
        storeZcentroid = zCentroid
        yCentroid = 0.0
        zCentroid = 0.0
    return rotated, edges

def _warp_open(finalCells, paths, vertices)->list:
    # Walk through the first path
    omega_trial = [[paths[0][0], 0.0]]
    for k in range(len(paths[0]) - 1):
        # Coordinates of end 1
        end1 = paths[0][k]
        y1 = vertices[int(end1) - 1, 0]
        z1 = vertices[int(end1) - 1, 1]

        # Coordinates of end 2
        end2 = paths[0][k + 1]
        y2 = vertices[int(end2) - 1, 0]
        z2 = vertices[int(end2) - 1, 1]

        # Cross-product
        cross = y2 * z1 - y1 * z2

        # Set omega values
        omega_trial.append([end2, (omega_trial[k])[1] + cross])

    # Walk through the rest of the paths
    for i in range(1, len(paths)):

        # Find the first node in this path that already has an omega value
        omega_found = False
        for k in range(len(paths[i])):

            # Check if there is an omega value at this node
            node = (paths[i])[k]
            omegaIndex, value = _omega_lookup(omega_trial, node)

            if omegaIndex != 0:
                omega_found = True
                locationInPath = k

                # Walk forward from first-found-node-with-omega-value
                for kk in range(locationInPath, len(paths[i])-1):

                    # Nothing to do if next node has omega value too
                    end2 = (paths[i])[kk+1]
                    omegaIndex, value = _omega_lookup(omega_trial, end2)
                    if omegaIndex == 0:

                        # Coordinates of end 1
                        end1 = (paths[i])[kk]
                        y1 = vertices[int(end1)-1, 0]
                        z1 = vertices[int(end1)-1, 1]

                        # Omega value at end 1
                        omegaIndex, value = _omega_lookup(omega_trial, end1)

                        # Coordinates of end 2
                        y2 = vertices[int(end2)-1, 0]
                        z2 = vertices[int(end2)-1, 1]

                        # Cross-product
                        cross = y2 * z1 - y1 * z2

                        # Amend omega diagram
                        omega_trial.append([end2, value + cross])

                # Walk backward from first-found-node-with-omega-value
                for kk in range(len(paths[i])-1, locationInPath-1, -1):

                    # Nothing to do if next node has omega value too
                    end2 = paths[i][kk-1]
                    omegaIndex, value = _omega_lookup(omega_trial, end2)
                    if omegaIndex == 0:

                        # Coordinates of end 1
                        end1 = (paths[i])[kk]
                        y1 = vertices[int(end1)-1, 0]
                        z1 = vertices[int(end1)-1, 1]

                        # Omega value at end 1
                        omegaIndex, value = _omega_lookup(omega_trial, end1)

                        # Coordinates of end 2
                        y2 = vertices[int(end2)-1, 0]
                        z2 = vertices[int(end2)-1, 1]

                        # Cross-product
                        cross = y2 * z1 - y1 * z2

                        # Amend omega diagram
                        omega_trial.append([end2, value + cross])

        if omega_found == False:
            raise ValueError("Found a path without any omega values")
    return omega_trial


def _warp_single(cells, flanges, vertices, centroid):
    # Figure out if the path around the cell is clockwise
    J, cellFractionOfJ, cellArea = _venant_single(cells, flanges, centroid)

    clockwise = False
    ySum = 0
    zSum = 0
    for i in range(len(cells[0])-1):
        node = cells[0][i]
        ySum += vertices[node-1, 0]
        zSum += vertices[node-1, 1]
    yCentre = ySum/float(len(cells[0])-1)
    zCentre = zSum/float(len(cells[0])-1)
    y1 = vertices[cells[0][0]-1, 0]
    y2 = vertices[cells[0][1]-1, 0]
    z1 = vertices[cells[0][0]-1, 1]
    z2 = vertices[cells[0][1]-1, 1]
    cross = (y1-yCentre) * (z2-zCentre) - (y2-yCentre) * (z1-zCentre)
    if cross < 0:
        clockwise = True

    # Walk around the cell
    omega_trial = [[cells[0][0], 0.0]]
    for i in range(len(cells[0])-2):

        # Length of this part
        end1 = cells[0][i]
        end2 = cells[0][i+1]

        y1 = vertices[end1-1, 0]
        z1 = vertices[end1-1, 1]
        y2 = vertices[end2-1, 0]
        z2 = vertices[end2-1, 1]
        delta_y = y2 - y1
        delta_z = z2 - z1
        length = np.sqrt( delta_y**2 + delta_z**2)

        # Thickness of this part
        for k in range(edges.shape[0]):
            node1 = int(edges[k, 0])
            node2 = int(edges[k, 1])

            if node1 == end1 and node2 == end2:
                thickness = edges[k, 2]
            elif node1 == end2 and node2 == end1:
                thickness = edges[k, 2]

        # Cross-product
        cross = y1*z2 - y2*z1

        # Correction for cell walls
        theSign = -1
        if clockwise:
            theSign = 1
        correction = theSign * cellFractionOfJ * J * length / (2.0 * thickness * cellArea)

        # Set omega values
        omega_trial.append([end2, omega_trial[i][1] + cross + correction])

    # Now address free flanges obeying kinematic compatibility with the cell they are connected to
    for i in range(len(flanges)):

        # Get the two points that are connected to this flange
        end1 = flanges[i][0]
        end2 = flanges[i][1]
        y1 = vertices[end1-1, 0]
        z1 = vertices[end1-1, 1]
        y2 = vertices[end2-1, 0]
        z2 = vertices[end2-1, 1]

        # Find the cell node that is connected to this flange
        for j in range(len(cells[0])-1):

            node = cells[0][j]
            yNode = vertices[node - 1, 0]
            zNode = vertices[node - 1, 1]
            if end1 == node:

                # Must now find the cell-node next to "node" that is in the same plane as End 1
                forwardNode = 0
                backNode = 0
                if j == 0:
                    backNode = cells[0][len(cells[0])-2]
                    forwardNode = cells[0][j + 1]
                else:
                    forwardNode = cells[0][j+1]
                    backNode = cells[0][j-1]
                yBack = vertices[backNode - 1, 0]
                zBack = vertices[backNode - 1, 1]
                yForward = vertices[forwardNode - 1, 0]
                zForward = vertices[forwardNode - 1, 1]

                if zBack == z1 and z1 == z2:
                    _, nodeValue = _omega_lookup(omega_trial, node)
                    _, backValue = _omega_lookup(omega_trial, backNode)
                    change = (nodeValue-backValue)/(yNode-yBack)*(y2-y1)
                    omega_trial.append([end2, nodeValue + change])
                    break

                elif yBack == y1 and y1 == y2:
                    _, nodeValue = _omega_lookup(omega_trial, node)
                    _, backValue = _omega_lookup(omega_trial, backNode)
                    change = (nodeValue-backValue)/(zNode-zBack)*(z2-z1)
                    omega_trial.append([end2, nodeValue + change])
                    break

                elif zForward == z1 and z1 ==  z2:
                    _, nodeValue = _omega_lookup(omega_trial, node)
                    _, forwardValue = _omega_lookup(omega_trial, forwardNode)
                    change = (nodeValue-forwardValue)/(yNode-yForward)*(y2-y1)
                    omega_trial.append([end2, nodeValue + change])
                    break

                elif yForward == y1 and y1 == y2:
                    _, nodeValue = _omega_lookup(omega_trial, node)
                    _, forwardValue = _omega_lookup(omega_trial, forwardNode)
                    change = (nodeValue-forwardValue)/(zNode-zForward)*(z2-z1)
                    omega_trial.append([end2, nodeValue + change])
                    break

            elif end2 == node:

                # Now find the cell-node next to "node" that is in the same plane as End 2
                forwardNode = 0
                backNode = 0
                if j == 0:
                    backNode = cells[0][len(cells[0]) - 2]
                    forwardNode = cells[0][j + 1]
                else:
                    forwardNode = cells[0][j + 1]
                    backNode = cells[0][j - 1]
                yBack = vertices[backNode - 1, 0]
                zBack = vertices[backNode - 1, 1]
                yForward = vertices[forwardNode - 1, 0]
                zForward = vertices[forwardNode - 1, 1]

                if zBack == z1 and z1 == z2:
                    _, nodeValue = _omega_lookup(omega_trial, node)
                    _, backValue = _omega_lookup(omega_trial, backNode)
                    change = (nodeValue - backValue) / (yNode - yBack) * (y1 - y2)
                    omega_trial.append([end1, nodeValue + change])
                    break

                elif yBack == y1 and y1 == y2:
                    _, nodeValue = _omega_lookup(omega_trial, node)
                    _, backValue = _omega_lookup(omega_trial, backNode)
                    change = (nodeValue - backValue) / (zNode - zBack) * (z1 - z2)
                    omega_trial.append([end1, nodeValue + change])
                    break

                elif zForward == z1 and z1 == z2:
                    _, nodeValue = _omega_lookup(omega_trial, node)
                    _, forwardValue = _omega_lookup(omega_trial, forwardNode)
                    change = (nodeValue - forwardValue) / (yNode - yForward) * (y1 - y2)
                    omega_trial.append([end1, nodeValue + change])
                    break

                elif yForward == y1 and y1 == y2:
                    _, nodeValue = _omega_lookup(omega_trial, node)
                    _, forwardValue = _omega_lookup(omega_trial, forwardNode)
                    change = (nodeValue - forwardValue) / (zNode - zForward) * (z1 - z2)
                    omega_trial.append([end1, nodeValue + change])
                    break

    # Flip the sign of the omega diagram if the walk was counterclockwise
    for i in range(len(omega_trial)):
        omega_trial[i][1] = theSign * omega_trial[i][1]

    return omega_trial


def warp(edges, vertices):
    # ------------------------------------------------------------------------
    # WARPING DIAGRAM PLUS ysc, zsc, Cw
    # ------------------------------------------------------------------------
    paths, cells, flanges = find_cells(edges, vertices)
    A, areas, thicknesses, lengths, thetas, centroids, centroid = geometry(edges, vertices)
    (yCentroid, zCentroid) = centroid
    thetaPrincipal, (Iy, Iz), (IyRotated, IzRotated) = moi(edges, lengths, thicknesses, thetas, centroid, centroids, areas)
    vertices, edges = _rotate(vertices, edges, thetaPrincipal)

    if len(cells) > 1:
        print('\n'"Multi-cell cross-sections not considered in warping torsion calculations")

    # Open cross-section
    elif len(cells) == 0:
        omega_trial = _warp_open(cells, paths, vertices)

    # Single-cell cross-sections
    elif len(cells) == 1:
        omega_trial = _warp_single(cells, flanges, vertices, centroid)

    # Calculate integrals needed to find final omega diagram
    intOmega  = 0.0
    intyOmega = 0.0
    intzOmega = 0.0
    for i in range(len(edges)):
        # End 1 of this edge
        end1 = edges[i, 0]
        y1 = vertices[int(end1) - 1, 0] - yCentroid
        z1 = vertices[int(end1) - 1, 1] - zCentroid

        # End 2 of this edge
        end2 = edges[i, 1]
        y2 = vertices[int(end2) - 1, 0] - yCentroid
        z2 = vertices[int(end2) - 1, 1] - zCentroid

        # Omega values for this edge
        _, w1 = _omega_lookup(omega_trial, end1)
        _, w2 = _omega_lookup(omega_trial, end2)

        # Pick up stored length and thickness
        dA = lengths[i]*thicknesses[i]

        # Integrals
        intOmega  += 0.5*(w1 + w2)*dA
        intyOmega += dA/6.0*(y1*(2.0*w1 + w2) + y2*(2.0*w2 + w1))
        intzOmega += dA/6.0*(z1*(2.0*w1 + w2) + z2*(2.0*w2 + w1))

    # Shear centre coordinates and normalizing constant
    if np.abs(thetaPrincipal) > 1e-10:
        ysc = -intzOmega/IyRotated
        zsc =  intyOmega/IzRotated
    else:
        ysc = -intzOmega / Iy
        zsc =  intyOmega / Iz

    # Finalize omega diagram
    C  = -intOmega/A
    omegaFinal = []
    minWarp = 0
    maxWarp = 0

    for i in range(len(omega_trial)):
        node  = omega_trial[i][0]
        yHere = vertices[node-1, 0] - yCentroid
        zHere = vertices[node-1, 1] - zCentroid
        trialValue = omega_trial[i][1]
        finalValue = trialValue + C + ysc * zHere - zsc * yHere

        omegaFinal.append([node, finalValue])
        if finalValue < minWarp:
            minWarp = finalValue
        if finalValue > maxWarp:
            maxWarp = finalValue


    return omegaFinal, (zsc, ysc), (minWarp, maxWarp)


def warping_constant(edges, vertices, omega):
    A, areas, __, lengths, *_ = geometry(edges, vertices)
    Cw = 0.0
    # Cross-section constant for warping torsion
    for edge, length in zip(edges, lengths):
        # Nodes of this segment
        i,j,thick = edge
        # Area of this segment
        dA = length*thick

        # Omega values at segment nodes
        _, w1 = _omega_lookup(omega, i)
        _, w2 = _omega_lookup(omega, j)

        # Accumulate the integral
        if np.sign(w1) == np.sign(w2):
            Cw += dA/3.0*(w1**2 + w2**2 + w1*w2)
        else:
            Cw += dA/3.0*(w1**2 + w2**2 - np.abs(w1*w2))
    return Cw


    # ------------------------------------------------------------------------
    # ROTATE AXES BACK
    # ------------------------------------------------------------------------

    if np.abs(thetaPrincipal) > 1e-10:
        yCentroid = storeYcentroid
        zCentroid = storeZcentroid
        y = ysc
        z = zsc
        ysc = y * np.cos(thetaPrincipal) - z * np.sin(thetaPrincipal) + yCentroid
        zsc = y * np.sin(thetaPrincipal) + z * np.cos(thetaPrincipal) + zCentroid
        for i in range(len(vertices)):
            y = vertices[i-1, 0]
            z = vertices[i-1, 1]
            vertices[i - 1, 0] = y * np.cos(thetaPrincipal) - z * np.sin(thetaPrincipal) + yCentroid
            vertices[i - 1, 1] = y * np.sin(thetaPrincipal) + z * np.cos(thetaPrincipal) + zCentroid

def plot(edges, vertices, omegaFinal, Cw, sc, minmax):
    A, areas, thicknesses, lengths, thetas, centroids, centroid = geometry(edges, vertices)
    (yCentroid, zCentroid) = centroid
    thetaPrincipal, (Iy, Iz), (IyRotated, IzRotated) = \
            moi(edges, lengths, thicknesses, thetas, centroid, centroids, areas)
    zsc, ysc = sc
    (minWarp, maxWarp) = minmax

    for i in range(len(edges)):

        # End 1 of this edge
        end1 = edges[i, 0]
        y1 = vertices[int(end1)-1, 0]
        z1 = vertices[int(end1)-1, 1]

        # End 2 of this edge
        end2 = edges[i, 1]
        y2 = vertices[int(end2)-1, 0]
        z2 = vertices[int(end2)-1, 1]

        # Omega values for this edge
        if Cw != 0:
            _, w1 = _omega_lookup(omegaFinal, end1)
            _, w2 = _omega_lookup(omegaFinal, end2)

        # Discretize and plot the line
        segments = 10
        for j in range(segments):
            y1Here = y1 + j/float(segments)*(y2-y1)
            z1Here = z1 + j/float(segments)*(z2-z1)
            y2Here = y1 + (j+1)/float(segments)*(y2-y1)
            z2Here = z1 + (j+1)/float(segments)*(z2-z1)

            if Cw > 1e-5:
                warpValueHere = w1 + j/float(segments) * (w2 - w1)
                colorCode = (warpValueHere - minWarp) / (maxWarp - minWarp)
                plt.plot(np.array([y1Here, y2Here]), np.array([z1Here, z2Here]), color=(colorCode, 0, 1 - colorCode), linewidth=6.0, zorder=1)
            else:
                colorCode = 0.0
                plt.plot(np.array([y1Here, y2Here]), np.array([z1Here, z2Here]), color=(0.52, 0.52, 0.52), linewidth=6.0, zorder=1)

    # Create the plot
    max = np.amax(np.amax((vertices)))
    min = np.amin(np.amin(vertices))
    yMin = min - (max-min)*0.4
    border = (max - min) * 0.1
    plt.plot(np.array([yCentroid, max+border]), np.array([zCentroid, zCentroid+(np.abs(max)+border-yCentroid)*np.tan(thetaPrincipal)]), 'k-',linewidth=1.0)
    plt.plot(np.array([yCentroid, yMin-border]), np.array([zCentroid, zCentroid-(np.abs(yMin)+border+yCentroid)*np.tan(thetaPrincipal)]), 'k-',linewidth=1.0)
    plt.plot(np.array([yCentroid, yCentroid-(np.abs(max)+border-zCentroid)*np.tan(thetaPrincipal)]), np.array([zCentroid, max+border]), 'k-',linewidth=1.0)
    plt.plot(np.array([yCentroid, yCentroid+(np.abs(min)+border+zCentroid)*np.tan(thetaPrincipal)]), np.array([zCentroid, min-border]), 'k-',linewidth=1.0)
    plt.scatter([yCentroid], [zCentroid], s=100, c='k', marker='x', zorder=2)
    plt.scatter([ysc], [zsc], s=100, c='k', marker='x', zorder=2)
    plt.axis([yMin-border, max+border, min-border, max+border])
    plt.title("Centroid + Shear Centre + Warping")
    plt.show()


if __name__ == "__main__":

# ------------------------------------------------------------------------
# INPUT (The last one is used by the code)
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# S-SHAPE
    t = 10.0
    b = 200.0
    h = 200.0
    vertices = np.array([(-b, 0),
                         (0, 0),
                         (0, h),
                         (b, h)])
    edges = np.array([(1, 2, t),
                      (2, 3, t),
                      (3, 4, t)])

# ------------------------------------------------------------------------
# ONE CELL
    b = 500.0
    h = 200.0
    ttop = 10.0
    tbottom = 10.0
    tleft = 15.0
    tright = 20.0
    vertices = np.array([(0.0, 0.0),
                         (b,   0.0),
                         (0.0, h),
                         (b,   h)])
    edges = np.array([(1, 2, tbottom),
                      (3, 4, ttop),
                      (1, 3, tleft),
                      (2, 4, tright)])

# ------------------------------------------------------------------------
# LAMBDA-SHAPE
    t = 10.0
    b = 100.0
    h = 200.0
    vertices = np.array([(0, 0),
                         (0, h),
                         (b, h)])
    edges = np.array([(1, 2, t),
                      (2, 3, t)])

# ------------------------------------------------------------------------
# CELL WITH FLANGES
    b = 200.0
    h = 150.0
    a = 80.0
    t = 10.0
    vertices = np.array([(-a,  h),
                         ( 0,  h),
                         ( b,  h),
                         (b+a, h),
                         ( 0,  0),
                         ( b,  0)])

    edges = np.array([(1, 2, t),
                      (2, 3, t),
                      (3, 4, t),
                      (5, 2, t),
                      (6, 3, t),
                      (5, 6, t)])

    omega, sc, minmax = warp(edges, vertices)
    Cw = warping_constant(edges, vertices, omega)
    plot(edges, vertices, omega, Cw, sc, minmax)
    print(f"{Cw = }")


# ------------------------------------------------------------------------
# I-BEAM
    ttop = 15.0
    tbottom = 15.0
    tweb = 15.0
    btop = 100.0
    bbottom = 150.0
    h = 250.0
    vertices = np.array([(-btop/2, h),
                         (0.0,  h),
                         (btop/2,  h),
                         (-bbottom/2, 0.0),
                         (0.0,  0.0),
                         (bbottom/2,  0.0)])

    edges = np.array([(1, 2, ttop),
                      (2, 3, ttop),
                      (2, 5, tweb),
                      (4, 5, tbottom),
                      (5, 6, tbottom)])

# ------------------------------------------------------------------------
# L-SHAPE
    t = 10.0
    b = 100.0
    h = 200.0
    vertices = np.array([(0, 0),
                         (0, h),
                         (b, 0)])
    edges = np.array([(1, 2, t),
                      (1, 3, t)])

# ------------------------------------------------------------------------
# CHANNEL PROFILE
    t = 20.0
    b = 80.0
    h = 200.0
    vertices = np.array([(b, 0),
                         (0, 0),
                         (0, h),
                         (b, h)])

    edges = np.array([(1, 2, t),
                      (2, 3, t),
                      (3, 4, t)])

# ------------------------------------------------------------------------

    omega, sc, minmax = warp(edges, vertices)
    Cw = warping_constant(edges, vertices, omega)
    plot(edges, vertices, omega, Cw, sc, minmax)
    print(f"{Cw = }")

#   # ------------------------------------------------------------------------
#   # PRINT RESULTS
#   # ------------------------------------------------------------------------
#   s  = ''
#   s += f'\n'"Area: {A}'
#   s += ('\n'"Centroid y-coordinate: %.2f" % yCentroid)
#   s += ('\n'"Centroid z-coordinate: %.2f" % zCentroid)
#   s += ('\n'"Moment of inertia Iy: %.2f" % Iy)
#   s += ('\n'"Moment of inertia Iz: %.2f" % Iz)
#   if np.abs(thetaPrincipal) > 1e-10:
#       s += ('\n'"Principal axis rotation (degrees): %.2f" % (thetaPrincipal/np.pi*180.0))
#       s += ('\n'"Product of inertia Iyz: %.2f" % Iyz)
#       s += ('\n'"New principal moment of inertia Iy: %.2f" % IyRotated)
#       s += ('\n'"New principal moment of inertia Iz: %.2f" % IzRotated)
#       s += ('\n'"Shear centre y-coordinate in original coordinate system: %.2f" % ysc)
#       s += ('\n'"Shear centre z-coordinate in original coordinate system: %.2f" % zsc)
#   else:
#       s += ('\n'"Shear centre y-coordinate: %.2f" % ysc)
#       s += ('\n'"Shear centre z-coordinate: %.2f" % zsc)
#   s += ('\n'"St. Venant torsion constant J: %.2f" % J)
#   s += ('\n'"Warping torsion constant Cw: %.2f" % Cw)
