"""
This python file contains methods used by sort_all_relaxed_structres.py to classify coiled coil types in EGFR TM-JM structures
Created by Allison Walker, Schepartz Lab, Yale University
9/12/2017
"""
import math
import os

"""
This methods takes the filename of a PDB file and returns a list of coordinates of the Ca atoms
arguments: filename - the relative path for the PDB file
returns: a list of coordinates for the A and B chains in the PDB file
"""
def readPDBFile(filename):
    in_file = open(filename, 'r')
    ca_coords_A = []
    ca_coords_B = []
    for line in in_file:
        if "CA" in line and "ATOM" in line:
            x = float(line[31:38])
            y = float(line[39:46])
            z = float(line[47:56])
            if " A " in line:
                ca_coords_A.append((x, y, z))
            else:
                ca_coords_B.append((x, y, z))
    in_file.close()
    return (ca_coords_A, ca_coords_B)
 
"""
This methods takes a directory and returns the Ca coordinates for all files with the pdb extension
arguments: directory - the relative path for the directory
returns: a list of a list of coordinates
"""    
def readAllFilesInDirectory(directory):
    coordinates = {}
    for filename in os.listdir(directory):
        if ".pdb" in filename and filename[len(filename)-10] == "_":
            (A, B) = readPDBFile(directory + "/" + filename)
            if len(A) == 0:
                continue
            coordinates[filename[len(filename)-10:len(filename)-4]] = (A, B)
    return  coordinates

"""
This methods takes list of chain coordinates and returns the cross location and angle of those chains
arguments: coords = a 2D array where the first entry is a list of coordinates for one chain and th second entry is a list of coordinates for the second chain
returns: a tuple with three entries
"""       
def calcCrossLocationAndAngle(coords):
    A = coords[0]
    B = coords[1]
    min_distance = 1000
    min_distance_index = 0
    #find the cross location
    for j in range(0,len(A)):
        xA = A[j][0]
        yA = A[j][1]
        zA = A[j][2]
        xB = B[j][0]
        yB = B[j][1]
        zB = B[j][2]
        distance = math.sqrt(math.pow(xA-xB, 2) + math.pow(yA-yB, 2) + math.pow(zA-zB, 2))
        if distance < min_distance and j < 30 and j > 2:
            min_distance = distance
            min_distance_index = j
    if min_distance > 13:
        return (-1, -1,-1)
    cross_at = min_distance_index
    #find the cross angle
    prevXA = A[cross_at][0]
    prevYA = A[cross_at][1]
    prevZA = A[cross_at][2]
    prevXB = B[cross_at][0]
    prevYB = B[cross_at][1]
    prevZB = B[cross_at][2]
    vectorAX = 0.0
    vectorAY = 0.0
    vectorAZ = 0.0
    vectorBX = 0.0
    vectorBY = 0.0
    vectorBZ = 0.0
    for j in range(cross_at +4, cross_at +5):
        vectorAX += A[j][0] - prevXA
        vectorAY += A[j][1] - prevYA
        vectorAZ += A[j][2] - prevZA
        vectorBX += B[j][0] - prevXB
        vectorBY += B[j][1] - prevYB
        vectorBZ += B[j][2] - prevZB
        prevXA =A[j][0]
        prevYA =A[j][1]
        prevZA =A[j][2]
        prevXB =B[j][0]
        prevYB =B[j][1]
        prevZB =B[j][2]
    dot_product = vectorAX*vectorBX + vectorAY*vectorBY + vectorAZ*vectorBZ
    A_norm = math.sqrt(vectorAX*vectorAX + vectorAY*vectorAY + vectorAZ*vectorAZ)
    B_norm = math.sqrt(vectorBX*vectorBX + vectorBY *vectorBY + vectorBZ * vectorBZ)
    angle = math.acos(dot_product/(A_norm*B_norm))*180/math.pi
    #print "here!"
    return (cross_at, angle, min_distance)

"""
This methods takes list of chain coordinates and returns the coiled-coil classification
arguments: coords = a 2D array where the first entry is a list of coordinates for one chain and th second entry is a list of coordinates for the second chain
returns: a string denoting the JM coiled coil type, possible values are:
    'parallel' - if the angle between vectors representing the two coiled coils are less than 20 degrees
    'not antiparallel' - returned if the structure is between parallel and anti parallel (if the angle between them is greater than 20 degrees but less than 90 degrees)
    'not antiparallel term dist' - returned if the coiled coils are antiparallel because they are not lined up with each other (if the c-terminus of one coil is closer to the c-terminus of the other coil, than it is to the n-terminus of the other coil)
    'no jm contact' - assigned if the total distance between a pair of coiled-coil faces exceeds 40 Angstroms
    'egf' - assigned if the face1 to face1 distance is smaller than any other face-face distances
    'tgfa' - assigned if face3 to face3 distance is smaller than any other face-face distances
    'other' - assigned if the coiled coil is an antiparallel coiled coil that is not egf or tgfa 
"""           
def classifyJM(coords):
    A = coords[0]
    B = coords[1]
    #residues that define the faces
    face1 = [37, 41, 45] #cch-1/egf
    face2 = [36, 40, 44] #other
    face3 = [35, 39, 43] #cch5/tgfa
    face3_2 = [39, 43, 47] #cch5/tgfa (One turn down the JM from face3)
    face4 = [38, 42, 46] #cch9/both
    faces = [face1, face2, face3, face3_2, face4]
    #residues that define the c and n termini of the JM coil
    c_term = 46
    n_term = 35
    c_term_c_term =math.sqrt(math.pow(A[c_term][0] - B[c_term][0], 2) + math.pow(A[c_term][1] - B[c_term][1], 2) + math.pow(A[c_term][2] - B[c_term][2], 2))
    n_term_c_term = math.sqrt(math.pow(A[n_term][0] - B[c_term][0], 2) + math.pow(A[n_term][1] - B[c_term][1], 2) + math.pow(A[n_term][2] - B[c_term][2], 2))
    
    A_n_to_c = (A[n_term][0] - A[c_term][0], A[n_term][1] - A[c_term][1], A[n_term][2] - A[c_term][2])
    B_n_to_c = (B[n_term][0] - B[c_term][0], B[n_term][1] - B[c_term][1], B[n_term][2] - B[c_term][2])
    dot_prod = A_n_to_c[0]*B_n_to_c[0] + A_n_to_c[1]*B_n_to_c[1] + A_n_to_c[2]*B_n_to_c[2]
    A_norm = math.sqrt(A_n_to_c[0]*A_n_to_c[0] + A_n_to_c[1]*A_n_to_c[1] + A_n_to_c[2]*A_n_to_c[2])
    B_norm = math.sqrt(B_n_to_c[0]*B_n_to_c[0] + B_n_to_c[1]*B_n_to_c[1] + B_n_to_c[2]*B_n_to_c[2])
    angle = math.acos(dot_prod/(A_norm*B_norm))*180/math.pi
    if angle < 20:
        return "parallel"
    if angle < 90:
        return "not antiparallel"
    if c_term_c_term < n_term_c_term:
        return "not antiparallel term dist"
    distances_for_current_struct = []
    min_dist = 1000
    closest_pair = ()
    distances_for_current_struct = []
    for i in range (0, 5):
        faceA = faces[i]
        for j in range(0, 5):
            faceB = faces[j]
            distance = math.sqrt(math.pow(A[faceA[0]][0] - B[faceB[2]][0], 2) + math.pow(A[faceA[0]][1] - B[faceB[2]][1], 2) + math.pow(A[faceA[0]][2] - B[faceB[2]][2], 2))
            distance += math.sqrt(math.pow(A[faceA[1]][0] - B[faceB[1]][0], 2) + math.pow(A[faceA[1]][1] - B[faceB[1]][1], 2) + math.pow(A[faceA[1]][2] - B[faceB[1]][2], 2))
            distance += math.sqrt(math.pow(A[faceA[2]][0] - B[faceB[0]][0], 2) + math.pow(A[faceA[2]][1] - B[faceB[0]][1], 2) + math.pow(A[faceA[2]][2] - B[faceB[0]][2], 2))
            if distance < min_dist:
                closest_pair = (i+1, j+1)
                min_dist = distance
            distances_for_current_struct.append(distance)
    if min_dist > 40:
        return "no jm contact"
    if closest_pair == (1, 1):
        return "egf"
    elif closest_pair == (4,4) or closest_pair == (3, 4) or closest_pair == (4, 3):
        return "tgfa"
    else:
        return "other"
    