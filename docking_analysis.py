"""
This script reads the output of a Rosetta dock run (file name prefix given by pdb_filename)
Created by Allison Walker, Schepartz Lab, Yale University
9/12/2017
"""
import math

pdb_filename = "" #Set this variable to the file name prefix for the Rosetta dock output
n_decoys = 1000 #number of decoys produced by the docking run
score_file = open("score.sc", 'r') #name of the file that has the docked structures score
directory_name = "" #name of the directory where the docking run is located
outfile = open(directory_name +"cross_and_angle_cluster.txt", 'w') #outfile to store the results of the docking analysis

#read the scores from the score file
scores = {}
i =0
for line in score_file:
    if "SEQUENCE" in line or "total_score" in line:
        continue
    scores[i] = float(line[10:18])
    i += 1
score_file.close()

#read the decoy files
filenames = []
for i in range(1, n_decoys + 1):
    if i >= 1000:
        num_str = str(i)
        filenames.append(pdb_filename + num_str + ".pdb")
    elif i >= 100:
        num_str = "0" + str(i)
        filenames.append(pdb_filename + num_str + ".pdb")
    elif i >= 10:
        num_str = "00" + str(i)
        filenames.append(pdb_filename + num_str + ".pdb")
    else:
        num_str = "000" + str(i)
        filenames.append(pdb_filename + num_str + ".pdb")

all_coords = {}
i = 0
for f in filenames:
    in_file = open(f, 'r')
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
    all_coords[i] =(ca_coords_A, ca_coords_B)
    i+= 1
    in_file.close()

#find the cross location
cross_location = {}
distinct_cross_location = []
for i in all_coords:
    A = all_coords[i][0]
    B = all_coords[i][1]
    min_distance = 1000
    min_distance_index = 0
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
        cross_location[i] = -1
        continue
    else:
        cross_location[i] = min_distance_index
    if min_distance_index not in distinct_cross_location:
        distinct_cross_location.append(min_distance_index)


#calculate cross angles
cross_angles = {}
for i in all_coords:
    A = all_coords[i][0]
    B = all_coords[i][1]
    cross_at = cross_location[i]
    if cross_location == -1:
        cross_angles[i] == -1
        continue
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
    
    angle = math.acos(dot_product/(A_norm*B_norm))*180/3.14
    cross_angles[i] = angle



#group by cross location and then by angles
groups = {}
for cl in distinct_cross_location:
    list_members = []
    for i in cross_location:
        if cross_location[i] == cl:
            list_members.append(i)
    if len(list_members) < 5:
        continue
    groups[cl] = {}
    min_angle = 1000
    max_angle = 0
    for i in list_members:
        angle = cross_angles[i]
        if angle < min_angle:
            min_angle = angle
        if angle > max_angle:
            max_angle = angle
    #bin the angles
    num_bins = int(math.ceil(math.sqrt(len(list_members))))
    for i in range(0, num_bins):
        bin_min = i*(max_angle - min_angle)/num_bins + min_angle
        bin_max = (i+1)*(max_angle - min_angle)/num_bins + min_angle
        groups[cl][(bin_min, bin_max)] = []
        for j in list_members:
            angle= cross_angles[j]
            if angle < bin_max and angle >= bin_min:
                groups[cl][(bin_min, bin_max)].append((j, scores[j], angle))
            if angle == max_angle and i == num_bins -1:
                groups[cl][(bin_min, bin_max)].append((j, scores[j], angle))

#write output file
#rank groups by energy, find representitive minimum energy structure for each cross angle
for cl in groups:
    cross_angle_vs_score_file = open(directory_name +"cross_angle_vs_score_cross_at_" + str(cl) +".txt", 'w')
    outfile.write("crossing at amino acid: " + str(cl) + "\n")
    for angles in groups[cl]:
        outfile.write("between angles " + str(angles[0]) +" and " + str(angles[1]) + "\n")
        sorted_by_score = sorted(groups[cl][angles], key=lambda tup: tup[1])
        for s in sorted_by_score:
            outfile.write("decoy number: " + str(s[0]+1) + " score: " + str(s[1]) + " angle: " + str(s[2]) + "\n")
            cross_angle_vs_score_file.write(str(s[2]) + ", " +str(s[1]) + "\n")
        outfile.write("\n")
    cross_angle_vs_score_file.close()
    outfile.write("\n")
outfile.close()