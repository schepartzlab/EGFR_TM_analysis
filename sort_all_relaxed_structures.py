"""
This script reads the output of a Rosetta relax run or runs, the script will walk through all folders in the directory to and apply relax analysis to those files
Created by Allison Walker, Schepartz Lab, Yale University
9/12/2017
"""
import  egf_tools_classification as egf_tools
import os
import math

egf_outfile = open("all_egf_angles.txt", 'w') #out file for cross angles for egf-type
tgfa_outfile = open("all_tgfa_angles.txt", 'w') #outfile for cross angles for tgfa-type
no_jm_contact_out = open("no_jm_contact.txt", 'w') #outfile for cross angles for structures that have
parallel_out = open("parallel.txt", 'w') #outfile for cross angles for structures that have parallel jms
not_anti_parallel_out = open("not_antiparallel.txt", "w") #outfile for cross angle for structures that have jms that are not parallel or antiparallel
ambiguous_out = open("ambigous_new.txt", "w") #outfile for cross angles for structures with ambiguous jms
other_jm_out = open("other_jm.txt", 'w') #outfile for cross angles for structures that do not fit any of the above classifications
cc_type_counts_out = open("cc_counts.txt", 'w') #counts of how many of each type of jm were present in relax results

directory_list = [x[0] for x in os.walk('.')] #get a list of all directories in the file

#lists of cross angles for each cross location and coiled coil type                  
cross_at_625_angles_egf = []
cross_at_629_angles_egf = []
cross_at_625_angles_tgfa = []
cross_at_629_angles_tgfa = []

#counts of different coiled coil type for each mp relax run
egf_dir_counts = {}
tgfa_dir_counts = {}
parallel_dir_counts = {}
not_antiparallel_dir_counts = {}
no_jm_contact_dir_counts = {}
other_jm_dir_counts = {}
ambiguous_dir_counts = {}

#analyze all structures and write all cross angles for each type of coiled coil to an output file
for directory in directory_list:
    egf_dir_counts[directory] = 0
    tgfa_dir_counts[directory] = 0
    parallel_dir_counts[directory] = 0
    not_antiparallel_dir_counts[directory] = 0
    no_jm_contact_dir_counts[directory] = 0
    other_jm_dir_counts[directory] = 0
    ambiguous_dir_counts[directory] = 0
    coords_dict = egf_tools.readAllFilesInDirectory(directory)
    for f in coords_dict:
        coords = coords_dict[f]
        jm_class = egf_tools.classifyJM(coords)
        
        (cross_loc, cross_angle, min_distance) = egf_tools.calcCrossLocationAndAngle(coords)
        if jm_class == "parallel":
            parallel_out.write(directory[2:len(directory)] + "," + f[1:len(f)] + "," + str(cross_loc) + "," + str(cross_angle)+ "," + str(min_distance) +"\n")
            parallel_dir_counts[directory] += 1
        if jm_class == "not antiparallel" or jm_class == "not antiparallel term dist":
            not_antiparallel_dir_counts[directory] += 1
            not_anti_parallel_out.write(directory[2:len(directory)] + "," + f[1:len(f)] + "," + str(cross_loc) + "," + str(cross_angle)+ "," + str(min_distance) +"\n")
        if jm_class == "no jm contact":
            no_jm_contact_dir_counts[directory] += 1
            no_jm_contact_out.write(directory[2:len(directory)] + "," + f[1:len(f)] + "," + str(cross_loc) + "," + str(cross_angle)+ "," + str(min_distance) + "\n")
        if jm_class == "ambiguous":
            ambiguous_dir_counts[directory] += 1
            ambiguous_out.write(directory[2:len(directory)] + "," + f[1:len(f)] + "," + str(cross_loc) + "," + str(cross_angle)+ "," + str(min_distance) +"\n")
        if jm_class == "other":
            other_jm_dir_counts[directory] += 1
            other_jm_out.write(directory[2:len(directory)] + "," + f[1:len(f)] + "," + str(cross_loc) + "," + str(cross_angle)+ "," + str(min_distance) +"\n")
        if jm_class != "egf" and jm_class != "tgfa":
            continue
        #if the structure does not make contact in the TM helix, do not include it in the counts for the egf and tgfa-type structures
        if cross_loc == -1:
            continue
        if jm_class == "egf":
            egf_dir_counts[directory] += 1
            egf_outfile.write(directory[2:len(directory)] + "," + f[1:len(f)] + "," + str(cross_loc) + "," + str(cross_angle)+ "," + str(min_distance)  + "\n")
            if cross_loc == 7:
                cross_at_625_angles_egf.append(cross_angle)
            if cross_loc == 11:
                cross_at_629_angles_egf.append(cross_angle)
                
        else:
            tgfa_dir_counts[directory] += 1
            tgfa_outfile.write(directory[2:len(directory)] + "," + f[1:len(f)] + "," + str(cross_loc) + "," + str(cross_angle)+ "," + str(min_distance)  + "\n")
            if cross_loc == 7:
                cross_at_625_angles_tgfa.append(cross_angle)
            if cross_loc == 11:
                cross_at_629_angles_tgfa.append(cross_angle)
                
#write the results of the coiled coil counts for each MPrelax simulation
for directory in egf_dir_counts:
    cc_type_counts_out.write(directory + ",")
    cc_type_counts_out.write(str(egf_dir_counts[directory]) + ",")
    cc_type_counts_out.write(str(tgfa_dir_counts[directory]) + ",")
    cc_type_counts_out.write(str(parallel_dir_counts[directory]) + ",")
    cc_type_counts_out.write(str(not_antiparallel_dir_counts[directory]) + ",")
    cc_type_counts_out.write(str(no_jm_contact_dir_counts[directory]) + ",")
    cc_type_counts_out.write(str(other_jm_dir_counts[directory]) + ",")
    cc_type_counts_out.write(str(ambiguous_dir_counts[directory]) + ",")
    cc_type_counts_out.write("\n")
cc_type_counts_out.close()
    
egf_outfile.close()
tgfa_outfile.close()
parallel_out.close()
not_anti_parallel_out.close()
no_jm_contact_out.close()
ambiguous_out.close()
other_jm_out.close()

#calculate average and standard deviation for egf and tgfa cross angles at each major cross location
egf_sum_625 = 0.0
egf_sum_629 = 0.0
tgfa_sum_625 = 0.0
tgfa_sum_629 = 0.0
for a in cross_at_625_angles_egf:
    egf_sum_625 += a
for a in cross_at_629_angles_egf:
    egf_sum_629 += a
for a in cross_at_625_angles_tgfa:
    tgfa_sum_625 += a
for a in cross_at_629_angles_tgfa:
    tgfa_sum_629 += a
print "sum"
print tgfa_sum_629
egf_sd_625 = 0.0
egf_sd_629 = 0.0
tgfa_sd_625 = 0.0
tgfa_sd_629 = 0.0


#calculate standard deviations
#print averages and standard deviation
if len(cross_at_625_angles_egf) > 0:
    avg = egf_sum_625/len(cross_at_625_angles_egf)
    print "avg egf 625: " + str(egf_sum_625/len(cross_at_625_angles_egf))
    for a in cross_at_625_angles_egf:
        egf_sd_625 += (a-avg)*(a-avg)
    print "sd: " + str(math.sqrt(egf_sd_625/len(cross_at_625_angles_egf)))
    print "len: " + str(len(cross_at_625_angles_egf))
else:
    print "no egf at 625"
    
if len(cross_at_629_angles_egf) > 0:
    avg = egf_sum_629/len(cross_at_629_angles_egf)
    print "avg egf 629: " + str(egf_sum_629/len(cross_at_629_angles_egf))
    for a in cross_at_629_angles_egf:
        egf_sd_629 += (a-avg)*(a-avg)
    print "sd: " + str(math.sqrt(egf_sd_629/len(cross_at_629_angles_egf)))
    print "len: " + str(len(cross_at_629_angles_egf))
else:
    print "no egf at 629"
    
if len(cross_at_625_angles_tgfa) > 0:
    avg = tgfa_sum_625/len(cross_at_625_angles_tgfa)
    print "avg tgfa 625: " + str(tgfa_sum_625/len(cross_at_625_angles_tgfa))
    for a in cross_at_625_angles_tgfa:
        tgfa_sd_625 += (a-avg)*(a-avg)
    print "sd: " + str(math.sqrt(tgfa_sd_625/len(cross_at_625_angles_tgfa)))
    print "len: " + str(len(cross_at_625_angles_tgfa))
else:
    print "no tgfa at 625"
    
if len(cross_at_629_angles_tgfa) > 0:
    avg = tgfa_sum_629/len(cross_at_629_angles_tgfa)
    print "avg tgfa 629: " + str(tgfa_sum_629/len(cross_at_629_angles_tgfa))
    for a in cross_at_629_angles_tgfa:
        tgfa_sd_629 += (a-avg)*(a-avg)
    print "sd: " + str(math.sqrt(tgfa_sd_629/len(cross_at_629_angles_tgfa)))
    print "len: " + str(len(cross_at_629_angles_tgfa))
else:
    print "no tgfa at 629"