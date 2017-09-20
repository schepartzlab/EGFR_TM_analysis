# EGFR_TM_analysis
Python scripts for analyzing rosetta runs on the TM-JM region of EGFR (taken from the Endres et al NMR structure PDBID: 2M20).
docking_analysis.py analyzes the results of a Rosetta MP docking_analysis
sort_all_relaxed_structures.py analyzes the results of a Rosetta MP relax simulation, sorting the structures by JM type and calculating the cross angles for these structures
egf_tools_classification.py contains methods used by sort_all_relaxed_structures.py 
