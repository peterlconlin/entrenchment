#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 10:28:20 2021

@author: peterconlin

Reading file line by line (avoids loading the whole file)
Before using this:
I have to know how columns are separated and how lines (rows) are separated
This case is a 2 column separated by ´\t´ and each line is separated by "\n"

"""

import os
import glob

my_files = [f for f in glob.glob("/Users/peterconlin/Desktop/Entrenchment_data_05MAY21/details_final/*unicell_detail.dat")]

for input_file in my_files:
    
    # Open input file
    file_in = input_file
    pf=open(file_in,'r') # this opens a pointer in memory to a file
    filename = os.path.basename(file_in)
    
    # Open output file
    file_out = open("%s_output_v2.txt" % file_in, "w")
    file_out.write("strategy_rep\t")
    file_out.write("timepoint\t")
    file_out.write("genome_location\t")
    file_out.write("mutation\t")
    file_out.write("iteration\n")
    
    # Read header, set counter
    pf.readline() # this reads the header
    count = 0
    
    # Create object to track last line saved to output
    multicelled=False
    replicated=False
    line_out=None
    
    for line in pf: # line is any string between\n
        count += 1
        temp = [ line.strip().split('"')[0].split(' ')[i] for i in [0, 1, 2, 3, 4] ]
        #temp0 = line.strip().split('"')[0]
        #temp = [ temp0.split(' ')[i] for i in [0, 1, 2] ]     #this strips out the '\n' from the line string
                                                                    #and splits the string into 2 strings separated by '\t
                                                                    #and selects the relevant elements of the list
        
        if count==1:
            # Start tracking a new mutant and continue to the next line
            multicelled=False
            replicated=False
            line_past=temp # Update line_past
            line_out=temp[0:4]
            #print(line_past)
            
        else:
            # Check for new mutant
            # If mutation doesn't match: 
            if temp[0:4] != line_past[0:4]:
                
                #First, save data about previous mutant
                if line_out != None:
                    
                    if multicelled == False and replicated == True:
                    
                        # Write line to output file
                        print("Saving ", line_out)
                        file_out.write("%s\t" % filename.split("_")[0])
                        k = 0
                        for item in line_out:
                            k += 1
                            if k < len(line_out):
                                file_out.write("%s\t" % item)
                            elif k == len(line_out):
                                file_out.write("%s\n" % item)
                            elif k > len(line_out):
                                pass
                        
                        # Then start tracking a new mutant and continue to the next line
                        multicelled=False
                        replicated=False
                        line_past=temp # Update line_past
                        line_out=temp[0:4]
                        #print(line_past)
                
                    else:
                        # Then start tracking a new mutant and continue to the next line
                        multicelled=False
                        replicated=False
                        line_past=temp # Update line_past
                        line_out=temp[0:4]
                        #print(line_past)
                
            
            # If same mutant, update org    
            elif temp[0:4]==line_past[0:4] and temp[4] == line_past[4]:
                line_past=temp
                multicelled=True
        
                        
            elif temp[0:4]==line_past[0:4] and temp[4] != line_past[4]:
                line_past=temp
                replicated=True
        
    pf.close() # closing/releasing the pointer in memory.
    file_out.close()

