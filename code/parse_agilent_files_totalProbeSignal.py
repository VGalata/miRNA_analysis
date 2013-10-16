#!/usr/bin/python
# -*- coding: utf-8 -*-


import os, re, sys, collections
import numpy as np



def parseAgilentFileMiRNA(filename, resultlist, uniq_mirnas, spezies="hsa"):
  input = open(filename, "r")
  headerFound = False
  header = []
  nametoindex = {}
  tmp_dict = {}
  mirna_probes_values=collections.defaultdict(lambda: collections.defaultdict(list))

  for line in input:
    if not headerFound:
      if line.find("SystematicName") != -1:
        headerFound = True
        header = line.strip().split("\t")
        nametoindex = dict((header[i],i ) for i in range(len(header)) )
        #print nametoindex
    else: # parse data
      entries = line.strip().split("\t")
      mirna = entries[nametoindex['SystematicName']]
      #print mirna
      control_type = entries[nametoindex['ControlType']]
      probe_name = entries[nametoindex['ProbeName']]
      probe_signal = float(entries[nametoindex['gTotalProbeSignal']])
      probe_error = float(entries[nametoindex['gTotalProbeError']])
      is_detected = int(entries[nametoindex['gIsGeneDetected']])
      #handle detected genes as agilent would:
      if is_detected:
        #check if probe_signal > 3x probe_error
        if probe_signal < (3 * probe_error):
            continue
         
      #do not parse blanks and controls
      if control_type == "0": #1: positive control, -1: negative control
        #print "no control"
        if spezies: #if spezies string is empty we parse all miRNAs
            #print spezies
            if mirna.find(spezies + "-") != -1:#only parse spezies specific mirnas
                uniq_mirnas.add(mirna)
                mirna_probes_values[mirna][probe_name].append(probe_signal)
                #print mirna
        else:#we also parse viral mirnas
            #print "HIER"
            uniq_mirnas.add(mirna)
            mirna_probes_values[mirna][probe_name].append(probe_signal)
      
  input.close()
  #print mirna_probes_values
  resultlist.append(mirna_probes_values)
  return resultlist, uniq_mirnas


def parseAgilentFileMRNA(filename, resultlist, uniq_mrnas):
  input = open(filename, "r")
  headerFound = False
  header = []
  nametoindex = {}
  tmp_dict = {}
  mrna_probes_values=collections.defaultdict(lambda: collections.defaultdict(list))

  for line in input:
    if not headerFound:
      if line.find("SystematicName") != -1:
        headerFound = True
        header = line.strip().split("\t")
        nametoindex = dict((header[i],i ) for i in range(len(header)) )
    else: # parse data
      entries = line.strip().split("\t")
      mrna = entries[nametoindex['SystematicName']]
      control_type = entries[nametoindex['ControlType']]
      probe_name = entries[nametoindex['ProbeName']]
      median_signal = float(entries[nametoindex['gMedianSignal']])

         
      #do not parse blanks and controls
      if control_type == "0": #1: positive control, -1: negative control
        uniq_mrnas.add(mrna)
        mrna_probes_values[mrna][probe_name].append(median_signal)
      
  input.close()
  resultlist.append(mrna_probes_values)
  return resultlist, uniq_mrnas



def parseAgilentFiles(list_of_files, spezies="hsa", array_type = "miRNA"):
  resultlist = []
  uniq_mirnas = set()

  for filename in list_of_files:
    print "parsing: " + filename
    if array_type == "miRNA":
        #print "miRNA"
        #print spezies
        resultlist, uniq_mirnas = parseAgilentFileMiRNA(filename, resultlist, uniq_mirnas, spezies)
    else: #mRNA
        #print "mRNA"
        resultlist, uniq_mirnas = parseAgilentFileMRNA(filename, resultlist, uniq_mirnas)
    
  return resultlist, uniq_mirnas


#parse logfile
def parse_logfile(filename, array_name = "miRNA"):
    input = open(log_file, "r")
    samplenames = []
    filenames = []

    nametoix = {}
    #chip2info = collections.defaultdict(dict)
    #chip2arrays = collections.defaultdict(list)
    for line in input:
        if line.strip() != "" and line.find('#', 0, 1) == -1: #do not parse empty or comment lines
            if line.find("Chip") != -1:
                header = line.strip().split("\t")
                nametoix = dict((header[i],i ) for i in range(len(header)))
            else:  
                entries = line.strip().split("\t")
                chip = entries[nametoix["Chip"]]
                array = entries[nametoix["Array"]]
                sample = entries[nametoix["Probe"]]
                ## BEGIN: 2013.09.04 Valentina
                if len(entries)==3: # execute Christina's code
                    filename = "US11153896_" + chip + "_S*_"+array_name+"_1010_Sep10_"+array+".txt"#TODO abfangen mehr als eine datei pro sample!!!
                else: # execute modified code
                    us_number = entries[nametoix["US_Num"]]
                    s_id = entries[nametoix["S_ID"]]
                    number = entries[nametoix["Number"]]
                    date = entries[nametoix["Date"]]
                    array_name = entries[nametoix["Array_Name"]]
                    filename = us_number + "_" + chip + "_" + s_id + "_" + array_name + "_" + number + "_" + date + "_" + array + ".txt"#TODO abfangen mehr als eine datei pro sample!!!
                samplenames.append(sample)
                filenames.append(filename)
                ## END: 2013.09.04 Valentina
                ## WAS
                ##filename = chip + "_S01_" + array + "_GeneView.txt"
                #filename = "US11153896_" + chip + "_S*_"+array_name+"_1010_Sep10_"+array+".txt"#TODO abfangen mehr als eine datei pro sample!!!
                ##chip2info[chip + "_" + array]["filename"] = filename
                ##chip2info[chip + "_" + array]["samplename"] = sample
                ##chip2arrays[chip].append(array)
                #samplenames.append(sample)
                #filenames.append(filename)                
    input.close()

    #if not we have not unique names
    assert len(samplenames) == len(set(samplenames))
    assert len(filenames) == len(set(filenames))
    return filenames, samplenames


if __name__ == '__main__':
    if len(sys.argv) != 6:
        sys.exit("usage: " + sys.argv[0] + " <logfile_scanner> <basefolder_agilent_files> <organism [hsa]> <type [miRNA | GE1]> <output_expr_matrix>")
    
    log_file = sys.argv[1]
    basefolder = sys.argv[2]
    organism = sys.argv[3]
    array_name = sys.argv[4]
    outfile = sys.argv[5]
    #print "organism: *"+ organism+"*" 
    mRNA_arrays = ['GE1']#in case we have some other names for mRNA arrays
    if array_name in mRNA_arrays:
        array_type = "mRNA"
    else:
        array_type = "miRNA"
    filenames, samplenames = parse_logfile(log_file, array_name)
    print filenames
    print samplenames
    
    #find GeneView files
    full_filenames = []
    for filename in filenames:
        print "searching ... " + filename
        full_path_name = os.popen("find \"" + basefolder + "\" -name " + filename ).read().strip()
        #print full_path_name
        if not full_path_name:
            print "\tNOT FOUND!"
        else:
            #check if we have more than one file with this name
            full_path_names = full_path_name.split("\n")
            if len(full_path_names) > 1:
                print "warning: more than one possible file for chip found!"
                print full_path_names
            full_filenames.append(full_path_names[0]) #we take the first hit either way

    assert len(full_filenames) != 0, 'feature extraction files not found!'
    assert len(full_filenames) == len(filenames), 'Not all files found!'

    resultlist, uniq_mirnas = parseAgilentFiles(full_filenames, organism, array_type)


    output = open(outfile, "w")
    #print header
    for sample in xrange(len(samplenames)):
        output.write("\t" + samplenames[sample])
    output.write("\n")
    #output matrix
    if array_type == "miRNA":
        for mirna in uniq_mirnas:
            output.write(mirna)
            for sample in xrange(len(samplenames)):
                #compute a self-made totalGeneSignal from the sum of the means of the different totalProbeSignals per miRNA
                tgs = 0.0
                for probe, values in resultlist[sample][mirna].items():
                    probe_mean = np.mean(values)
                    tgs += probe_mean
                output.write("\t" + str(tgs))
            output.write("\n")
    else:#mRNA array -> don't summarize probes
        for mirna in uniq_mirnas:
            output.write(mirna)
            for sample in xrange(len(samplenames)):
            #probes = resultlist[0][mirna]#NOTE hope that we have the same number of probes for all mRNAs
            #print probes
                median_probes_per_mirna = []
                for probe, values in resultlist[sample][mirna].items():
                    #print probe, values
                    median_probes_per_mirna.append(np.median(values))
                output.write("\t" + str(np.median(median_probes_per_mirna)))
            output.write("\n")
                #output.write(mirna+ "\t" + probe)
                #for number in xrange(len(samplenames)):
                    #output.write("\t" + str(np.median(resultlist[number][mirna][probe])))
                #output.write("\n")

    output.close()
"""
        #output each probe, but summarize value per probe
        for mirna in uniq_mirnas:
            #output.write(mirna)
            
            #for sample in xrange(len(samplenames)):
            probes = resultlist[0][mirna]#NOTE hope that we have the same number of probes for all mRNAs
            #print probes
            for probe in probes:
                #print probe, values
                output.write(mirna+ "\t" + probe)
                for number in xrange(len(samplenames)):
                    output.write("\t" + str(np.median(resultlist[number][mirna][probe])))
                output.write("\n")
"""

