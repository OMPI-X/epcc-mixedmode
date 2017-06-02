#!/usr/bin/python

# SOURCE: http://www2.epcc.ed.ac.uk/~markb/mpiopenmpbench/mixedModeFileParser.py
#####################################################################
# Mixed mode MPI/OpenMP benchmark suite v1.0   
#
# Parser for output files.
#
# Writes that data size and time per rep columns in the output file 
# from a run of the benchmark suite to individual 
# files for each benchmark. 
# The name of the files take the form:
# "<benchmarkName>_<#MPIProcs>MPIProcs_<#Threads>Threads.dat"    
#
# Usage: python mixedModeFileParser.py <inputFile>
#####################################################################

import sys
import string

inFile = []
benchParamStrings = ["Number of MPI processes", "Number of OpenMP threads"]
benchParams = []
benchmarkName = []
ppType = []
dataSize = []
timePerRep = []

def openInFile():
    #Reads the file name from the command line and opens.
    global inFile #file handle for input file

    #Read filename from argument list
    if (len(sys.argv) != 2):
        print "ERROR: No filename supplied"
        print "Usage: python mixedModeFileParser.py <inputFileName>"
        raise SystemExit
        
    print "Attempting to open", sys.argv[1], "...",
    #Try to Open filename
    try:
        inFile = open(sys.argv[1],'r')
    except IOError:
        print "ERROR opening file."
        raise SystemExit
    else:
        print "File opened ok."
            
    return        
    
def getBenchmarkParams():
    #Reads the number of MPI processes & OpenMP threads from
    #the input file.
    global inFile, benchParamStrings
    
    for param in benchParamStrings:
        #Set paramFound to False
        paramFound = False
        
        while (paramFound != True):
            #Read a line from the input file
            line = inFile.readline()
            if (param in line):
                #Tokenise 'line', splitting at '='
                lineSplit = line.split('=')
                #Add the number to benchParams 
                #(and remove new line and white space).
                benchParams.append(lineSplit[1].strip('\n '))
                #set paramFound to True
                paramFound = True
    
    return
                
def parseinFile():
    #Parses the input file.
    #Implemented as a state machine with 5 states:
    #1. finding the name of the benchmark
    global inFile, ppType, dataSize, timePerRep
    
    #initialise state variables
    findingBenchName = True
    checkingForPPBench = False
    readingData = False
    
    #Start loop over each line of inFile
    for lines in inFile:
        
        # remove the newline from the end of 'lines'
        lines.strip('\n')
        #State 1 - benchmark name
        if (findingBenchName == True):
            findingBenchName = getBenchName(lines)
            if (findingBenchName == False):
                checkingForPPBench = True # move onto next state
        #State 2 - finding if benchmark is inter/intra node pingping/pingpong
        elif (checkingForPPBench == True):
            checkingForPPBench = getPPType(lines)
            if (checkingForPPBench == False):
                #print "PP type = " + ppType
                readingData = True
        #State 3 - read dataSize & timePerRep for benchmark
        elif (readingData == True):
            readingData = getBenchData(lines)
            if (readingData == False):
                #State 4 - write this data to file
                writeToFile()
                # and start from state 1 again.
                findingBenchName = True
                
    writeToFile()
        
    #Finished with input file - close it.
    inFile.close()
    
    return

            
def getBenchName(line):
    # Checks if the line contains the '#' character.
    #If yes, copies contents of line to benchmarkName
    global benchmarkName
    
    findingName = True
    
    #Remove leading and trailing whitespace
    line = line.strip(' \n')
    if (line.startswith('#')):
        #remove the '#'
        line = line.strip('#')
        #remove any whitespace
        line = line.replace(" ", "")
        benchmarkName = line
        findingName = False
    
    return findingName
        
def getPPType(line):
    #Checks if the benchmark is a pingping or a pingpong.
    #If yes, then it finds if it was inter or intra node.
    global benchmarkName, ppType
    
    findingPPType = True
    
    #For multi-PP benchmark not looking for intra or inter node info.
    if (benchmarkName.find("MultiPing") != -1):
        ppType = ''
        findingPPType = False
    #Check if benchmark is pingpong or pingping
    elif (benchmarkName.find("Ping") != -1):
        if (line.find("Inter") != -1): 
            ppType = "Inter"
            findingPPType = False
        elif (line.find("Intra") != -1):
            ppType = "Intra"
            findingPPType = False
        #need to search next line if not found
        else:
            findingPPType = True
    #benchmark isn't multi-pingping/pingpong or pingping/pingpong
    else:
        ppType = ''
        findingPPType = False
    
    return findingPPType

def getBenchData(line):
    #Checks if 'line' contains benchmark data (begins with 'd')
    #Splits these lines at whitespaces and fills the dataSize & timePerRep variables.
    
    global dataSize, timePerRep
    
    gettingBenchData = True
    
    #Remove leading and trailing whitespace
    line = line.strip(' ')
    #Check if line begins with 'd'
    if (line.startswith('d')):
        #Split the line at whitespace
        lineSplit = line.split()
        dataSize.append(lineSplit[2]) #Add the 3rd element of lineSplit to dataSize list
        timePerRep.append(lineSplit[5]) #Add 6th element of lineSplit to timePerRep list
    #If line doesn't begin with a 'd' and there are elements in dataSize list finish getBenchData
    elif (len(dataSize) > 0):
        gettingBenchData = False
    
    return gettingBenchData

def writeToFile():
    # Writes message size and time per rep data to 
    #individual files for each benchmark.
    global benchmarkName, ppType, dataSize, timePerRep, benchParams
    
    #Construct name of output file for benchmark.
    writeFileName = benchmarkName + ppType + "_" + benchParams[0] + "MPIProcs_" + benchParams[1] + "Threads.dat"
    
    #open the output file
    try:
        outFile = open(writeFileName,'w')
    except IOError:
        print "ERROR opening output file: " + writeFileName
    
    #Print a header to outFile
    headerString = "# " + benchmarkName + " benchmark.\n"
    headerString += "# Using " + benchParams[0] + " MPI processes and " \
    + benchParams[1] + " OpenMP threads \n"

    outFile.write(headerString)
    
    
    #Print benchmark data
    #1) For barrier benchmark
    if (benchmarkName == "Barrier"):
        #Print header for barrier
        outFile.write("# Procs \t Threads \t Time/Rep\n")
        #print one line of data
        outFile.write(benchParams[0] + "\t" + benchParams[1] + "\t" + timePerRep[0] + "\n")
    #2) For all other benchmarks
    else:
        #Print header
        outFile.write("# Message Size \t Time/Rep\n")
        #Loop over data and print
        for i in range(len(dataSize)):
            outFile.write(dataSize[i] + "\t" + timePerRep[i] + "\n")
    
    #Print write confirm message.
    print "Data written to " + writeFileName + "\n";
    #Close outFile
    outFile.close()
    
    #Clear the dataSize and timePerRep lists.
    dataSize = []
    timePerRep = []
                
    return
          
#Program starts executing here....
#1) Open input file
openInFile()
#2) Get number of MPI processes and OpenMP threads
getBenchmarkParams()
#3) Start parsing input file
parseinFile()
