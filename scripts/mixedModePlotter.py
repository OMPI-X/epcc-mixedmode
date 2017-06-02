#!/usr/bin/python

# SOURCE: http://www2.epcc.ed.ac.uk/~markb/mpiopenmpbench/mixedModePlotter.py
#####################################################################
# Mixed mode MPI/OpenMP benchmark suite v1.0   
#
# Plots mixedmode benchmark data using gnuplot.
#
# For barrier benchmark plots the execution time against
# the number of threads.
# For all other benchmarks plots the execution times are 
# normalised to the execution time with one OpenMP thread 
# per MPI process (i.e. the pure MPI run). 
# 
# The script expects the input files to be of the form produced by
# the mixedModeFileParser script - i.e:
# one file for each run of the benchmark with...
# ..3 columns for the barrier benchmark: (# MPI procs, # Threads, Time/rep)
# ...2 columns for all other benchmarks: (Data size, Time/rep)
#
# Usage: python mixedModeFilePlotter.py 
#
# Note: Program requires the gnuplot.py package to be installed.
#       Find it here: http://gnuplot-py.sourceforge.net/
#####################################################################

import sys
import os
import numpy
try:
    import Gnuplot, Gnuplot.PlotItems, Gnuplot.funcutils
except ImportError:
    print "ERROR importing gnuplot"
    print "gnuplot.py (http://gnuplot-py.sourceforge.net/) needs to be installed."
    raise SystemExit
    #kludge in case Gnuplot hasn't been installed as a module yet:
    #import __init__
    #Gnuplot = __init__
    
#Create a gnuplot object, gplotter
gplotter = Gnuplot.Gnuplot()


def findPlotType():
    #Want to plot for barrier or other benchmark?
    
    print "Enter the plot type: 0 for barrier; any key for other benchmarks..."
    plotType = sys.stdin.readline()
    
    #depending on plotType call barrierPlotter or ratioPlotter
    if (plotType.rstrip('\n') == '0'):
        print "...plotting a barrier benchmark." 
        barrierPlotter()
    else:
        print "...plotting non-barrier benchmark."
        ratioPlotter()


def barrierPlotter():
    # Reads the names of the files containing barrier benchmark data.
    # (File names need to be written in increasing thread count order)
    # Concatenates these files into a file called BarrierPlotter_catFiles.dat
    # Plots the time per reps against the number of threads to an eps file.
    global gplotter
    
    
    while(True):
        #read file name
        print "Enter the name of the file containing barrier data, 0 to exit..."
        fileName = sys.stdin.readline()
        #Remove \n from fileName
        fileName = fileName.rstrip('\n')
        if (fileName != '0'):
            #cat the file to BarrierPlotter_catFiles.dat
            os.system("cat " + fileName + ">> BarrierPlotter_catFiles.dat")
        else:
            #exit out of loop
            break
        
    
    #Find output file name for the plot     
    print("Enter file name for hardcopy of barrier plot...")
    outFileName = sys.stdin.readline()
    outFormat = ".eps"
    outFormat = outFormat.rstrip('\n')
    outFileName = outFileName.rstrip('\n') + outFormat
    
    #Plot barrier data...
    #1) Clear
    gplotter.clear()

    #2) Set axis labels
    gplotter('set border')
    gplotter.xlabel('Number of Threads')
    gplotter.ylabel('Time (seconds)')
    gplotter('set autoscale x')
    gplotter('set autoscale y')
    #gplotter('set format y "10^{%L}"')
  
    #3) Don't want a legend for barrier plot
    gplotter('set key off') #gap between legend entries
    
    #4) Plot BarrierPlotter_catFiles.dat using column 2 (number of threads) and column 3 (time per rep) 
    gplotter.plot(Gnuplot.File("BarrierPlotter_catFiles.dat", using=(2,3),with_='linespoints'))
  
  
    #5) Write hardcopy
    print 'Hardcopy written to:', outFileName
    gplotter.hardcopy(outFileName, eps=1, color=1, enhanced=1)
    
    #6) Uncomment if want to convert the eps output file to a pdf
    #gplotter.close()
    #convertCommand = "epstopdf " + outFileName 
    #os.system(convertCommand)

    #Remove BarrierPlotter_catFiles.dat
    os.system("rm BarrierPlotter_catFiles.dat")

def ratioPlotter():
    # For a particular benchmark, plots the ratio of time per rep 
    # of to a base file against the data size.
    # The base file is (usually) the file containing data for 
    # the 1 OpenMP thread per MPI process benchmark run.
    global gplotter
    
    #Find the output file name
    print("Enter file name for hardcopy of plot...")
    outFileName = sys.stdin.readline()
    #setup file name
    outFormat = ".eps"
    outFormat = outFormat.rstrip('\n')
    outFileName = outFileName.rstrip('\n') + outFormat
    
    
    #gplotter = Gnuplot.Gnuplot()  
    #Clear the gunplot object
    gplotter.clear()

    #Setup axis labels in plot
    gplotter('set border')
    gplotter.xlabel('Message Size (bytes)')
    gplotter.ylabel('Ratio')
    #Use a logscale along x axis..
    gplotter('set logscale x')
    # ..and format of numbers along x axis is 10^
    gplotter('set format x "10^{%L}"')
    
    #Uncomment next two lines for log scale along y axis
    #gplotter('set logscale y')
    #gplotter('set format y "10^{%L}"')
    
    #Setup legend
    gplotter('set key spacing 1.5') #gap between legend entries
    gplotter('set key height 2') #gap b/w top/bottom legend & frame
    gplotter('set key width 2') #gap b/w left/right legend & frame
    gplotter('set key box') #turn on frame around key
    #Line below sets position of the key in the plot - change this for different position
    gplotter('set key top left')
    
    #Start reading in names of data from base file
    #1) Read the name of base file.
    print "Enter name of the file to be used as the divisor (output for 1 OpenMP thread per MPI Process).."
    str = sys.stdin.readline()
    baseDataFileName = str.rstrip('\n') #remove newline from end
    #2) Open the base file
    baseFileHandle = open(baseDataFileName,'r')
    
    #3) Initialise the base file message size and time lists...
    baseFileMsgSize = []
    baseFileTime = []
            
    #4) Read data from baseFile
    str = baseFileHandle.readline()
    while (str != ''): # while not EOF
        str = str.lstrip() # Remove whitespace from start of line
        if ( not str.startswith('#')): # Check if line is a comment
            file1str = str.split() # Split line at whitespaces
            baseFileMsgSize.append(file1str[0]) # Add to message size list ..
            baseFileTime.append(file1str[1]) # ...and to time per rep list
        str = baseFileHandle.readline() # Read next line
    #Close base file after reading
    baseFileHandle.close()

    #Reading data from other files
    while(True):
        #Read file name
        print "Enter the name of file to get the ratio of, 0 to exit"
        ratioFileName = sys.stdin.readline()
        #Remove \n from fileName
        ratioFileName = ratioFileName.rstrip('\n') 
        if (ratioFileName == '0'): #if user eneteed a '0'
            #exit out of loop
            break
        else:
            #Open the file
            ratioFileHandle = open(ratioFileName,'r')
            
            #Initialise the ratio file message size and time lists...
            ratioListMsgSize = []
            ratioListTime = []
    
            
            #Need counter so that only get ratio for same message sizes...
            #..initialise msgSizePosCounter to 0.
            msgSizePosCounter = 0
            #Start reading the file
            line = ratioFileHandle.readline()
            while (line != ''): #while not EOF
                line = line.lstrip() #remove whitespaces from start of line
                if (not line.startswith('#')): #if not a comment
                    msgSizeReadComplete = False
                    lineSplit = line.split() #split the line at whitespace
                    msgSize = lineSplit[0] #first element of line split is message size
                    #Now need to find if the same message size is in baseFileMsgSize list...
                    while (msgSizeReadComplete == False):
                        # ...if we've read more message sizes than length of baseFileMsgSize list
                        if (msgSizePosCounter > len(baseFileMsgSize)-1):
                            # don't use that message size
                            msgSizeReadComplete = True
                        # ...if message size is equal baseFileMsgSize[msgSizePosCounter]
                        elif (int(msgSize) == int(baseFileMsgSize[msgSizePosCounter])):
                            #keep mesaage size and find ratio
                            ratioListMsgSize.append(int(msgSize))
                            ratioListTime.append(float(lineSplit[1]) / float(baseFileTime[msgSizePosCounter]))
                            msgSizeReadComplete = True # finsihed with this message size: go to next one.
                            msgSizePosCounter = msgSizePosCounter + 1
                        # ..if message size if greater than baseFileMsgSize[msgSizePosCounter]
                        elif (int(msgSize) > int(baseFileMsgSize[msgSizePosCounter])):
                            #move to next element in baseFileMsgSize list
                            msgSizePosCounter = msgSizePosCounter + 1
                            msgSizeReadComplete = False
                        # ...if message size if less than baseFileMsgSize[msgSizePosCounter]
                        elif (int(msgSize) < int(baseFileMsgSize[msgSizePosCounter])):
                            #this message file doesn't exist in baseFileMsgSize list: read to next line
                            msgSizeReadComplete = True
                #Read the next line from the file and loop again
                line = ratioFileHandle.readline() 
            #When finished reading file, close it.
            ratioFileHandle.close()
            
            #Now plot the ratios...
            #First need to convert the lists to numpy arrays
            ratioArrayMsgSize = numpy.array(ratioListMsgSize)   
            ratioArrayTime = numpy.array(ratioListTime)
    
            #Find title for the plot 
            splitName = ratioFileName.split('_') #Split at '_'
            procs = splitName[1].rstrip('MPIProcs') #Find number of MPI Procs (2nd element)
            threads = splitName[2].rstrip('Threads.dat') #and number of threads (3rd element)
            t = threads + ' Threads, ' + procs + ' Processes' #create title
    
            #Plot the data    
            gplotter.replot(Gnuplot.Data(ratioArrayMsgSize,ratioArrayTime,with_='linespoints',title=t))
            print "\"", t,"\"", "plotted"
              
    
    gplotter('set xrange[1:100000000] writeback')
    #Plot line along y=1 to represent 1 OpenMP thread per MPI process
    gplotter.replot(Gnuplot.Func('1', with_='lines lt -1',title=None))        
    #Write hardcopy
    print 'Hardcopy written to:', outFileName
    gplotter.hardcopy(outFileName, eps=1, color=1, enhanced=1)
    gplotter.close()
    #change yellow lines to brown
    command = "sed -i 's/^LT5/LT7/g' " + outFileName 
    os.system(command)
    # Ucomment if you want to creat a pdf of plot
    #command = "epstopdf " + outFileName 
    #os.system(command)            
            
            
            

#Program starts executing here....
#Find if plot is barrier plot or other benchmark plot...
#if barrier want to plot # threads vs. time per rep
#if other want to plot ratio time/rep to time/rep for 1 thread (pure MPI) for all other thread counts 
findPlotType()

