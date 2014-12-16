'''
Created on Aug 26, 2014

@author: eswens13
'''

import os, subprocess
from random import seed
from numpy import zeros
from copy import deepcopy

import Enumerator, Extractor, Structs2Poscar, JobManager, MakeUncleFiles, Fitter, GSS, Analyzer, DistanceInfo

def readSettingsFile():
    currDir = os.getcwd()
    infile = open(currDir + '/needed_files/settings.in', 'r')
    inlines = []
    for line in infile:
        firstPart = line.strip().split()[0]
        firstChar = list(firstPart)[0]
        if firstChar != '#':
            inlines.append(line.strip())
    infile.close()
    
    atoms = []
    volRange = []
    clusterNums = []
    trainStructs = 0
    fitStructs = 0
    fitSubsets = 0
    plotTitle = "title"
    xlabel = "xlabel"
    ylabel = "ylabel"
    case = 3
    
    for line in inlines:
        if line.split()[0] == 'ATOMS:':
            while len(line) > 0:
                low = line.find('[')
                high = line.find(']')
                if low == -1 or high == -1:
                    break
                atoms.append(line[low+1:high].translate(None, " []"))
                line = line[high+1:]
            if len(atoms) == 0:
                atoms = [line[line.find(':')+1:line.find('#')].strip().replace(' ', ',')]

        elif line.split()[0] == 'CASE:':
            case = int(line.split()[1])
            
        elif line.split()[0] == 'VOL_RANGE:':
            low = int(line.split()[1])
            high = int(line.split()[2])
            volRange = [low, high]
            
        elif line.split()[0] == 'CLUSTER_NUMS:':
            parts = line.split()
            for i in xrange(1, 11):
                clusterNums.append(int(parts[i]))
        
        elif line.split()[0] == 'TRAINING_STRUCTS:':
            trainStructs = int(line.split()[1])
          
        elif line.split()[0] == 'FITTING_STRUCTS:':
            fitStructs = int(line.split()[1])
            
        elif line.split()[0] == 'STRUCT_SUBSETS:':
            fitSubsets = int(line.split()[1])
        
        elif line.split()[0] == 'PLOT_TITLE:':
            plotTitle = line.split('\'')[1]
        
        elif line.split()[0] == 'XLAB:':
            xlabel = line.split('\'')[1]
        
        elif line.split()[0] == 'YLAB:':
            ylabel = line.split('\'')[1]


    if len(atoms[0].split(',')) != case:
        print "\nERROR: Number of atoms does not match case. Check settings.in\n"
        print "Setting case to " + str(len(atoms[0].split(','))) + "-nary. . .\n"
        case = len(atoms[0].split(','))
    
    return [atoms, volRange, clusterNums, trainStructs, fitStructs, fitSubsets, plotTitle, xlabel, ylabel, case]

def contains(struct, alist):
    for i in xrange(len(alist)):
        if str(struct) == str(alist[i]):
            return True
    
    return False

def equals(alist, blist):
    clist = deepcopy(alist)
    dlist = deepcopy(blist)
    while len(clist) > 0:
        if len(dlist) == 0:
            return False
        if not contains(clist[0], dlist):
            return False
        else:
            clist.remove(clist[0])
            dlist.remove(clist[0])
    
    if len(dlist) > 0:
        return False
    else:
        return True
          
if __name__ == '__main__':
    seed()
    
    [atomList, volRange, clusterNums, trainingStructs, fitStructs, fitSubsets, plotTitle, xlabel, ylabel,case] = readSettingsFile()
    uncleOutput = open('uncle_output.txt','w') # All output from UNCLE will be written to this file.
    
    enumerator = Enumerator.Enumerator(atomList, volRange, clusterNums, trainingStructs, uncleOutput, case)
    subprocess.call(['echo','\nEnumerating symmetrically unique structures. . .\n'])
    enumerator.enumerate()
    
    changed = True
    iteration = 1
    
    newStructs = []
    gssStructs = []
    lowestStructsFile = open('lowest_vasp.txt','w')
    lowestGssFile = open('lowest_gss.txt','w')
    failedFile = open('failed_vasp.txt','w')
    
    while changed:
        changed = False
        
        subprocess.call(['echo','\n========================================================'])
        subprocess.call(['echo','\t\tIteration ' + str(iteration)])
        subprocess.call(['echo','========================================================\n'])
        
        # Extract the pseudo-POSCARs from struct_enum.out
        extractor = Extractor.Extractor(atomList, uncleOutput,case)
        if iteration == 1:
            extractor.setTrainingStructs()
        elif iteration > 1:
            extractor.setStructList(newStructs)

        toCalculate = extractor.getStructList()
        extractor.extract()
    
        # Convert the extracted pseudo-POSCARs to VASP POSCAR files, make directories for them
        # and put the POSCARs in their corresponding directories.
        subprocess.call(['echo','\nConverting outputs to VASP inputs. . .\n'])
        toPoscar = Structs2Poscar.Structs2Poscar(atomList, toCalculate)
        toPoscar.convertOutputsToPoscar()
     
        # Start VASP jobs and wait until they all complete or time out.
        manager = JobManager.JobManager(atomList)
        manager.runLowJobs(toCalculate)
        manager.runNormalJobs(toCalculate)
        manager.runDOSJobs(toCalculate)
    
        # Create structures.in and structures.holdout files for each atom.
        uncleFileMaker = MakeUncleFiles.MakeUncleFiles(atomList)
        uncleFileMaker.makeUncleFiles()
        
        # Get all the structs that have been through VASP calculations for each atom. These
        # should be sorted by formation energy during the work done by makeUncleFiles()
        [vaspStructs, failedStructs] = uncleFileMaker.getStructureList()
        structuresInLengths = uncleFileMaker.getStructuresInLengths() 
        
        # Perform a fit to the VASP data in structures.in for each atom.
        fitter = Fitter.Fitter(atomList, fitStructs, fitSubsets, structuresInLengths, uncleOutput)
        fitter.makeFitDirectories()
        fitter.fitVASPData(iteration)
    
        # Perform a ground state search on the fit for each atom.    
        gss = GSS.GSS(atomList, volRange, plotTitle, xlabel, ylabel, uncleOutput)
        gss.makeGSSDirectories()
        gss.performGroundStateSearch(iteration)
        gss.makePlots(iteration)
        gssStructs = gss.getAllGSSStructures(iteration, failedStructs)
        
        # Check the lowest 100 hundred structs from VASP against the lowest 100 structs from UNCLE
        # for each atom.  If they match, then that atom has converged and we remove it from the 
        # lists.
        removeAtoms = []
        removeGss = []
        removeVasp = []
        for i in xrange(len(vaspStructs)):
            atomLength = len(vaspStructs[i])
            if atomLength >= 100:
                if equals(vaspStructs[i][:100], gssStructs[i][:100]):
                    removeAtoms.append(atomList[i])
                    removeGss.append(gssStructs[i])
                    removeVasp.append(vaspStructs[i])
            else:
                # If there are not yet 100 structs that have converged in VASP.
                if equals(vaspStructs[i][:atomLength], gssStructs[i][:atomLength]):
                    removeAtoms.append(atomList[i])
                    removeGss.append(gssStructs[i])
                    removeVasp.append(vaspStructs[i])
        
        for i in xrange(len(removeAtoms)):
            atomList.remove(removeAtoms[i])
        for i in xrange(len(removeGss)):
            gssStructs.remove(removeGss[i])
        for i in xrange(len(removeVasp)):
            vaspStructs.remove(removeVasp[i])
        
        # If all of the atoms have converged, exit the loop.  Else, keep going.
        if len(atomList) > 0:
            changed = True
                    
        if not changed:
            subprocess.call(['echo','\n----------------- The loop has converged! ---------------'])

        # Print the lowest energy structures that have been through VASP calculations to a file.
        try:
            lowestStructsFile.write('==============================================================\n')
            lowestStructsFile.write('\tIteration: ' + str(iteration) + '\n')
            lowestStructsFile.write('==============================================================\n')
            for i in xrange(len(vaspStructs)):
                lowestStructsFile.write('\n******************** ' + atomList[i] + ' ********************\n')
                atomLength = len(vaspStructs[i])
                if atomLength >= 100:
                    for j in xrange(len(vaspStructs[i][:100])):
                        if (j + 1) % 20 == 0 or j == 99:
                            lowestStructsFile.write(str(vaspStructs[i][j]) + '\n')
                        else:
                            lowestStructsFile.write(str(vaspStructs[i][j]) + ', ')
                else:
                    for j in xrange(len(vaspStructs[i][:atomLength])):
                        if (j + 1) % 20 == 0 or j == atomLength - 1:
                            lowestStructsFile.write(str(vaspStructs[i][j]) + '\n')
                        else:
                            lowestStructsFile.write(str(vaspStructs[i][j]) + ', ')
            lowestStructsFile.flush()
            os.fsync(lowestStructsFile.fileno())
        except IOError:
            subprocess.call(['echo','\n~~~~~~~~~~ Couldn\'t write to lowest_vasp file. ~~~~~~~~~~\n'])
        
        # Write the all the structures that have failed VASP calculations to a file.
        # TODO: Only write the failed structures that are unique to this iteration to the file.
        try:
            failedFile.write('==============================================================\n')
            failedFile.write('\tIteration: ' + str(iteration) + '\n')
            failedFile.write('==============================================================\n')
            for i in xrange(len(failedStructs)):
                failedFile.write('\n******************** ' + atomList[i] + ' ********************\n')
                for j in xrange(len(failedStructs[i])):
                    if (j + 1) % 20 == 0 or j == len(failedStructs[i]) - 1:
                        failedFile.write(str(failedStructs[i][j]) + '\n')
                    else:
                        failedFile.write(str(failedStructs[i][j]) + ', ')
            failedFile.flush()
            os.fsync(failedFile.fileno())
        except IOError:
            subprocess.call(['echo','\n~~~~~~~~~~ Couldn\'t write to failed_vasp file. ~~~~~~~~~~\n'])
                
        # Add the 100 structures with the lowest formation energy that have not been through VASP
        # calculations (converged or failed) to the newStructs list for each remaining atom.
        newStructs = []
        added = zeros(len(atomList))
        for i in xrange(len(gssStructs)):
            atomStructs = []
            for j in xrange(len(gssStructs[i])):
                if added[i] >= 100:
                    break
                elif not contains(gssStructs[i][j], vaspStructs[i]):
                    atomStructs.append(str(gssStructs[i][j]))
                    added[i] += 1
            newStructs.append(atomStructs)
        
        # Print the new GSS structures to a file.
        try:
            lowestGssFile.write('==================================================================\n')
            lowestGssFile.write('\tIteration: ' + str(iteration) + '\n')
            lowestGssFile.write('==================================================================\n')
            for i in xrange(len(newStructs)):
                lowestGssFile.write('\n***************** ' + atomList[i] + ' *******************\n')
                atomLength = len(newStructs[i])
                if atomLength >= 100:
                    for j in xrange(len(newStructs[i])):
                        if (j + 1) % 20 == 0 or j == 99:
                            lowestGssFile.write(str(newStructs[i][j]) + '\n')
                        else:
                            lowestGssFile.write(str(newStructs[i][j]) + ', ')
            lowestGssFile.flush()
            os.fsync(lowestGssFile.fileno())
        except IOError:
            subprocess.call(['echo','\n~~~~~~~~~~ Couldn\'t write to lowest_gss file. ~~~~~~~~~~\n'])
         
        # Keep track of which iteration we're on.
        iteration += 1
    
    uncleOutput.close()
    lowestStructsFile.close()
    lowestGssFile.close()
    # Should do some analysis after the loop has finished as well.
        

    
        
    
    
 
 
 
 
 
 
 
 
    
