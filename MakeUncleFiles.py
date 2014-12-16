'''
Created on Aug 29, 2014

@author: eswens13
'''
from numpy import zeros
import os, subprocess
from random import random


class MakeUncleFiles:


    def __init__(self, atoms):
        """ CONSTRUCTOR """
        
        self.atoms = atoms
        self.structuresInLengths = zeros(len(self.atoms))
    
        self.structList = []
        self.failedStructs = []
        self.case = len(atoms[0].split(','))
        self.pureEnergies = [0.0 for x in range(self.case)]
        
        self.infile = None
        self.holdoutFile = None
        self.outfile = self.infile
        
        self.inCount = 0.0
        self.holdoutCount = 0.0
        
        self.header = "peratom\nnoweights\nposcar\n"
        self.idString = ""
        
        self.vec1x = 0.0
        self.vec2x = 0.0
        self.vec3x = 0.0
        
        self.vec1y = 0.0
        self.vec2y = 0.0
        self.vec3y = 0.0
        
        self.vec1z = 0.0
        self.vec2z = 0.0
        self.vec3z = 0.0
        
        self.atomPositions = []
        self.xPositions = []
        self.yPositions = []
        self.zPositions = []
        
        self.atomCounts = []
        
        self.energy = 0.0
    
    def reinitialize(self):
        """ Re-initializes the class for a new metal atom. """
        self.infile = None
        self.holdoutFile = None
        self.outfile = self.infile
        
        self.inCount = 0.0
        self.holdoutCount = 0.0
        
        self.header = "peratom\nnoweights\nposcar\n"
        self.idString = ""
        
        self.vec1x = 0.0
        self.vec2x = 0.0
        self.vec3x = 0.0
        
        self.vec1y = 0.0
        self.vec2y = 0.0
        self.vec3y = 0.0
        
        self.vec1z = 0.0
        self.vec2z = 0.0
        self.vec3z = 0.0
        
        self.atomPositions = []
        self.xPositions = []
        self.yPositions = []
        self.zPositions = []
        
        self.atomCounts = []
        
        self.energy = 0.0      
    
    def FinishCheck(self, folder):
        """ Tests whether VASP is done by finding "Voluntary" in last line of OUTCAR.  The input
            parameter, folder, is the directory containing OUTCAR, not the OUTCAR file itself. """   
        lastfolder = os.getcwd()
        os.chdir(folder)
        
        proc = subprocess.Popen(['grep', 'Voluntary', 'OUTCAR'],stdout=subprocess.PIPE)
        newstring = proc.communicate()
        os.chdir(lastfolder)   
         
        return newstring[0].find('Voluntary') > -1 #True/False

    def convergeCheck(self, folder, NSW):
        """ Tests whether force convergence is done by whether the last line of OSZICAR (the last
            ionic relaxation step) is less than NSW."""
        try:
            value = self.getSteps(folder)
            return value < NSW #True/False
        except:
            return False #True/False

    def getSteps(self, folder):
        """ Returns the number of ionic relaxation steps that VASP performed as an integer. """
        lastfolder = os.getcwd()
        os.chdir(folder)
        if not os.path.exists('OSZICAR') or os.path.getsize('OSZICAR') == 0:
            os.chdir(lastfolder) 
            return -9999
        oszicar = open('OSZICAR','r')
        laststep = oszicar.readlines()[-1].split()[0]
        oszicar.close()
        os.chdir(lastfolder)  
        try:
            value = int(laststep)
            return value
        except:
            return 9999        

    def setStructureList(self):
        """ Initializes the list of structures to add to the structures.in and structures.holdout
            files by adding only the structures that VASP finished relaxing to the member
            structList. Sorts the list by concentration (N_M / N_total). Adds the structures that
            failed VASP calculations to the member 'failedList'. """   
        self.structList = []
        
        lastDir = os.getcwd()
        
        for atom in self.atoms:
            atomDir = lastDir + '/' + atom
            os.chdir(atomDir)
        
            conclist = []
            atomStructs = []
            failed = []
                        
            dirList = os.listdir(atomDir)
            for item in dirList:
                fullPath = os.path.abspath(item)
                if os.path.isdir(fullPath):
                    if os.path.isdir(fullPath + '/DOS'):
                        if self.FinishCheck(fullPath + '/DOS') and self.convergeCheck(fullPath + '/DOS', 2):
                        
                            # Check for concentration
                            self.setAtomCounts(fullPath)
                        
                            concentration = float(float(self.atomCounts[0]) / self.numberOfSpots)
                        
                            conclist.append([concentration, fullPath])
                        else:
                            failed.append(fullPath)
                    else:
                        # Don't add the 'gss' or 'fits' directories.
                        if not fullPath.split('/')[-1] == 'gss' and not fullPath.split('/')[-1] == 'fits':
                            failed.append(fullPath)
        
            self.failedStructs.append(failed)
            conclist.sort()
            
            for i in xrange(len(conclist)):
                atomStructs.append(conclist[i][1])
            self.structList.append(atomStructs)
            
            os.chdir(lastDir)
    
    def sortStructsByFormEnergy(self, atomInd):
        """ Sorts the list of structures by formation energy. """
        lastDir = os.getcwd()
        os.chdir(lastDir + '/' + self.atoms[atomInd])

        self.pureEnergies = []
        struct = 1
        for nextPureCase in range(self.case,0,-1):
            structDir = lastDir + '/' + self.atoms[atomInd] + '/' + str(struct)

            poscar = open(structDir + '/POSCAR', 'r')
            poscarLines = poscar.readlines()
            poscar.close()

            if poscarLines[0].find('PURE CASE') != -1:
                self.setAtomCounts(structDir)
                self.setEnergy(structDir)        
                self.pureEnergies.append(float(self.energy))

            struct = struct + nextPureCase

        formEnergyList = []
        sortedStructs = []
        for structDir in self.structList[atomInd]:
            self.setAtomCounts(structDir)
            self.setEnergy(structDir)
            structEnergy = float(self.energy)
        
            concentration = float(float(self.atomCounts[0]) / self.numberOfSpots)
                        
            formationEnergy = float(self.energy)
            for i, count in enumerate(self.atomCounts):
                formationEnergy -= float(count) / self.numberOfSpots * self.pureEnergies[i]

            formEnergyList.append([formationEnergy, structDir])
        
        formEnergyList.sort()
        
        for pair in formEnergyList:
            sortedStructs.append(pair[1])
        
        self.structList[atomInd] = sortedStructs
            
        os.chdir(lastDir)
    
    def getStructureList(self):
        """ Returns the list of usable structures along with the list of structures that failed
            VASP calculations. """
        returnList = []
        for i in xrange(len(self.structList)):
            subList = []
            for j in xrange(len(self.structList[i])):
                subList.append(self.structList[i][j].split('/')[-1])
            returnList.append(subList)
        
        failedList = []
        for i in xrange(len(self.failedStructs)):
            subList = []
            for j in xrange(len(self.failedStructs[i])):
                subList.append(self.failedStructs[i][j].split('/')[-1])
            failedList.append(subList)
        
        return returnList, failedList
    
    def getStructuresInLengths(self):
        """ Returns a list of the lengths of the structures.in file for each atom. """
        lengths = zeros(len(self.structuresInLengths))
        for i in xrange(len(self.structuresInLengths)):
            lengths[i] = self.structuresInLengths[i]
        
        return lengths
    
    def setLatticeVectors(self, structFile):
        """ Gets the lattice vectors from the first structure in the structList and sets the 
            corresponding member components. """
        vecFile = open(structFile + '/POSCAR','r')
        vecFileLines = vecFile.readlines()
        vecFile.close()
        
        vec1 = vecFileLines[2].strip().split()
        vec2 = vecFileLines[3].strip().split()
        vec3 = vecFileLines[4].strip().split()
        
        vec1comps = [float(comp) for comp in vec1]
        vec2comps = [float(comp) for comp in vec2]
        vec3comps = [float(comp) for comp in vec3]
            
        
        self.vec1x = vec1comps[0]
        self.vec1y = vec1comps[1]
        if vec1comps[2] == 15.0:
            self.vec1z = 1000.0
        else:
            self.vec1z = vec1comps[2]
        
        self.vec2x = vec2comps[0]
        self.vec2y = vec2comps[1]
        if vec2comps[2] == 15.0:
            self.vec2z = 1000.0
        else:
            self.vec2z = vec2comps[2]
        
        self.vec3x = vec3comps[0]
        self.vec3y = vec3comps[1]
        if vec3comps[2] == 15.0:
            self.vec3z = 1000.0
        else:
            self.vec3z = vec3comps[2]
    
    def closeOutFiles(self):
        """ Closes both the structures.in and structures.holdout files. """
        self.infile.close()
        self.holdoutFile.close()

    def setIDString(self, poscarDir):
        """ Sets the first written line of each structure to the form:
                C H (Metal Atom)  Structure:  #(Decimal) (#(Binary)) """
        poscar = open(poscarDir + '/POSCAR', 'r')
        ID = poscar.readlines()[0].strip()
        
        self.idString = ID

    def setAtomPositions(self, poscarDir):
        """ Retrieves the positions of each of the atoms.  Appends the x-coordinate to the 
            xPositions list, the y-coordinate to the yPositions list.  For a surface in UNCLE the 
            z-coordinate is always zero. """
        poscar = open(poscarDir + '/struct', 'r')
        poscarLines = poscar.readlines()
        poscar.close()
        
        self.atomPositions = []
        self.xPositions = []
        self.yPositions = []
        self.zPositions = []
        
        self.atomPositions = poscarLines[int(7 + self.numberOfC):int(7 + self.numberOfC + self.numberOfSpots)]
        self.atomPositions = [line.strip().split() for line in self.atomPositions]
           
        element = 0
        count = 1
        for pos in self.atomPositions:
            self.xPositions.append(float(pos[0]))
            self.yPositions.append(float(pos[1]))
            self.zPositions.append(0.0)

    def setAtomCounts(self, poscarDir):
        """ Retrieves the number of H atoms and the number of M atoms from the POSCAR file and sets 
            the corresponding members. """
        self.atomCounts = []
        self.vacancyNum = -1
        self.vacancies = 0

        poscar = open(poscarDir + '/POSCAR', 'r')
        poscarLines = poscar.readlines()
        poscar.close()
        
        counts = poscarLines[5].strip().split()
        present = poscarLines[0][poscarLines[0].find(':') + 1:poscarLines[0].find('-')].strip()
        countnum = 1

        for i in range(self.case):
            if present.find(str(i+1)) != -1:
                self.atomCounts.append(int(counts[countnum]))
                countnum = countnum + 1
            else:
                self.atomCounts.append(0)

        if poscarLines[0].find('Vacancies') != -1:
            elements = poscarDir.strip().split('/')[-2].split(',')
            for i, element in enumerate(elements):
                if element.find('Vc') != -1:
                    self.vacancyNum = i
            self.vacancies = int(poscarLines[0].strip().split()[-1])
            self.atomCounts[self.vacancyNum] = self.vacancies

        self.numberOfC = int(counts[0])
        self.numberOfSpots = sum(self.atomCounts)

    def setEnergy(self, directory):  
        """ Retrieves the energy of the structure from the OSZICAR file and sets the corresponding 
            member. """   
        try:
            oszicar = open(directory + '/DOS/OSZICAR','r')
            energy = oszicar.readlines()[-1].split()[2]
            oszicar.close()
        except:
            energy = 0
        
        energy = float(energy)
        peratom = energy / self.numberOfSpots
        
        self.energy = str(peratom)

    def writeHeader(self):
        """ Writes the headers of the structures.in and structures.holdout files. """
        self.infile.write(self.header)
        self.holdoutFile.write(self.header)
    
    def writeDashedLine(self):
        """ Writes a dashed line in the structures.in/.holdout files as a separator between
            different structures. """
        self.outfile.write("#------------------------------------------------\n")
    
    def writeIDString(self):
        """ Writes the ID string of the current structure to either the structures.in or 
            structures.holdout file. """
        concentration = float(float(self.atomCounts[0]) / self.numberOfSpots)
        formationEnergy = float(self.energy)
        for i, count in enumerate(self.atomCounts):
            formationEnergy -= float(count) / self.numberOfSpots * self.pureEnergies[i]

        self.outfile.write(self.idString + " FE = " + str(formationEnergy) + ", Concentration = " + str(concentration) + "\n")
        
    def writeLatticeVecs(self):
        """ Writes the lattice vectors of the current structure to the structures.in or
            structures.holdout file. """
        self.outfile.write("1.0\n")
        self.outfile.write("  %12.8f  %12.8f  %12.8f\n" % (self.vec1z, self.vec1x, self.vec1y))
        self.outfile.write("  %12.8f  %12.8f  %12.8f\n" % (self.vec2z, self.vec2x, self.vec2y))
        self.outfile.write("  %12.8f  %12.8f  %12.8f\n" % (self.vec3z, self.vec3x, self.vec3y))
    
    def writeAtomCounts(self):
        """ Writes the number of H atoms and the number of M atoms to the structures.in or 
            structures.holdout file. """
        for i, count in enumerate(self.atomCounts):
            self.outfile.write(str(count) + " ")
        self.outfile.write("\n")
    
    def writeAtomPositions(self):
        """ Writes the positions of the atoms in the current structure to the structures.in
            or structures.holdout file.  The positions are written in the form:
                z-coord   x-coord   y-coord
            because this is the convention that UNCLE uses. """
        self.outfile.write("Cartesian\n")
        for i in xrange(len(self.atomPositions)):
            self.outfile.write("  %12.8f  %12.8f  %12.8f\n" % 
                               (self.zPositions[i], self.xPositions[i], self.yPositions[i]))
     
    def writeEnergy(self):
        """ Writes the energy of the current structure to the structures.in or structures.holdout
            file. """
        self.outfile.write("#Energy:\n")
        self.outfile.write(str(self.energy) + "\n") 
    
    def writePOSCAR(self, poscarDir, atomInd):
        """ Calls all the methods needed to write all the needed information about the current
            structure to the structures.in or structures.holdout files.  Puts a maximum of 
            approximately 10% of the structures in the structures.holdout file. """

        if self.holdoutCount / float(len(self.structList[atomInd])) < .10 and random() < 0.15: #random() < .15 to get a holdout
            self.outfile = self.holdoutFile
            self.holdoutCount += 1
        elif self.inCount == float(len(self.structList[atomInd])) - 1 and self.holdoutCount == 0:
            self.outfile = self.holdoutFile
            self.holdoutCount += 1
        else:
            self.outfile = self.infile
            self.inCount += 1
        
        #self.outfile = self.infile
        self.setIDString(poscarDir)
        self.setLatticeVectors(poscarDir)
        self.setAtomCounts(poscarDir)
        self.setAtomPositions(poscarDir)
        self.setEnergy(poscarDir)
        
        # Make sure the pure structures go in structures.in
        if self.idString.find('PURE') != -1 and self.outfile == self.holdoutFile:
            self.outfile = self.infile
            self.inCount += 1            
            self.holdoutCount += -1
        
        self.writeDashedLine()
        self.writeIDString()
        self.writeLatticeVecs()
        self.writeAtomCounts()
        self.writeAtomPositions()
        self.writeEnergy()
        
        if self.outfile.name.split('.')[-1] =='holdout':
            return 'holdout'
        else:
            return 'in'

    def makeUncleFiles(self):
        """ Runs through the whole process of creating structures.in and structures.holdout files
            for each metal atom. """
        self.setStructureList()
        
        for i in xrange(len(self.atoms)):
            atomDir = os.getcwd() + '/' + self.atoms[i]
            if os.path.isdir(atomDir):
                subprocess.call(['echo', '\nCreating structures.in and structures.holdout files for ' + self.atoms[i] + '\n'])
                self.reinitialize()
                self.infile = open(atomDir + '/structures.in','w')
                self.holdoutFile = open(atomDir + '/structures.holdout','w')
                self.sortStructsByFormEnergy(i)
                self.writeHeader()
                
                num = 0
                structuresInCount = 0
                for structure in self.structList[i]:
                    if structuresInCount >= 500:    # Write a maximum of 500 structures to the file
                        break                       # for any given atom.
                    whichFile = self.writePOSCAR(structure, i)
                    if whichFile == 'in':
                        structuresInCount += 1
                    num += 1
                
                self.structuresInLengths[i] = structuresInCount
                self.closeOutFiles()









