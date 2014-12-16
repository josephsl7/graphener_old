'''
Created on Aug 13, 2014

@author: eswens13
'''
from glob import glob
from math import sqrt
import os, subprocess, shutil


class Structs2Poscar:
    """ Converts the pseudo-POSCAR files generated by the Extractor class to standard POSCAR files
        that VASP can use.  The Converter class below is the one that does most of the work in 
        conversion.  This class is the administrative class, managing the files and directories.  
        It creates directories for each of the metal atoms specified by the user and, within each 
        atom's directory, creates directories for each of the structures to be run through VASP
        calculations.  It populates these directories with the converted POSCAR file corresponding 
        to that structure. """
        
    def __init__(self, atoms, structList):
        """ CONSTRUCTOR """
        self.atoms = atoms
        self.structList = structList
    
    def makeAtomDirectories(self):
        """ Creates a directory for each atom in the atom list specified in settings.in.  All the 
            VASP and UNCLE files for the atom will be placed in this directory. """
        for atom in self.atoms:
            atomDir = os.getcwd() + '/' + atom
            if not os.path.isdir(atomDir):
                subprocess.call(['mkdir',atomDir])
    
    def makePlots(self, plotDir):
        """ NOT WORKING CURRENTLY.  This method puts a visual representation of each structure in 
            its directory. """
        toPlot = os.listdir(self.getExportDir())
        for struct in toPlot:
            subprocess.call(['python', 'PlotGraphene.py', self.getPlotDir(), struct, '-u'])
    
    def convertOutputsToPoscar(self):
        """ Converts all the pseudo-POSCARs created by the Extractor class into POSCAR files that
            VASP can use to run calculations and puts them in a directory that will contain all
            the information for that structure. """
        self.makeAtomDirectories()
        
        for name in glob('enum/vasp.0*'):
            structNum = str(self.retrieveStructNum(name))
            self.changeToPoscar(name)
            
            for i in xrange(len(self.structList)):
                if self.contains(structNum, self.structList[i]):
                    structDir = os.getcwd() + '/' + self.atoms[i] + '/' + structNum
                    if os.path.isdir(structDir):
                        shutil.rmtree(structDir)
                    subprocess.call(['mkdir', structDir])                    
            
                    infile = open('POSCAR', 'r')
                    inLines = infile.readlines()
                    infile.close()

                    poscar = open(structDir + '/POSCAR','w')

                    line = inLines[0]
                    elements = ['C'] + self.atoms[i].split(',')

                    poscar.write(line[:line.find(':') + 1] + ' ')

                    vacancies = 0

                    for num, element in enumerate(elements):
                        if self.atomCounts[num] != 0 and element.find('Vc') == -1 and num != 0:
                            poscar.write(str(num) + ' ')
                        elif element.find('Vc') != -1:
                            vacancies = self.atomCounts[num]

                    if vacancies == 0:
                        poscar.write(line[line.find('-'):])
                    else:
                        poscar.write(line[line.find('-'):-1].strip() + ' -- ' + 'Vacancies: ' + str(vacancies) + '\n')

                    poscar.write(inLines[1] + inLines[2] + inLines[3] + inLines[4])

                    countLine = ''        

                    for num, count in enumerate(self.atomCounts):
                        if count != 0 and elements[num].find('Vc') == -1:
                            countLine = countLine + str(count) + ' '
                    poscar.write(countLine.strip() + '\n')

                    poscar.write(inLines[6])

                    linenum = 7
                    for num, count in enumerate(self.atomCounts):
                        if count != 0:
                            for line in range(count):
                                if elements[num].find('Vc') == -1:
                                    pos = inLines[linenum].strip().split()
                                    pos = [float(comp) for comp in pos]
                                    z = pos[2]
                                    if elements[num] == 'H':
                                        if z >= 0:
                                            z += 1.1
                                        else:
                                            z += -1.1
                                    elif elements[num] != 'C':
                                        if z >= 0:
                                            z += 2.2
                                        else:
                                            z += -2.2
                                    poscar.write('%12.8f %12.8f %12.8f\n' % (pos[0], pos[1], z))
                                linenum += 1
                           
                    poscar.close()
                    subprocess.call(['cp', 'POSCAR', structDir + '/struct'])
            subprocess.call(['rm', name])
            subprocess.call(['rm','POSCAR'])
            
    def retrieveStructNum(self, structFile):
        """ Returns the structure number from the name of the pseudo-POSCAR file.  For example,
            the file "vasp.000478" would return 478 as the structure number. """
        fileChars = list(structFile)
        i = len(fileChars) - 1
        while fileChars[i] != '.':
            i -= 1

        return int(''.join(fileChars[i + 1:]))

    def contains(self, struct, alist):
        """ Retuns true if the list 'alist' contains the structure 'struct', false otherwise. """
        for i in xrange(len(alist)):
            if struct == alist[i]:
                return True
        
        return False

    def changeToPoscar(self, structFile):
        """ The hook for the Converter class.  Uses the Converter class to convert a pseudo-POSCAR
            file into a standard VASP POSCAR file. """
        infile = open(structFile, 'r')
        inLines = infile.readlines()
        infile.close()
        
        converter = Converter(structFile)
        converter.convert()

        self.atomCounts = converter.getAtomCounts()
        
        poscar = open('POSCAR','w')
        
        poscar.write("Contains atoms: ")
        for i, count in enumerate(self.atomCounts):
            if count != 0 and i != 0:
                poscar.write(str(i) + " ")
        
        if converter.isPure():
            poscar.write("- (PURE CASE) ")

        poscar.write("- " + inLines[0])
            
        poscar.write('1.0\n')
        
        # These are in the 'scrambled', (z, x, y) form from UNCLE.  We need to unscramble them.
        latticevecs = converter.getLatticeVectors()
        
        for vec in latticevecs:
            poscar.write('%12.8f %12.8f %12.8f\n' % (vec[1], vec[2], vec[0]))
        
        countLine = ''        

        for count in self.atomCounts:
            if count != 0:
                countLine = countLine + str(count) + ' '
        poscar.write(countLine.strip() + '\n')
        
        poscar.write('Cartesian\n')
        
        Cpos = converter.getCPositions()
        for pos in Cpos:
            poscar.write('%12.8f %12.8f %12.8f\n' % (pos[0], pos[1], pos[2]))
        
        Pos = converter.getPositions()
        for pos in Pos:
            poscar.write('%12.8f %12.8f %12.8f\n' % (pos[0], pos[1], pos[2]))
        
        poscar.close()
        

class Converter:
    """ Converts one of the files made from the 'makestr.x' routine in UNCLE to a POSCAR file
        that is ready to be used by VASP. """
        
    def __init__(self, uncleFile):
        """ CONSTRUCTOR """
        self.uncleFile = uncleFile
        self.atomCounts = []
        
        self.lattVec1 = []
        self.lattVec2 = []
        self.lattVec3 = []
        
        # 2D distance between hexagonal C atoms
        self.distance = 1.42085901
        
        # Distance from the plane when buckled
        self.dC = .22856
        
        # Bond distances from C atoms
        self.dH = 1.1
        self.dM = 2.2
        
        self._2dCpos = []
        self._3dCpos = []
        
        self._2dHpos = []
        self._3dHpos = []
        
        self._2dMpos = []
        self._3dMpos = []
        
        self.pure = False
    
    def convert(self):
        """ Runs through the whole conversion process for a single POSCAR file. """
        uncleFile = open(self.uncleFile, 'r')
        uncleLines = uncleFile.readlines()
        uncleFile.close()
        
        # Extract the lattice vectors.
        vectorLines = [line.strip().split() for line in uncleLines[2:5]]
        for i in xrange(len(vectorLines)):
            for j in xrange(len(vectorLines[i])):
                if list(vectorLines[i][j])[0] == '*':
                    vectorLines[i][j] = '15.0'
                    break
            
        self.lattVec1 = [float(vectorLines[0][0]), float(vectorLines[0][1]), float(vectorLines[0][2])]
        self.lattVec2 = [float(vectorLines[1][0]), float(vectorLines[1][1]), float(vectorLines[1][2])]
        self.lattVec3 = [float(vectorLines[2][0]), float(vectorLines[2][1]), float(vectorLines[2][2])]
        
        # Get the number of each type of atom.
        nonCatomCounts = uncleLines[5].strip().split()
        nonCatomCounts = [int(count) for count in nonCatomCounts]
        
        self.atomCounts = [sum(nonCatomCounts)] + nonCatomCounts
        if sum([count != 0 for count in nonCatomCounts]) == 1:
            self.pure = True
                
        positionLines = uncleLines[7:]
        
        self.set3dCpositionsFromDirectCoordinates(positionLines)
        self.set3dpositionsFromDirectCoordinates(positionLines)
                      
    def get2DDistance(self, atom1, atom2):
        """ Returns the two-dimensional distance between two atoms. """
        xcomp = atom2[0] - atom1[0]
        ycomp = atom2[1] - atom1[1]
        
        return sqrt(pow(xcomp, 2) + pow(ycomp, 2))

    def set3dCpositionsFromDirectCoordinates(self, positionLines):
        """ Sets the 3D positions of the C atoms given the 2D 'surface' positions from UNCLE.
            Another way to think of this is that it introduces the "buckling" into the sheet of
            C atoms. """
        self._3dCpos = []
        for line in positionLines:
            position = line.strip().split()
            position = [float(comp) for comp in position]
            
            comp1 = [position[0] * self.lattVec1[0], position[0] * self.lattVec1[1], position[0] * self.lattVec1[2]]
            comp2 = [position[1] * self.lattVec2[0], position[1] * self.lattVec2[1], position[1] * self.lattVec2[2]]
            comp3 = [position[2] * self.lattVec3[0], position[2] * self.lattVec3[1], position[2] * self.lattVec3[2]]
            
            old_z = comp1[0] + comp2[0] + comp3[0]
            x = comp1[1] + comp2[1] + comp3[1]
            y = comp1[2] + comp2[2] + comp3[2]
            
            z = 0.0
            if old_z > 0:
                z = self.dC
            else:
                z = -self.dC
            
            self._3dCpos.append([x, y, z])
            
    def set3dpositionsFromDirectCoordinates(self, positionLines):
        """ Sets the 3D positions of the H atoms in the system. """
        self._3dpos = []
        for line in positionLines[0:]:
            position = line.strip().split()
            position = [float(comp) for comp in position]
            
            comp1 = [position[0] * self.lattVec1[0], position[0] * self.lattVec1[1], position[0] * self.lattVec1[2]]
            comp2 = [position[1] * self.lattVec2[0], position[1] * self.lattVec2[1], position[1] * self.lattVec2[2]]
            comp3 = [position[2] * self.lattVec3[0], position[2] * self.lattVec3[1], position[2] * self.lattVec3[2]]
            
            old_z = comp1[0] + comp2[0] + comp3[0]
            x = comp1[1] + comp2[1] + comp3[1]
            y = comp1[2] + comp2[2] + comp3[2]
            
            z = 0.0
            if old_z > 0:
                z = self.dC
            else:
                z = -self.dC
            
            self._3dpos.append([x, y, z])
        
    def getLatticeVectors(self):
        """ Returns the lattice vectors of the structure. """
        return [self.lattVec1, self.lattVec2, self.lattVec3]

    def getAtomCounts(self):
        """ Returns the number of each atom in the structure. """
        return self.atomCounts

    def getCPositions(self):
        """ Returns the list of C positions in the structure. """
        return self._3dCpos
    
    def getPositions(self):
        """ Returns the list of H positions in the structure. """
        return self._3dpos

    def isPure(self):
        """ Returns true if the structure is a pure structure, false otherwise. """
        return self.pure
    
    



    