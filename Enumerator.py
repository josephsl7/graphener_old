'''
Created on Aug 26, 2014

@author: eswens13
'''
import os, subprocess, shutil


class Enumerator:
    """ This class enumerates symmetrically unique structures in a given volume range using UNCLE.  
        It then builds the clusters necessary to perform a cluster expansion and chooses a 
        specified number of "training structures" to perform a first fit on.  After this class 
        finishes its work, the Extractor class will take over and extract pseudo-POSCAR files from 
        the struct_enum.out file that is produced. The methods in this class are only needed for 
        the first iteration of the main convergence loop. """
  
    def __init__(self, atoms, volRange, clusterNums, trainStructNum, uncleOutput, case):
        """ CONSTRUCTOR """
        self.atoms = atoms
        self.volRange = volRange
        
        self.clusterNums = clusterNums
        self.trainStructNum = trainStructNum
        
        self.uncleExec = os.path.abspath('needed_files/uncle.x')
        self.enumFile = 'enum/struct_enum.out'
        self.enumExec = os.path.abspath('needed_files/enum.x')
        self.uncleOut = uncleOutput
        self.case = case

    def changeEnumFile(self):
        """ In order to build the clusters that will be used in the cluster expansion correctly, 
            we have to change the 'surf' setting in struct_enum.out (from UNCLE enumeration) to 
            'bulk'.  It changes the name of the old 'surf' version to 'struct_enum.out_OLD'. """
        subprocess.call(['mv',self.enumFile, self.enumFile + '_OLD'])
        
        oldfile = open(self.enumFile + '_OLD','r')
        oldlines = [line for line in oldfile]
        oldfile.close()
        
        newfile = open(self.enumFile, 'w')
        for i in xrange(len(oldlines)):
            if i == 1:
                newfile.write('bulk\n')
            else:
                newfile.write(oldlines[i])
        
        newfile.close()
    
    def buildClusters(self):
        """ Uses UNCLE to build the number of each n-body clusters specified in the settings.in
            file. """
        oldLatFile = 'needed_files/lat.in'
        oldFile = open(oldLatFile, 'r')
        oldLines = [line for line in oldFile]
        oldFile.close()
        
        newFile = open('enum/lat.in','w')
        for i in xrange(len(oldLines)):
            if i == 42:
                for num in self.clusterNums:
                    newFile.write(str(num) + " ")
                newFile.write("\n")
            elif i == 13:
                newFile.write(str(self.case) + "  #case\n")
            else:
                newFile.write(oldLines[i])
        newFile.close()
        
        lastDir = os.getcwd()
        os.chdir(lastDir + '/enum')
        
        subprocess.call([self.uncleExec, '10'], stdout=self.uncleOut)
        #subprocess.call(['cp', 'struct_enum.out', '../needed_files/struct_enum.out'])
        
        os.chdir(lastDir)

    def chooseTrainingStructures(self):
        """ Chooses a list of i.i.d. structures from struct_enum.out. The length of the list 
            is determined by the TRAINING_STRUCTS setting in settings.in. """
        lastDir = os.getcwd()
        os.chdir(lastDir + '/enum')
        
        subprocess.call([self.uncleExec, '42', str(self.trainStructNum)], stdout=self.uncleOut)
        
        os.chdir(lastDir)
    
    def enumerate(self):
        """ Runs through the whole process of enumeration, cluster building, and choosing an
            i.i.d. set of training structures. """
        if os.path.isdir('enum'):
            shutil.rmtree('enum')
        subprocess.call(['mkdir','enum'])

        infile = open('needed_files/struct_enum.in','r')
        inlines = []
        for line in infile:
            inlines.append(line)
        infile.close()
        
        structFile = open('enum/struct_enum.in','w')
        for i in xrange(len(inlines)):
            if i == 9:
                structFile.write(str(self.volRange[0]) + " " + str(self.volRange[1]) + " ")
                structFile.write("# Starting and ending cell sizes for search\n")
            elif i == 5:
                structFile.write(" " + str(self.case) + " -nary case\n")
            elif i == 7 or i == 8:
                structFile.write(inlines[i].split()[0]+" ")
                structFile.write(inlines[i].split()[1]+" ")
                structFile.write(inlines[i].split()[2]+"    ")
                for num in range(self.case):
                    structFile.write(str(num))
                    if num < self.case-1:
                        structFile.write("/")
                structFile.write("   # d0" + str(i-6) + " d-vector\n")
            else:
                structFile.write(inlines[i])
        structFile.close()
        
        lastDir = os.getcwd()
        os.chdir(lastDir + '/enum')
        subprocess.call([self.enumExec,'struct_enum.in'], stdout=self.uncleOut)
        
        os.chdir(lastDir)
        
        self.changeEnumFile()
        subprocess.call(['echo','\nGenerating clusters. . .\n'])
        self.buildClusters()
        subprocess.call(['echo','\nChoosing i.i.d. structures. . .\n'])
        self.chooseTrainingStructures()
        
            
            
        
