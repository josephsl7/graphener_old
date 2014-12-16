'''
Created on Aug 27, 2014

@author: eswens13
'''
import os, subprocess

class Extractor:
    """ This class is responsible for creating "pseudo-POSCAR" files for each structure in the
        set of structures that we want to run through VASP calculations.  It does this by 
        "extracting" the information about the structure from struct_enum.out using the makestr.x
        routine from the enumlib library in UNCLE.  The Structs2Poscar class will then take this 
        set of pseudo-POSCARs and prepare them for VASP calculations. """

    def __init__(self, atoms, uncleOutput, case):
        """ CONSTRUCTOR """
        
        self.atoms = atoms
        self.extractExec = os.path.abspath('needed_files/makestr.x')
        self.uncleOut = uncleOutput
        self.structList = []
        self.case = case

    def setTrainingStructs(self):
        """ This method sets the structure list for each atom to the list of training structures
            generated by UNCLE i.i.d. selection.  It DOES NOT choose the set of traning structures,
            it just sets the member list of this class.  We always include the pure structures in 
            the list. """
        self.structList = []
        trainStructs = []
        trainFile = open('enum/training_set_structures.dat','r')
        for line in trainFile:
            trainStructs.append(line.strip().split()[1])
        trainFile.close()
        
        # Make sure the pure structures are in the training structures list.
        struct = 1        
        for i, nextPureCase in enumerate(range(self.case,0,-1)):
            if not self.contains(str(struct), trainStructs):
                trainStructs[i] = str(struct)
            struct = struct + nextPureCase
        
        for atom in self.atoms:
            self.structList.append(trainStructs)
    
    def contains(self, struct, alist):
        """ Returns True if 'alist' contains the item 'struct'. """
        if len(alist) == 0:
            return False
        
        for i in xrange(len(alist)):
            if struct == alist[i]:
                return True
        
        return False
    
    def setStructList(self, alist):
        """ Sets the list of structures that we want to run through VASP calculations.  The list 
            being passed to this method (alist) will actually be a list of lists.  It will have a 
            structure list for each atom that still has not finished the main convergence loop. """
        self.structList = []
        for atomStructs in alist:
            self.structList.append(atomStructs)
    
    def getStructList(self):
        """ Returns the list of structures for each atom that has not finished the convergence
            loop. """
        structs = []
        for atomStructs in self.structList:
            structs.append(atomStructs)
        
        return structs
          
    def extract(self):
        """ You must call either the setTrainingStructs() or setStructList() functions before calling
            this method.  This method uses the makestr.x executable from the enumlib in UNCLE to 
            create the pseudo-POSCAR files for each structure in self.structList. These files are 
            generally called something like "vasp.000241" indicating the structure number in 
            struct_enum.out.  We only want to extract the union of all the lists in 
            self.structList. """
        subprocess.call(['echo','\nExtracting structures from struct_enum.out\n'])
        lastDir = os.getcwd()
        os.chdir(lastDir + '/enum')
        uniqueSet = set(self.structList[0])
        
        # Only extract the union of all the sets of structures.  (No duplicates)
        for i in xrange(1,len(self.structList)):
            uniqueSet = uniqueSet.union(self.structList[i])
        
        for struct in uniqueSet:
            subprocess.call([self.extractExec, 'struct_enum.out', struct], stdout=self.uncleOut)
       
        os.chdir(lastDir)
        

        