import os, subprocess

class RunVasp:
    """ This class is responsible for preparing the directories, retrieving the needed files, and
        submitting VASP jobs to the supercomputer.  It keeps track of the SLURM job ids of all the
        jobs that are currently running from a particular instance of the class. """
    
    def __init__(self, atomList):
        """ CONSTRUCTOR """
        
        self.atomList = atomList
        self.neededFilesFolder = os.getcwd() + '/needed_files/'
        
        self.currJobIds = []
    
    def prepareForVasp(self, structList):
        """ Makes all of the files that could be copied to a first, low-precision VASP run for any 
            given structure.  This includes concatenating the POTCARS for the pure and non-pure
            cases. """
        self.makeLowINCARs()
        self.makePurePOTCARs()
        self.makePOTCARs()
        self.makeKPOINTS(8, 8)
        self.makeJobFiles()
        self.copyVaspExec()
        
        self.fillDirectories(structList)

    def makeLowINCARs(self):
        """ Creates a standard INCAR file and puts it in each different structure's top 
            directory. """
        dirList = self.atomList
        
        for direc in dirList:
            incar = open(direc + '/INCAR','w')
    
            incar.write("IBRION=2\n")
            incar.write("ISIF=4\n")
            incar.write("NSW=400\n")
            incar.write("Algo=VeryFast\n")
            incar.write("PREC=Low\n") #bch NORMAL for first and only run...was Low!
            incar.write("EDIFF=5E-4\n")
            incar.write("EDIFFG=5E-4\n")
            incar.write("ISMEAR=0\n")
            incar.write("ISPIN=2\n")
            incar.write("LREAL=Auto\n") #Usually "LREAL=Auto\n" but FALSE for very small cells
            incar.write("SIGMA=0.1\n")
            incar.write("LWAVE=.TRUE.\n")
            incar.write("LCHARG=.TRUE.\n")
    
            incar.close()
  
    def makePurePOTCARs(self):
        """ Some of the structures that need to be submitted to VASP for relaxation are what we 
            call "pure" structures.  This means that (other than carbon atoms) the structure only 
            contains one other type of atom.  In the binary representation of the structure, this 
            means that it is either all '1's or all '0's.  VASP gives a segmentation fault when we 
            try to run a pure structure with a POTCAR that contains atoms that are not in the pure 
            structure, even if we tell it that there are zero of one of the kinds of atoms.  It 
            needs a POTCAR file that doesn't even mention the element that is not a part of the 
            "pure" structure. This method creates these POTCAR files. """       
        for atom in self.atomList:
            atomPotcarDir = "/fslhome/josephsl/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/C"
            
            if os.path.isdir(atomPotcarDir):
                purePotcar = open(atom + "/POTCAR_C", "w")
                
                CPotcar = open(atomPotcarDir + "/POTCAR", "r")
                CLines = CPotcar.readlines()
                CPotcar.close()
             
                for line in CLines:
                    purePotcar.write(line)
                
                purePotcar.close()

    def makePOTCARs(self):
        """ Creates a POTCAR file for each atom in the member 'atomList'. Concatenates the 
            individual POTCAR files to make a single POTCAR file for the multi-atom structure. """
        for atomgroup in self.atomList:
            atoms = atomgroup.split(',')
            for atom in atoms:
                atomPotcarDir = "/fslhome/josephsl/fsl_groups/hessgroup/vaspfiles/src/potpaw_PBE/" + atom + "/POTCAR"
                if os.path.exists(atomPotcarDir):
                    atomPotcar = open(atomPotcarDir,'r')
                    atomLines = atomPotcar.readlines()
                    atomPotcar.close()
            
                    potcar = open(atomgroup + '/POTCAR_' + atom, 'w')
     
                    for line in atomLines:
                        potcar.write(line)

                    potcar.close()
                        
                elif atom.find('Vc') == -1:
                    subprocess.call(['echo','ERROR: Could not find a POTCAR file for \'' + atom + '\''])
                    subprocess.call(['echo','Removing POTCAR . . .'])
                    subprocess.call(['rm','POTCAR'])
                    return

    def makeKPOINTS(self, num1, num2):
        """ Creates a KPOINTS file based on the input parameters num1 and num2. It specifies that 
            the job will have num1 x num2 kpoints. For example, if we wanted to specify an 8x8 
            kpoints mesh, we would call makeKPOINTS(8, 8). """
        dirList = self.atomList
        
        for direc in dirList:
            kpoints = open(direc + '/KPOINTS','w')
    
            kpoints.write("Automatic mesh\n")
            kpoints.write("0\n")
            kpoints.write("Gamma\n")
            kpoints.write(str(num1) + ' ' + str(num2) + ' 1\n')
            kpoints.write('0 0 0')
    
            kpoints.close()

    def makeJobFiles(self):
        """ Creates a standard job file for submitting a VASP job on the supercomputer. """
        dirList = self.atomList
        
        for direc in dirList:
            jobFile = open(direc + '/job','w')
    
            jobFile.write("#!/bin/bash\n\n")
            jobFile.write("#SBATCH --time=01:00:00\n")
            jobFile.write("#SBATCH --ntasks=16\n")  # Don't change this!
            jobFile.write("#SBATCH --mem-per-cpu=1024M\n")  # Don't change this!
            jobFile.write("#SBATCH --mail-user=hess.byu@gmail.com\n")
            jobFile.write("#SBATCH --mail-type=FAIL\n")
            jobFile.write("\nmpiexec vasp533 > vasp.out\n")
    
            jobFile.close()

    def copyVaspExec(self):
        """ Copies the vasp executable file to the current working directory. """
        dirList = self.atomList
        
        for direc in dirList:
            subprocess.call(['cp','/fslhome/josephsl/bin/vasp', direc + '/vasp533'])
            subprocess.call(['chmod','777', direc + '/vasp533']) #This makes the file executable

    def fillDirectories(self, structList):
        """ Fills all the directories in 'structList' with the needed files for VASP to run, namely
            POSCAR, POTCAR, KPOINTS, INCAR, a SLURM job file, and the VASP executable file. """
        for i in xrange(len(self.atomList)):
            lastDir = os.getcwd()
            atomDir = lastDir + '/' + self.atomList[i]
            
            os.chdir(atomDir)
            structures = []
            for item in structList[i]:
                if os.path.isdir(item):
                    structures.append(item)

            atoms = self.atomList[i].split(',')
            
            for structure in structures:
                structureDir = os.path.abspath(structure)
                subprocess.call(['cp', 'KPOINTS', 'INCAR', 'job', 'vasp533', structureDir])
                
                poscar = open(structureDir + '/POSCAR','r')
                poscarLines = [line.strip() for line in poscar]
                poscar.close()
                
                line = poscarLines[0]
                present = line[line.find(':')+1:line.find('-')].strip().split()

                potcar = open(structureDir + '/POTCAR', 'w')

                atomPotcarDir = atomDir + '/POTCAR_C'

                if os.path.exists(atomPotcarDir):
                    atomPotcar = open(atomPotcarDir,'r')
                    atomLines = atomPotcar.readlines()
                    atomPotcar.close()                    
                
                    for writeline in atomLines:
                        potcar.write(writeline)
                else:
                    subprocess.call(['echo','ERROR: Could not find '+ atomPotcarDir + '\'']) 

                for atomnum in present:
                    atomPotcarDir = atomDir + '/POTCAR_' + atoms[int(atomnum)-1]

                    if os.path.exists(atomPotcarDir):
                        atomPotcar = open(atomPotcarDir,'r')
                        atomLines = atomPotcar.readlines()
                        atomPotcar.close()                    
                
                        for writeline in atomLines:
                            potcar.write(writeline)
                    else:
                        subprocess.call(['echo','ERROR: Could not find '+ atomPotcarDir + '\'']) 
                        
                potcar.close()
            
            os.chdir(lastDir)
    
    def run(self, runNum, structList):
        """ Starts the VASP runs (specified by 'runNum') for each of the structures in
            'structList'. For runNum = 1, starts a low-precision run, runNum = 2, starts a 
            normal-precision run, runNum = 3 starts a DOS run. """
        if runNum == 1:
            self.startJobs(structList)
    
        elif runNum == 2:
            self.makeNormalDirectories(structList)
            self.startNormalJobs(structList)
           
        elif runNum == 3:
            self.makeDOSDirectories(structList)
            self.startDOSJobs(structList)

    def startJobs(self, structList):
        """ Submits all the VASP jobs for structures in 'structList' to the supercomputer for 
            low-precision relaxation and records their job IDs. """
        self.clearCurrentJobIds()
        
        for i in xrange(len(self.atomList)):
            lastDir = os.getcwd()
            atomDir = lastDir + '/' + self.atomList[i]
            
            os.chdir(atomDir)
            
            structures = []
            for item in structList[i]:
                if os.path.isdir(item):
                    structures.append(item)
            
            for structure in structures:
                os.chdir(structure)
                proc = subprocess.Popen(['sbatch','job'], stdout=subprocess.PIPE)
                jobid = proc.communicate()[0].split()[3]
                subprocess.call(['echo', 'Submitted job ' + jobid])
                self.currJobIds.append(jobid)
                os.chdir(atomDir)
            
            os.chdir(lastDir)

    def clearCurrentJobIds(self):
        """ Clears the list of current job IDs. """
        self.currJobIds = []

    def makeNormalDirectories(self, structList):
        """ After a structure has finished low-precision VASP relaxation, makes a directory for 
            normal-precision relaxation and populates it with the files from the low-precision run.
            Copies the low CONTCAR to the normal POSCAR. """
        topDir = os.getcwd()
        for i in xrange(len(self.atomList)):
            elementDir = topDir + '/' + self.atomList[i]
            if os.path.isdir(elementDir):
                os.chdir(elementDir)
                for structure in structList[i]:
                    structDir = elementDir + '/' + structure
                    if os.path.isdir(structDir) and self.finishCheck(structDir) and self.convergeCheck(structDir, 400):
                        os.chdir(structDir)
                        subprocess.call(['mkdir', 'normal'])
                        subprocess.call(['cp','CONTCAR','DOSCAR','EIGENVAL',
                                         'IBZKPT','KPOINTS','vasp533',
                                         'OSZICAR','OUTCAR','PCDAT',
                                         'POSCAR','POTCAR','REPORT','struct',
                                         'vasprun.xml','job','XDATCAR','normal'])
                        self.makeNormalINCAR()
                        subprocess.call(['cp','normal/CONTCAR','normal/POSCAR'])
                        os.chdir(elementDir)
            else:
                subprocess.call(['echo','The directory ' + elementDir + ' does not exist.'])
            
            os.chdir(topDir)      

    def makeNormalINCAR(self):
        """ Creates a standard INCAR file for normal-precision relaxation. """ 
        incar = open('normal/INCAR','w')
    
        incar.write("IBRION=2\n")
        incar.write("ISIF=4\n")
        incar.write("NSW=400\n")
        incar.write("Algo=VeryFast\n")
        incar.write("PREC=Normal\n")
        incar.write("EDIFF=5E-4\n")
        incar.write("EDIFFG=5E-4\n")
        incar.write("ISMEAR=0\n")
        incar.write("ISPIN=2\n")
        incar.write("LREAL=Auto\n") #Usually "LREAL=Auto\n" but FALSE for very small cells
        incar.write("SIGMA=0.1\n")
        incar.write("LWAVE=.TRUE.\n")
        incar.write("LCHARG=.TRUE.\n")
    
        incar.close()

    def startNormalJobs(self, structList):
        """ Submits all the VASP jobs for structures in 'structList' to the supercomputer for 
            normal-precision relaxation and records their job IDs. """
        self.clearCurrentJobIds()
        
        for i in xrange(len(self.atomList)):
            lastDir = os.getcwd()
            atomDir = lastDir + '/' + self.atomList[i]
            
            os.chdir(atomDir)
            
            structures = []
            for item in structList[i]:
                if os.path.isdir(item + '/normal'):
                    structures.append(item)
            
            for structure in structures:
                os.chdir(structure + '/normal')
                proc = subprocess.Popen(['sbatch','job'], stdout=subprocess.PIPE)
                jobid = proc.communicate()[0].split()[3]
                subprocess.call(['echo','Submitted job ' + jobid])
                self.currJobIds.append(jobid)
                os.chdir(atomDir)
            
            os.chdir(lastDir)

    def makeDOSDirectories(self, structList):
        """ After the normal-precision relaxation, creates a directory for the Density of States
            run and populates it with the files from the normal-precision run. Copies the normal
            CONTCAR to the DOS POSCAR. """  
        topDir = os.getcwd()
        for i in xrange(len(self.atomList)):
            elementDir = topDir + '/' + self.atomList[i]
            if os.path.isdir(elementDir):
                os.chdir(elementDir)
                for structure in structList[i]:
                    structDir = elementDir + '/' + structure
                    if os.path.isdir(structDir):
                        normalDir = structDir + '/normal'
                        if os.path.isdir(normalDir) and self.finishCheck(normalDir) and self.convergeCheck(normalDir, 400):
                            os.chdir(structDir)
                            subprocess.call(['mkdir', 'DOS'])
                            subprocess.call(['cp','normal/CONTCAR','normal/DOSCAR','normal/EIGENVAL',
                                             'normal/IBZKPT','normal/KPOINTS','normal/vasp533',
                                             'normal/OSZICAR','normal/OUTCAR','normal/PCDAT',
                                             'normal/POSCAR','normal/POTCAR','normal/REPORT','normal/struct',
                                             'normal/vasprun.xml','normal/XDATCAR','DOS'])
                            self.makeDOS_INCAR()
                            self.makeDOSJobFile()
                            subprocess.call(['cp','DOS/CONTCAR','DOS/POSCAR'])
                            os.chdir(elementDir)
            else:
                subprocess.call(['echo','The directory ' + elementDir + ' does not exist.'])
            
            os.chdir(topDir)

    def makeDOS_INCAR(self):
        """ Creates an INCAR file for the Density of States run in VASP.  The notable changes are:
                IBRION=-1 -- This tells VASP not to move the ions.
                NSW=0     -- This tells VASP that there will be no ionic relaxation steps.
                LORBIT=10 -- This creates the PROCAR file which can be used to project onto the 
                             C, H, and M atoms. """
                
        incar = open('DOS/INCAR','w')
        
        incar.write("IBRION=-1\n")
        incar.write("ISIF=4\n")
        incar.write("NSW=0\n")
        incar.write("Algo=Normal\n")
        incar.write("PREC=Normal\n")
        incar.write("EDIFF=5E-4\n")
        incar.write("EDIFFG=5E-4\n")
        incar.write("ISMEAR=0\n")
        incar.write("ISPIN=2\n")
        incar.write("LREAL=Auto\n") #Usually "LREAL=Auto\n" but FALSE for very small cells
        incar.write("LORBIT=10\n")
        incar.write("SIGMA=0.1\n")
        incar.write("LWAVE=.TRUE.\n")
        incar.write("LCHARG=.TRUE.\n")
    
        incar.close()   

    def makeDOSJobFile(self):
        """ Creates a SLURM job file for a VASP Density of States run. """
        jobFile = open('DOS/job','w')
        
        jobFile.write("#!/bin/bash\n\n")
        jobFile.write("#SBATCH --time=01:00:00\n")
        jobFile.write("#SBATCH --ntasks=16\n")  # Don't change this!
        jobFile.write("#SBATCH --mem-per-cpu=1024M\n")  # Don't change this!
        jobFile.write("#SBATCH --mail-user=hess.byu@gmail.com\n")
        jobFile.write("#SBATCH --mail-type=FAIL\n")
        jobFile.write("\nmpiexec vasp533 > vasp.out\n")
    
        jobFile.close()

    def startDOSJobs(self, structList):
        """ Submits all the VASP jobs for structures in 'structList' to the supercomputer for 
            Density of States calculations. Records their SLURM job IDs. """
        topDir = os.getcwd()
        self.clearCurrentJobIds()
        
        for i in xrange(len(self.atomList)):
            elementDir = topDir + '/' + self.atomList[i]
            if os.path.isdir(elementDir):
                os.chdir(elementDir)
                
                for structure in structList[i]:
                    structDir = elementDir + '/' + structure
                    if os.path.isdir(structDir):
                        os.chdir(structDir)
                        
                        dosDir = structDir + '/DOS'
                        if os.path.isdir(dosDir):
                            os.chdir(dosDir)
                            proc = subprocess.Popen(['sbatch','job'], stdout=subprocess.PIPE)
                            jobid = proc.communicate()[0].split()[3]
                            subprocess.call(['echo','Submitted job ' + jobid])
                            self.currJobIds.append(jobid)
                        
                        os.chdir(structDir)
                    
                    os.chdir(elementDir)
            else:
                subprocess.call(['echo','The directory ' + elementDir + ' does not exist.'])
            
            os.chdir(topDir)

    def getCurrentJobIds(self):
        """ Returns the list of current job IDs. """
        return self.currJobIds

    def hasFinished(self, folder):
        """ Determines whether the structure has finished and converged VASP calculations by 
            looking at its folder. """
        if self.finishCheck(folder) and self.convergeCheck(folder, 400):
            return True
        else:
            return False
        
    def finishCheck(self, folder):
        """ Tests whether VASP is done by finding "Voluntary" in last line of OUTCAR.  The input
            parameter, folder, is the directory containing OUTCAR, not the OUTCAR file itself. """
        lastfolder = os.getcwd()
        os.chdir(os.path.abspath(folder))
        proc = subprocess.Popen(['grep', 'Voluntary', 'OUTCAR'],stdout=subprocess.PIPE)
        newstring = proc.communicate()
        os.chdir(lastfolder)    
        return newstring[0].find('Voluntary') > -1 #True/False
    
    def convergeCheck(self, folder, NSW):
        """ Tests whether force convergence is done by whether the last line of OSZICAR (the total
            number of ionic relaxation steps) is less than NSW."""
        try:
            value = self.getSteps(folder)
            return value < NSW #True/False
        except:
            return False #True/False

    def getSteps(self, folder):
        """ Returns the number of ionic relaxation steps the structure has gone through as an
            integer. """
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








