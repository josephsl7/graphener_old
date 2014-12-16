'''
Created on Aug 29, 2014

@author: eswens13
'''
import os, subprocess
from matplotlib.patches import Patch
from pylab import *

class GSS:
    """ This class performs the UNCLE ground state search for the lowest formation energy
        structures. From this, we get a list of every structure that was enumerated and its
        corresponding formation energy as predicted by UNCLE. We use this list to decide which new
        structures to add into the model. We also keep track of the list and the plots of the list
        from iteration to iteration. """

    def __init__(self, atoms, volRange, plotTitle, xlabel, ylabel, uncleOutput):
        """ CONSTRUCTOR """
        
        self.atoms = atoms
        self.volRange = volRange
        self.plotTitle = plotTitle
        self.xlabel = xlabel
        self.ylabel = ylabel
        
        self.enumFolder = os.getcwd() + '/enum/'
        self.fitsDir = None
        self.neededFilesDir = os.getcwd() + '/needed_files/'
        self.uncleExec = os.getcwd() + '/needed_files/uncle.x'
        self.uncleOut = uncleOutput
        self.case = len(self.atoms[0].split(","))

        singleEFile = open(os.getcwd() + '/needed_files/single_atoms/single_atom_energies', 'r')
        inlines = singleEFile.readlines()
        singleEFile.close()

        self.singleE = []
        for line in inlines:
            self.singleE.append([line.split()[0], line.split()[2]])

        hexEFile = open(os.getcwd() + '/needed_files/hex_monolayer_refs/hex_energies', 'r')
        inlines = hexEFile.readlines()
        hexEFile.close()

        self.hexE = []
        for line in inlines:
            self.hexE.append([line.split()[0], line.split()[4]])
    
    def makeGSSDirectories(self):
        """ Creates the 'gss' directories for each different metal atom and populates them with 
            the files that UNCLE needs to perform a ground state search.  These files are 
            struct_enum.out, lat.in, J.1.out, groundstatesearch.in, and gss_plot.gp. """
        for atom in self.atoms:
            atomDir = os.getcwd() + '/' + atom
            if os.path.isdir(atomDir):
                self.fitsDir = atomDir + '/fits/'
                gssDir = atomDir + '/gss/'
                subprocess.call(['mkdir',gssDir])
                subprocess.call(['cp',self.enumFolder + 'struct_enum.out', gssDir])
                subprocess.call(['cp',self.enumFolder + 'lat.in', gssDir])
                subprocess.call(['cp',self.fitsDir + 'J.1.out', gssDir])
                
                infile = open(self.neededFilesDir + 'groundstatesearch.in','r')
                inlines = [line for line in infile]
                infile.close()
                
                outfile = open(gssDir + 'groundstatesearch.in','w')
                for i in xrange(len(inlines)):
                    if i == 4:
                        outfile.write(str(self.volRange[0]) + " " + str(self.volRange[1]) + "\n")
                    else:
                        outfile.write(inlines[i])
                outfile.close()
                
                infile = open(self.neededFilesDir + 'gss_plot.gp','r')
                inlines = [line for line in infile]
                infile.close()
                
                outfile = open(gssDir + 'gss_plot.gp','w')
                for i in xrange(len(inlines)):
                    if i == 3:
                        outfile.write("set xlabel \"" + self.xlabel + "\"\n")
                    elif i == 4:
                        outfile.write("set ylabel \"" + self.ylabel + "\"\n")
                    elif i == 5:
                        outfile.write("set title \"" + self.plotTitle + " (" + atom + ")\"\n")
                    elif i == 7:
                        outfile.write('plot "gss.out" u 1:2 lt 3, "vaspgraph.out" using 1:2 lt 1 pt 6\n')
                    else:
                        outfile.write(inlines[i])
                outfile.close()
    
    def performGroundStateSearch(self, iteration):
        """ Performs the ground state search with the current fit from UNCLE. """
        lastDir = os.getcwd()
        for atom in self.atoms:
            gssDir = os.getcwd() + '/' + atom + '/gss/'
            if os.path.isdir(gssDir):
                subprocess.call(['echo','\nPerforming ground state search for ' + atom + '. . .\n'])
                os.chdir(gssDir)
                subprocess.call([self.uncleExec, '21'], stdout=self.uncleOut)
                os.chdir(lastDir)

    def makePlots(self, iteration):
        """ Creates the plots of the predicted energies of all the structures that have been 
            enumerated. Adds the iteration number onto the end of the filenames for the plots and
            the lists. """
        lastDir = os.getcwd()
        for atom in self.atoms:
            gssDir = os.getcwd() + '/' + atom + '/gss/'
            if os.path.isdir(gssDir):
                subprocess.call(['echo', '\nMaking plot of ground states for ' + atom + '. . .\n'])
                os.chdir(gssDir)
                
                infile = open('gss.out','r')
                inlines = [line for line in infile]
                infile.close()

                subprocess.call(['mv','gss.out','gss_' + str(iteration) + '.out'])

                #Save VASP formation energies to a file
                structures = []

                structsfile = open('../fits/structures.in','r')
                structlines = structsfile.readlines()
                structsfile.close
                for line in structlines:
                    structures.append(line)

                structsfile = open('../fits/structures.holdout','r')
                structlines = structsfile.readlines()
                structsfile.close
                for line in structlines:
                    structures.append(line)

                if os.path.isdir('vaspFE.out'):
                    subprocess.call(['rm', 'vaspFE.out'])

                vaspfile = open('vaspFE.out', 'w')

                eIsolatedH = -1.115 
                eIsolatedC = -1.3179 
                eH2 = -6.7591696/2.0 
                energyGraphene = -18.456521 #for 2 C atoms 

                for lineNum, line in enumerate(structures):
                    if line.find('FE') != -1:
                        energyIndex = [i for i, x in enumerate(line.split()) if x == 'FE'][0] + 2
                        formEnergy = line.split()[energyIndex].strip(',')
                        atomCounts = [float(x) for x in structures[lineNum + 5].split()]
                        totalCounts = sum(atomCounts)

                        for searchLineNum, searchLine in enumerate(structures[lineNum + 7:]):
                            if searchLine.find('#Energy:') != -1:
                                energy = float(structures[lineNum + searchLineNum + 8])
                                break

                        bindEnergy = energy - energyGraphene / 2
                        for element, count in zip(atom.split(','), atomCounts):
                            if element == 'Vc':
                                pass
                            elif element == 'H':
                                bindEnergy -= count / totalCounts * eIsolatedH
                            else:
                                bindEnergy -= count / totalCounts * self.getSingleE(element)

                        hexEnergy = energy - energyGraphene / 2
                        for element, count in zip(atom.split(','), atomCounts):
                            if element == 'Vc':
                                pass
                            elif element == 'H':
                                hexEnergy -= count / totalCounts * eH2
                            else:
                                hexEnergy -= count / totalCounts * self.getHexE(element)

                        for count in atomCounts:
                            vaspfile.write(str(float(count) / sum(atomCounts)) + ' ')
                        vaspfile.write(formEnergy + ' ' + str(bindEnergy) + ' ' + str(hexEnergy) + '\n')
                vaspfile.close()

                vaspfile = open('vaspFE.out', 'r')
                vasplines = vaspfile.readlines()
                vaspfile.close()

                #Make the plots
                        
                for plotNum in range(1, self.case + 1):
                    if self.case == 2 and plotNum != 2:
                        continue

                    outfile = open('gss.out','w')
                    for line in inlines[2:]:
                        if float(line.split()[plotNum]) == 0.0 or self.case == 2:
                            outfile.write(line.split()[1 + (plotNum == 1)] + ' ' + line.split()[3 + 2*self.case] + '\n')
                    outfile.close()

                    vaspoutfile = open('vaspgraph.out','w')
                    for line in vasplines:
                        if float(line.split()[plotNum - 1]) == 0.0 or self.case == 2:
                            vaspoutfile.write(line.split()[0 + (plotNum == 1)] + ' ' + line.split()[self.case] + '\n')
                    vaspoutfile.close()

                    gssinfile = open('gss_plot.gp','r')
                    gssinlines = [line for line in gssinfile]
                    gssinfile.close()

                    gssoutfile = open('gss_plot.gp','w')
                    for i in xrange(len(gssinlines)):
                        if i == 3 and plotNum == 1:
                            gssoutfile.write("set xlabel \"" + atom.split(',')[1] + " Concentration\"\n")
                        elif i == 3 and plotNum != 1:
                            gssoutfile.write("set xlabel \"" + atom.split(',')[0] + " Concentration\"\n")
                        elif i == 5 and self.case != 2:
                            gssoutfile.write("set title \"" + self.plotTitle + " (" + atom + ") No " + atom.split(',')[plotNum - 1] + "\"\n")
                        else:
                            gssoutfile.write(gssinlines[i])
                    gssoutfile.close()

                    subprocess.call(['gnuplot', 'gss_plot.gp'])
                    if self.case == 2:
                        subprocess.call(['mv','gss.pdf','gss_' + str(iteration) + '.pdf'])
                    else:
                        elements = '_'.join([x.split('_')[0] for i, x in enumerate(atom.split(',')) if i+1 != plotNum])
                        subprocess.call(['mv','gss.pdf','gss_' + str(iteration) + '_' + elements + '.pdf'])

                    subprocess.call(['rm','gss.out'])
                    subprocess.call(['rm','vaspgraph.out'])

                #Ternary plots
                if self.case == 3:
                    #GSS data
                    ternaryPoints = []
                    for line in inlines[2:]:
                        x1 = float(line.split()[1])
                        x2 = float(line.split()[2])
                        x3 = float(line.split()[3])
                        ternaryPoints.append([float(line.split()[9]),(x3+2*x1)/(2*(x1+x2+x3)),1.7320508*x3/(2*(x1+x2+x3))])
                    ternaryPoints.sort()
                    ternaryPoints = ternaryPoints[::-1]
                    xdata = [x[1] for x in ternaryPoints]
                    ydata = [x[2] for x in ternaryPoints]
                    energydata = [x[0] for x in ternaryPoints]
                    self.ternaryPlot(xdata, ydata, energydata, atom, "UNCLE Formation Energies Per Concentration (eV)", 'gss_' + str(iteration) + '_ternary.pdf')

                    #VASP data
                    ternaryPoints = []
                    for line in vasplines:
                        x1 = float(line.split()[0])
                        x2 = float(line.split()[1])
                        x3 = float(line.split()[2])      
                        ternaryPoints.append([float(line.split()[3]),(x3+2*x1)/(2*(x1+x2+x3)),1.7320508*x3/(2*(x1+x2+x3))])        
                    ternaryPoints.sort()
                    ternaryPoints = ternaryPoints[::-1]
                    xdata = [x[1] for x in ternaryPoints]
                    ydata = [x[2] for x in ternaryPoints]
                    energydata = [x[0] for x in ternaryPoints]
                    self.ternaryPlot(xdata, ydata, energydata, atom, "VASP Formation Energies Per Concentration (eV)", 'vasp_' + str(iteration) + '_ternary.pdf')

                    #VASP binding energies
                    ternaryPoints = []
                    for line in vasplines:
                        x1 = float(line.split()[0])
                        x2 = float(line.split()[1])
                        x3 = float(line.split()[2])      
                        ternaryPoints.append([float(line.split()[4]),(x3+2*x1)/(2*(x1+x2+x3)),1.7320508*x3/(2*(x1+x2+x3))])        
                    ternaryPoints.sort()
                    ternaryPoints = ternaryPoints[::-1]
                    xdata = [x[1] for x in ternaryPoints]
                    ydata = [x[2] for x in ternaryPoints]
                    energydata = [x[0] for x in ternaryPoints]
                    self.ternaryPlot(xdata, ydata, energydata, atom, "VASP Binding Energies Per Concentration (eV)", 'vasp_' + str(iteration) + '_bind.pdf')

                    #VASP hex energies
                    ternaryPoints = []
                    for line in vasplines:
                        x1 = float(line.split()[0])
                        x2 = float(line.split()[1])
                        x3 = float(line.split()[2])      
                        ternaryPoints.append([float(line.split()[5]),(x3+2*x1)/(2*(x1+x2+x3)),1.7320508*x3/(2*(x1+x2+x3))])        
                    ternaryPoints.sort()
                    ternaryPoints = ternaryPoints[::-1]
                    xdata = [x[1] for x in ternaryPoints]
                    ydata = [x[2] for x in ternaryPoints]
                    energydata = [x[0] for x in ternaryPoints]
                    self.ternaryPlot(xdata, ydata, energydata, atom, "VASP Hex Energies Per Concentration (eV)", 'vasp_' + str(iteration) + '_hex.pdf')

                os.chdir(lastDir)

    def ternaryPlot(self, xdata, ydata, energydata, atom, plotTitle, outfile):
        fig = figure(figsize=[12,8])
        if plotTitle.find('Formation') != -1:
            scatter(xdata, ydata, s=35, c=energydata, cmap=cm.hot, linewidths=0, vmax = 0.2)
        else:
            scatter(xdata, ydata, s=35, c=energydata, cmap=cm.hot, linewidths=0)
        colorbar()
        title(plotTitle, fontsize=20)
        xticks([])
        yticks([])
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.85)
        ax.text(-0.08,-0.06,atom.split(',')[1].split("_")[0],fontsize=20)
        ax.text(1.025,-0.06,atom.split(',')[0].split("_")[0],fontsize=20)
        ax.text(0.467,0.9,atom.split(',')[2].split("_")[0],fontsize=20)
        fig.savefig(outfile)
    
    def contains(self, struct, alist):
        """ Returns true if 'struct' is found in 'alist', false otherwise. """
        for i in xrange(len(alist)):
            if str(struct) == str(alist[i]):
                return True
    
        return False

    def getSingleE(self, atom):
        for lineAtom, energy in self.singleE:
            if lineAtom == atom:
                return float(energy)
        return 100

    def getHexE(self, atom):
        for lineAtom, energy in self.hexE:
            if lineAtom == atom:
                return float(energy)
        return 100
    
    def getAllGSSStructures(self, iteration, failedStructs):
        """ Returns a list of all the structures sorted by their predicted formation energies.
            It does this for each metal atom that has been specified by the user so this will 
            actually return a list of lists--a list for each atom. """
        allStructs = []
        for n in xrange(len(self.atoms)):
            atomStructs = []
            structsByEnergy = []
            gssFile = os.getcwd() + '/' + self.atoms[n] + '/gss/gss_' + str(iteration) + '.out'
            infile = open(gssFile, 'r')
            
            i = 0
            for line in infile:
                if i >= 2:
                    formEnergy = float(line.strip().split()[3 + 2*self.case])
                    struct = int(line.strip().split()[0])
                    
                    # Do not include structures that failed VASP calculations.
                    if not self.contains(struct, failedStructs[n]):
                        structsByEnergy.append([formEnergy, struct])
                i += 1
            infile.close()
            
            structsByEnergy.sort()
            
            for struct in structsByEnergy:
                atomStructs.append(str(struct[1]))
            
            allStructs.append(atomStructs)
            
        return allStructs
            
                







        
