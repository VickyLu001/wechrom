from email import header
import pkg_resources
import os
import numpy as np
from collections import defaultdict
from tqdm import tqdm
from .prepare_memory import generate_naked_memory_files, get_naked_memory_template_bonds
from .utils import num_intra_memory_bonds, num_inter_memory_bonds
from .utils import _DNA_S1_TYPE, _DNA_S2_TYPE, _DNA_RES

try:
    # for openmm version >= 7.6
    from openmm import LangevinIntegrator, HarmonicBondForce, CustomNonbondedForce, CustomBondForce, Platform
    from openmm.app import CheckpointReporter, Simulation, DCDReporter, PDBxFile, ForceField, DCDFile, PDBFile
    from openmm.unit import *
except:
    # for openmm version <= 7.5
    from simtk.openmm import LangevinIntegrator, HarmonicBondForce, CustomNonbondedForce, CustomBondForce, Platform
    from simtk.openmm.app import Simulation, DCDReporter, CheckpointReporter, DCDFile, PDBxFile, ForceField, PDBFile
    from simtk.unit import *


_DNA_XML = pkg_resources.resource_filename(__name__, 'data/DNA.xml')

class WechromSystem:
    """ Build a WEChroM system with topology and system information for openmm simulations
    The vanilla version should only have DNA beads
    
    Attributes:
        positions (np.array, N x 3):
            positions of all atoms in the initial pdbx file
        atoms (list of openmm Atoms):
            index of all atoms (DNA)
        residues (list of openmm Residues):
            index of all residues (base pairs)
        dnaStrands (dict of list of int):
            index of all atoms in the two strands I and J
        dnaRes (list of int):
            index of all DNA residues
        nBp (int):
            number of DNA residues
        dna_Atom2Res, dna_Res2Atom (dict of dict):
            atom2res and res2atom relation for the two dna strands
        dnaTypes (list):
            name of the two strands I and J
        pdbx (openmm PDBxFile):
        forcefield (openmm forcefield):
            built on the xml file
        system (openmm system):
            built on topology from pdbx and forcefield

    """
    def __init__(self, pdbxFile, xmlFile='DEFAULT',  memoryRange = 4, verbose = False):
        """Constructs all the necessary attributes to a WEChroM system

        Args:
            pdbxFile (str or file): the pdbx file containing system topology and initial positions
            xmlFile (str or file, optional): the xml file containing openmm forcefield. Defaults to data/DNA.xml.
            verbose(bool, optional): whether to print the information during the simulation. Defaults to True
        """
        # Construct a Topology and atom positions from the pdbx file
        
        # allow pdb filename beyond pdbx/mmcif file
        if isinstance(pdbxFile, str) and pdbxFile.endswith("pdb"):
            self.pdbx = PDBFile(pdbxFile)
        else:
            self.pdbx = PDBxFile(pdbxFile)
        self.topology = self.pdbx.topology
        self.positions = self.pdbx.getPositions(asNumpy=True)
        self.atoms = list(self.topology.atoms())
        self.residues = list(self.topology.residues())
        self.memoryRange = memoryRange

        # Read the xml file to construct a system
        if xmlFile == 'DEFAULT':
            xmlFile = _DNA_XML
        self.forcefield = ForceField(xmlFile)
        self.system = self.forcefield.createSystem(self.topology)
        self.force_Group2Name = {}

        # Construct indexing for DNA particles
        self.constructDNAIndex()

        self.verbose = verbose

    def constructDNAIndex(self):
        """Construct DNA atom and residue index, as well as their relations
        """
        self.dnaTypes = [_DNA_S1_TYPE, _DNA_S2_TYPE]

        # Atom index for each DNA strand
        self.dnaStrands = {}
        for dnaType in self.dnaTypes:
            self.dnaStrands[dnaType] = [atom.index for atom in self.atoms if atom.name == dnaType]
        
        # DNA residues (bairpairs)
        self.dnaRes = [residue.index for residue in self.residues if residue.name == _DNA_RES]
        self.nBp = len(self.dnaRes)

        # Res2Atom and Atom2Res relation for each type
        self.dna_Atom2Res = {}
        self.dna_Res2Atom = {}
        for dnaType in self.dnaTypes:
            self.dna_Atom2Res[dnaType] = {atom:self.atoms[atom].residue.index for atom in self.dnaStrands[dnaType]}
            self.dna_Res2Atom[dnaType] = {self.dna_Atom2Res[dnaType][atom]:atom for atom in self.dna_Atom2Res[dnaType]}
    
    def addForce(self, force, forceName):
        """Add a force to the system, set the force group and store the corresponding str name

        Args:
            forces (openmm force)
            forceName (str)
        """
        nForce = self.system.getNumForces()
        forceGroup = nForce
        force.setForceGroup(forceGroup)

        self.force_Group2Name[forceGroup] = forceName
        self.system.addForce(force)

    def addConnetivityForce(self, kCon=3000.0, bondLength=2.0, forceName = 'Con'):
        """Add a connectivity force to the system

        Args:
            kCon (float, optional): k in harmonic spring 0.5k(r-r0)^2. Defaults to 3000.0.
            bondLength (float, optional): r0 in harmonic spring 0.5k(r-r0)^2. Defaults to 2.0.
            forceName (str, optional): Defaults to 'Connectivity'.
        """
        if self.verbose:
            print("Building connectivity terms......", end = ' ')

        # set up the connectivity force between nearest neighbors on the same strand
        con = HarmonicBondForce()
        for dnaType in self.dnaTypes:
            for i in range(self.nBp - 1):
                con.addBond(self.dnaStrands[dnaType][i], self.dnaStrands[dnaType][i+1], bondLength, kCon)

        # add the connectivity force to our system
        self.addForce(con, forceName)
        if self.verbose:
            print("done")

    def addExcludVolumeForce(self, kExcl=5856.0, rExcl=2.07, forceName = 'Ex-vol'):
        """Add the connectivity term to the system"""
        if self.verbose:
            print("Building excluded volume force......", end = ' ')
        
        # set up the excluded volume force
        excl = CustomNonbondedForce("epsilon*step(sigma-r)*(r-sigma)^2; sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)")
        excl.addPerParticleParameter("sigma")
        excl.addPerParticleParameter("epsilon")

        # add all DNA atoms
        # note in naked DNA system, all atoms are DNA
        for atom in self.atoms:
            if atom.name in self.dnaTypes:
                excl.addParticle([rExcl, kExcl])
            else:
                raise Exception("None DNA atom detected in vanilla WeChroM system")

        # exclude nearest neighbors on the same strand
        # those connected by the connectivity term    
        for dnaType in self.dnaTypes:
            for i in range(self.nBp - 1):
                excl.addExclusion(self.dnaStrands[dnaType][i], self.dnaStrands[dnaType][i+1])    
        
        # Set cutoff distance and method
        excl.setCutoffDistance(rExcl)
        excl.setNonbondedMethod(excl.CutoffNonPeriodic)

        # Add to the system
        self.addForce(excl, forceName)
        if self.verbose:
            print("done")

    def addIntraStrandMemoryForce(self, temPos, depth = 0.3, width= 0.2, forceName = 'Intra_m'):
        """Add the intra-strand associative memory force to the system

        Args:
            temPos (dict): dict, 10 x 2 x memory_range x 3: positions of templates
            depth (float, optional): depth of the Gaussian Well. Defaults to 0.3. unit = kT
            width (float, optional): width of the Gaussian Well. Defaults to 0.2. unit = reduced length
            Available (depth, width) pairs: (0.15, 0.1),(0.2, 0.15), (0.25, 0.17), (0.3, 0.2)
        """
        # Apply the memory template repetitively along the DNA
        if self.verbose:
            print("Building DNA intra-strand associative memory force......", end = ' ')
        # Init the force in Gaussian Well
        intra = CustomBondForce(f'-depth_intra*exp((r-r0)^2/(-2.0*sigma^2))')
        intra.addGlobalParameter('depth_intra', depth)
        intra.addPerBondParameter('r0')
        intra.addPerBondParameter('sigma')

        # Initialize template distances temDist
        # temDist[dnaType][resDiff] should have size nTemplate * (memoryRange - resDiff)
        # in this way, it's equivalent to slide the memory along the DNA
        temDists = {}        
        for dnaType in self.dnaTypes:
            temDists[dnaType] = defaultdict(list)
        # fill in temDist
        for temId in temPos:
            for dnaType in temPos[temId]:
                positions = temPos[temId][dnaType]
                for i, pos1 in enumerate(positions):
                    for j in range(i + 1, len(positions)):
                        temDists[dnaType][j-i].append(np.linalg.norm(pos1 - positions[j]))
        
        # build associative memory bonds from templates
        for dnaType in self.dnaTypes:
            for i in range(self.nBp - 1):
                for j in range(i + 1, min(self.nBp, i + self.memoryRange)):
                    for temDist in temDists[dnaType][j-i]:
                        r0 = temDist
                        sigma = width * (j - i)**0.5
                        intra.addBond(self.dnaStrands[dnaType][i], self.dnaStrands[dnaType][j], [r0, sigma])
        
        # Sanity Check: number of bonds
        nBonds = intra.getNumBonds()
        nBondsExpect = num_intra_memory_bonds(self.nBp, self.memoryRange, len(temPos))
        if nBondsExpect != nBonds:
            print(f"\nALERT! Expect {nBondsExpect} bonds but got {nBonds} bonds for {forceName}")

        # Add to the system
        self.addForce(intra, forceName)

        if self.verbose:
            print("done")

    def addInterStrandMemoryForce(self, temPos, depth = 0.3, width= 0.2, forceName = 'Inter_m'):
        """Add the inter-strand associative memory force to the system

        Args:
            temPos (dict): dict, 10 x 2 x memory_range x 3: positions of templates
            depth (float, optional): depth of the Gaussian Well. Defaults to 0.3. unit = kT
            width (float, optional): width of the Gaussian Well, has a factor of 2 compared with the intra-strand width. Defaults to 0.2. unit = reduced length
            Available (depth, width) pairs: (0.15, 0.1),(0.2, 0.15), (0.25, 0.17), (0.3, 0.2) should keep the same as intra-strand
        """
        # Apply the memory template repetitively along the DNA
        if self.verbose:
            print("Building DNA inter-strand associative memory force......", end = ' ')
        # Init the force in Gaussian Well
        inter = CustomBondForce(f'-depth_inter*exp((r-r0)^2/(-2.0*sigma^2))')
        inter.addGlobalParameter('depth_inter', depth)
        inter.addPerBondParameter('r0')
        inter.addPerBondParameter('sigma')
        
        # inter-strand width factor
        width *= 2

        # Initialize template distances temDist
        # temDist[dnaType][resDiff] should have size nTemplate * (memoryRange - resDiff)
        # in this way, it's equivalent to slide the memory along the DNA
        temDists = {}        
        for dnaType in self.dnaTypes:
            temDists[dnaType] = defaultdict(list)

        # fill in temDist
        I, J = self.dnaTypes
        for temId in temPos:
            iPos = temPos[temId][I]
            jPos = temPos[temId][J]
            for id1 in range(len(iPos)):
                for id2 in range(id1, len(iPos)):
                    temDists[I][id2 - id1].append(np.linalg.norm(iPos[id1] - jPos[id2]))
                    temDists[J][id2 - id1].append(np.linalg.norm(jPos[id1] - iPos[id2]))
        
        # build associative memory bonds from templates
        for i in range(self.nBp - 1):
            for j in range(i + 1, min(self.nBp, i + self.memoryRange)):
                for temDist in temDists[I][j-i]:
                    r0 = temDist
                    sigma = width * (j - i)**0.5
                    inter.addBond(self.dnaStrands[I][i], self.dnaStrands[J][j], [r0, sigma])
                for temDist in temDists[J][j-i]:
                    r0 = temDist
                    sigma = width * (j - i)**0.5
                    inter.addBond(self.dnaStrands[J][i], self.dnaStrands[I][j], [r0, sigma])

        # fill in the same base-pair bonds
        for i in range(self.nBp):
            for temDist in temDists[I][0]:
                r0 = temDist
                sigma = width 
                inter.addBond(self.dnaStrands[I][i], self.dnaStrands[J][i], [r0, sigma])
        
        # Sanity Check: number of bonds
        nBonds = inter.getNumBonds()
        nBondsExpect = num_inter_memory_bonds(self.nBp, self.memoryRange, len(temPos))
        if nBondsExpect != nBonds:
            print(f"\nALERT! Expect {nBondsExpect} bonds but got {nBonds} bonds for {forceName}")

        # Add to the system
        self.addForce(inter, forceName)
        
        if self.verbose:
            print("done")

    def addDefaultForces(self):
        """Add default forces designed for the wechrom system, including the connectivity term, the excluded-volume term, generating associative memory templates and adding inter- and intra- strand associative memory terms
        """
        self.addConnetivityForce()
        self.addExcludVolumeForce()
        generate_naked_memory_files()
        self.temPos = get_naked_memory_template_bonds()
        self.addIntraStrandMemoryForce(self.temPos)
        self.addInterStrandMemoryForce(self.temPos)
            
    def initializeSimulation(self, platform = 'CPU', temperature = 300, collisionRate = 1.0, timeStep = 20.0):
        """Initialize an openmm simulation. 

        Args:
            platform (str, optional): platform. Defaults to 'CPU'. Can be 'CPU', 'CUDA' or 'OpenCL'
            temperature (int, optional): temperature in kelvin. Defaults to 300.
            collisionRate (float, optional): collision rate for the langevin integrator in 1/picosecond. Defaults to 1.0.
            timeStep (float, optional): simulation time step in femtosecond. Defaults to 20.0.

        Attributes:
            platform (openmm platform): CPU, CUDA or OpenCL. Defaults to CPU
            integrator (openmm integrator): a LangevinIntegrator with given temperature, collision rate and timestep
            simulation (openmm simulation): a simulation with wechrom topology, system, instance's integrator and given platform. The initial positions are set by wechrom pdbx file and the initial velocities are set to temperature. Local energy is minimized.
        """


        # set up integrator 
        self.platform = Platform.getPlatformByName(platform)
        self.temperature = temperature / 300 * 120.27  * kelvin # reduced temperature
        self.collisionRate = collisionRate * 2.75 / picosecond # reduced versed time
        self.timestep = timeStep * 0.001 / 2.75  * picoseconds # reduced time
        self.integrator = LangevinIntegrator(self.temperature, self.collisionRate , self.timestep)

        # set up simulation, initial positions and velocities
        self.simulation = Simulation(self.topology, self.system, self.integrator, self.platform)
        self.simulation.context.setPositions(self.positions)
        self.simulation.context.setVelocitiesToTemperature(self.temperature)
        self.simulation.minimizeEnergy() 
        if self.verbose:
            print("Langevin integrator and simulation initialized")

    def runSteps(self, steps = 1000, reportFreq = 100, append = False, outputDir = 'DEFAULT', energyFilename = 'energy.txt', dcdFilename = 'movie.dcd', chkFilename ='checkpt.chk' ):
        """Run simulation on openmm. Should be called after initializeSimulation(). Report energy for each force and trajectories in dcd format at the given report frequency. Save the last checkpoint

        Args:        
            outputDir (str, optional): output directory. Defaults to 'DEFAULT'.
            steps (int, optional): total steps for the simulation. Defaults to 1000.
            reportFreq (int, optional): report frequency. Defaults to 100.
            append (bool, optional): append to an old trajectory and energy or not. Defaults to False.
            energyFilename (str, optional): file name to report energy. Defaults to 'energy.txt'.
            dcdFilename (str, optional): file name to report trajectory. Defaults to 'movie.dcd'.
            chkFilename (str, optional): file name for the last checkpoint. Defaults to 'checkpt.chk'.
        """
        
        if outputDir == 'DEFAULT':
            outputDir = os.getcwd()
        self.outputDir = outputDir

        # convert steps and reportFreq to int
        steps = int(steps)
        reportFreq = int(reportFreq)
        assert steps % reportFreq == 0, "Steps must be devisible by reportFreq"

        if self.verbose:
            print(f"Simulation will take {steps} steps and get reported every {reportFreq} steps")
        # Set up reporters
        # a DCD reporter to store the trajectory 
        # dcd_file = os.path.join(self.outputDir, dcdFilename)
        # self.simulation.reporters.append(DCDReporter(dcd_file, reportFreq, append = append))

        # a checkpoint reporter to store the last checkpoint
        chk_file = os.path.join(self.outputDir, chkFilename)
        self.simulation.reporters.append(CheckpointReporter(chk_file, reportFreq))

        # initialize the energy reporter
        if not append:
            with open(energyFilename, 'w') as fe:
                header ='Steps ' + ' '.join(['{0:<8s}'.format(i) for i in self.force_Group2Name.values()]) + ' Total Kinetics(kT)\n'
                fe.write(header)
        
        if self.verbose:
            print("----------------Simulation Starts----------------")

        # add a dcd writer manually to get it properly closed at the end of simulation
        dcdFilepath = os.path.join(self.outputDir, dcdFilename)
        dcdFilehandle = open(dcdFilepath, 'bw')
        dcdWriter = DCDFile(dcdFilehandle, self.topology, self.timestep)

        # run the simulation
        nSavedFrames = steps // reportFreq
        for i in tqdm(range(nSavedFrames)):
            self.simulation.step(reportFreq)

            # potential energy from each group
            groupEnergy = []
            for group in self.force_Group2Name:
                state = self.simulation.context.getState(getEnergy=True, groups={group})
                groupEnergy.append(state.getPotentialEnergy().value_in_unit(kilojoule_per_mole))
            # total potential energy
            state = self.simulation.context.getState(getEnergy=True, getPositions=True)
            groupEnergy.append(state.getPotentialEnergy().value_in_unit(kilojoule_per_mole))
            # kinetic energy
            groupEnergy.append(state.getKineticEnergy().value_in_unit(kilojoule_per_mole))
            # positions
            positions = state.getPositions(True)
            dcdWriter.writeModel(positions)
            
            with open(energyFilename, 'a') as fe:
                curStep = (i+1)*reportFreq
                line = f'{curStep:<8}' + ' '.join(['{0:<8.2f}'.format(_) for _ in groupEnergy]) + '\n'
                fe.write(line)
        if self.verbose:
            print("\nSimulation done.")
            print(f"Please check your trajectory file {dcdFilename}, energy file {energyFilename}", end = ' ')
            print(f"at your output directory {self.outputDir}")
        # close the dcd file 
        dcdFilehandle.close()
