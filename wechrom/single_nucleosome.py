import pandas as pd
from .wechrom import *
from .prepare_memory import generate_nucleo_memory_files, get_nucleo_memory_template_bonds
from.utils import _PRO_RES, _PRO_TYPE, _VIRTUAL_TYPE

_NUCLEO_XML = pkg_resources.resource_filename(__name__, 'data/nucleo.xml')

class SingleNucleoSystem(WechromSystem):
    """Build a single-nucleosome WEChroM system with topology and system information for openmm simulations

    Args:
        WechromSystem (object): DNA only system

    Attributes:
        proTypes (list of str): names of protein atoms
        protein (list of int): index of protein atoms
        proRes (list of int): index of protein residues
        pro_Atom2Res (dict): atom id to res id relation for proteins
        pro_Res2Atom (dict): res id to atom id relation for proteins
        nucZero (int): res id of the zero point DNA bp of the nucleosome
    """
    def __init__(self, pdbxFile, xmlFile='DEFAULT', memoryRange=4, nucZero = -1, verbose=False):
        if xmlFile == 'DEFAULT':
            xmlFile = _NUCLEO_XML

        # inherent the parent initialization
        super().__init__(pdbxFile, xmlFile, memoryRange, verbose)

        # atom index for protein
        self.proTypes = [_PRO_TYPE]
        self.protein = [atom.index for atom in self.atoms if atom.name == _PRO_TYPE]
        if verbose:
            print(f"{len(self.protein)} histone particle(s) detected")

        # residues for protein
        self.proRes = [residue.index for residue in self.residues if residue.name == _PRO_RES]

        self.pro_Atom2Res = {atom:self.atoms[atom].residue.index for atom in self.protein}
        self.pro_Res2Atom = {self.pro_Atom2Res[atom]:atom for atom in self.pro_Atom2Res}

        # if zero point of the nucleosome DNA is not specified, use the central bp of the DNA molecule
        if nucZero == -1:
            self.nucZero = self.nBp//2
        else:
            self.nucZero = nucZero

        self.virtualSites = [atom.index for atom in self.atoms if atom.name == _VIRTUAL_TYPE]
        # sanity check
        for virtualSite in self.virtualSites:
            assert self.system.isVirtualSite(virtualSite), f"Virtual Sites {virtualSite} not aligned correctly, please check the index in your pdbx file"


    def addExcludVolumeForce(self, kExcl=5856.0, rExclDna=2.07, rExclPro = 26.0, forceName = 'Ex-vol'):
        """Add the connectivity term to the system.
        Overide the original method in WechromSystem"""
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
                excl.addParticle([rExclDna, kExcl])
            elif atom.name in self.proTypes:
                excl.addParticle([rExclPro, kExcl])
            elif atom.name == _VIRTUAL_TYPE:
                excl.addParticle([0.0, 0.0])
            else:
                raise Exception(f"Atom name {atom.name} not recognized")

        # exclude nearest neighbors on the same strand
        # those connected by the connectivity term    
        for dnaType in self.dnaTypes:
            for i in range(self.nBp - 1):
                excl.addExclusion(self.dnaStrands[dnaType][i], self.dnaStrands[dnaType][i+1])    
        
        # Set cutoff distance and method
        excl.setCutoffDistance((rExclDna + rExclPro)/2)
        excl.setNonbondedMethod(excl.CutoffNonPeriodic)

        # Add to the system
        self.addForce(excl, forceName)
        if self.verbose:
            print("done")

    def addNucleoCenterMemoryForce(self, centerBonds, centerDepth = 0.4, centerWidth = 3.0, forceName = "cen_nuc"):
        """Add the nucleosome center associative memory term to the system.
        """
        if self.verbose:
            print("Building nucleosome center associative memory force......", end = ' ')
        
        # initiliaze center and neighbor forces
        center = CustomBondForce('-depth_center*exp((r-r0)^2/(-2.0*sigma_center^2))')
        center.addGlobalParameter('depth_center', centerDepth)
        center.addGlobalParameter('sigma_center', centerWidth)
        center.addPerBondParameter('r0')

        proteinAtom = self.protein[0]
        for dnaType in centerBonds:
            for dnaRes in centerBonds[dnaType]:
                targetRes = dnaRes + self.nucZero
                dnaAtom = self.dna_Res2Atom[dnaType][targetRes]
                center.addBond(proteinAtom, dnaAtom, [centerBonds[dnaType][dnaRes]])

        self.addForce(center, forceName)
        if self.verbose:
            print("done")

    def addNucleoNeighborMemoryForce(self, neighborBonds, neighborDepth = 0.2, neighborWidth = 1.0, forceName = "ngb_nuc"):
        """Add the nucleosome neighbor associative memory term to the system.
        """
        if self.verbose:
            print("Building nucleosome neighbor associative memory force......", end = ' ')
        neighbor = CustomBondForce('-depth_neighbor*exp((r-r0)^2/(-2.0*sigma_neighbor^2))')
        neighbor.addGlobalParameter('depth_neighbor', neighborDepth)
        neighbor.addGlobalParameter('sigma_neighbor', neighborWidth)
        neighbor.addPerBondParameter('r0')

        for type1 in neighborBonds:
            for type2 in neighborBonds[type1]:
                for res1 in neighborBonds[type1][type2]:
                    for res2 in neighborBonds[type1][type2][res1]:
                        targetRes1 = res1 + self.nucZero
                        targetRes2 = res2 + self.nucZero
                        dnaAtom1 = self.dna_Res2Atom[type1][targetRes1]
                        dnaAtom2 = self.dna_Res2Atom[type2][targetRes2]
                        neighbor.addBond(dnaAtom1, dnaAtom2, [neighborBonds[type1][type2][res1][res2]])
        self.addForce(neighbor, forceName)
        if self.verbose:
            print("done")
        
    def addDefaultForces(self):
        """Add default forces designed for the wechrom system, including all DNA forces, updated excluded-volume forces and nucleosome associative memory forces
        """
        super().addDefaultForces()

        generate_nucleo_memory_files()
        centerBonds, neighborBonds = get_nucleo_memory_template_bonds()
        self.addNucleoCenterMemoryForce(centerBonds)
        self.addNucleoNeighborMemoryForce(neighborBonds)
