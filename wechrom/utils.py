import os
from MDAnalysis import Universe
import numpy as np
from simtk.openmm.app import pdbxfile, topology, element, pdbfile

def coarse_grain_atompdb_2_wechromcif(pdb_filename, out_dir = os.getcwd(), out_file_prefix = 'coarse_grain', out_pdb = False):
    """Coarse grain an atomistic DNA pdb file into a wechrom cif file.
    ***YOU MAY NEED TO MODIFY THIS FUNCTION DEPENDING ON YOUR PDB FILE***

    Args:
        pdb_filename (str or file): pdb file
        out_dir (str, optional): output directory. Defaults to os.getcwd().
        out_file_prefix (str, optional): output file prefix. Defaults to 'coarse_grain'.
        out_pdb (bool, optional): whether to generate a corase grained pdb file. Defaults to False.
    """
    u = Universe(pdb_filename)

    # read DNA and corase grain 
    # ***You may need to modify this section depending on your pdb file***
    nucleic_atoms = u.select_atoms("nucleic")
    I = nucleic_atoms.segments[0].atoms
    J = nucleic_atoms.segments[1].atoms
    I_pos = np.array([res.atoms.center_of_mass() for res in I.residues]) 
    J_pos = np.array([res.atoms.center_of_mass() for res in reversed(J.residues)]) # reverse order
    
    # convert to reduced units
    I_pos /= 0.241
    J_pos /= 0.241

    # sanity check
    n_bp = I_pos.shape[0]
    assert n_bp == J_pos.shape[0], "I and J strands should have same number of beads"  

    # build wechrom topology 
    topo = topology.Topology()
    ele_DNA = element.Element.getBySymbol("N")

    # add DNA atoms to chain 0
    # each DNA residue has name "BP"
    # each DNA atom has name "I"/"J", element "N"
    chain = topo.addChain(0) 
    for i in range(n_bp):
        residue = topo.addResidue(name = "BP", chain = chain) 
        atom = topo.addAtom(name = "I", element=ele_DNA, residue=residue)
        atom = topo.addAtom(name = "J", element=ele_DNA, residue=residue)

    positions = [None] * (2*n_bp)
    positions[::2] = I_pos
    positions[1::2] = J_pos
    
    # set unit cell dimensions
    dim = np.absolute(np.amax(positions, axis = 0) - np.amin(positions, axis=0))
    topo.setUnitCellDimensions(np.round(dim*2))

    # write corase grained files
    with open(os.join(out_dir, out_file_prefix+".cif"), "w") as fo:
        pdbxfile.PDBxFile.writeFile(topo, positions, fo)
    if out_pdb:
        with open(os.join(out_dir, out_file_prefix+".pdb"), "w") as fo:
            pdbfile.PDBFile.writeFile(topo, positions, fo)

def num_intra_memory_bonds(n_bp, memory_range, n_template):
    """Calculate the expected number of bonds for intra-strand associative memory force

    Args:
        n_bp (int): number of basepairs in the system
        memory_range (int): size of the memory
        n_template (int): number of templates used as memories

    Returns:
        int: number of expected intra-strand bonds
    """
    n_bonds = 0
    for i in range(1, memory_range):
        n_bonds += (n_bp - i) * (memory_range - i)
    n_bonds *= 2 * n_template
    return n_bonds
    
def num_inter_memory_bonds(n_bp, memory_range, n_template):
    """Calculate the expected number of bonds for inter-strand associative memory force

    Args:
        n_bp (int): number of basepairs in the system
        memory_range (int): size of the memory
        n_template (int): number of templates used as memories

    Returns:
        int: number of expected inter-strand bonds
    """
    n_bonds = 0
    for i in range(1, memory_range):
        n_bonds += (n_bp - i) * (memory_range - i)
    n_bonds *= 2 * n_template
    n_bonds += n_bp * memory_range * n_template
    return n_bonds