import os
import pkg_resources
import pandas as pd
from MDAnalysis import Universe
from simtk.openmm.app import PDBxFile

_DNA_TYPES = ['I', 'J']

def read_memory_template(template_list_dir = 'DEFAULT', template_list_filename = 'data/naked_database_list'):
    """Read the list of associative memory template filenames

    Args:
        template_list_dir (str, optional): directory of the template list file. Defaults to 'DEFAULT'.
        template_list (str, optional): name of the template list file. Defaults to 'data/naked_database_list'.

    Raises:
        Exception: custom template list is not supported right now

    Returns:
        pandas.DataFrame: fist column 'filename' and second column 'n_bp'
    """
    if template_list_dir == 'DEFAULT':
        template_location = pkg_resources.resource_filename(__name__, template_list_filename)
    else:
        raise Exception("Custom memory template not implemented yet")
    memory_list = pd.read_csv(template_location, sep = '\s+')
    return memory_list

def generate_memory_files(memory_dir = 'DEFAULT', memory_filename = 'naked.mem', memory_range = 4):
    """generate the file to store associative memory fragments

    Args:
        memory_dir (str, optional): Defaults to current workin directory.
        memory_filename (str, optional): Defaults to 'naked.mem'.
        memory_range (int, optional): Defaults to 4.
    """
    print("Preparing the associative memory files......", end = ' ')
    if memory_dir == 'DEFAULT':
        memory_dir = os.getcwd()
    # Save the memory fragments in desired file
    memory_location = os.path.join(memory_dir, memory_filename)
    with open(memory_location, 'w') as mo:
        # TO-DO: update query with system name
        mo.write('[Target]\nquery\n[Memories]\n')
        mo.write('filename template_res memory_range\n')

    # Read the template dataframe
    memory_list = read_memory_template()

    # Generate memory fragments from the memory templates
    with open(memory_location, 'a') as mo:
        for _, row in memory_list.iterrows():
            fragment = (row['n_bp'])//2
            if fragment + memory_range < row['n_bp']:
                mo.write(f"{row['filename']} {fragment} {memory_range}\n")

    print("done")

def get_memory_template_bonds(memory_list_dir = 'DEFAULT', memory_list_filename = 'naked.mem', memory_location = 'data/memory_template/'):
    """Extract the positions of template associative memory fragments

    Args:
        memory_list_dir (str, optional): directory of the memory list. Defaults to current workin directory.
        memory_list_filename (str, optional): filename of the memory list. Defaults to 'naked.mem'.
        memory_location (str, optional): location of the database. Defaults to 'data/memory_template/'.

    Returns:
        dict, 10 x 2 x memory_range x 3: positions of templates
            layer 1(10): template, keys 0-9
            layer 2(2): dna type (strand), keys 'I', 'J'
            layer 3(np.array, memory_range x 3): positions of template residues in the range of memory_range
    """
    if memory_list_dir == 'DEFAULT':
        memory_list_dir = os.getcwd()
    # read the associative memory list files
    memory_list_location = os.path.join(memory_list_dir, memory_list_filename)
    memory_list = pd.read_csv(memory_list_location, skiprows=4, sep='\s+', names=['filename', 'template_res', 'memory_range'])
    
    template_positions = {}
    # read each associative memory fragment
    for index, row in memory_list.iterrows():
        fragment_filename = pkg_resources.resource_filename(__name__, memory_location + row['filename'])

        # read the memory file
        template_positions[index] = {}
        frag = Universe(PDBxFile(fragment_filename))

        # Read the bps within the memmory_length for each strand
        for dna_type in _DNA_TYPES:
            template_query = f"resid {row['template_res']}:{row['template_res'] + row['memory_range'] - 1} and name {dna_type}"     
            template_atoms = frag.select_atoms(template_query)
            # get the positions
            template_positions[index][dna_type] = template_atoms.positions / 10 # MDanalysis default units Angstrom, convert to openmm default units nm

    return template_positions