import os
import pkg_resources
import pandas as pd
import numpy as np
from collections import defaultdict
from MDAnalysis import Universe
try:
    # for openmm version >= 7.6
    from openmm.app import PDBxFile
except:
    # for openmm version <= 7.5
    from simtk.openmm.app import PDBxFile
from .utils import _DNA_S1_TYPE, _DNA_S2_TYPE, _PRO_TYPE

_DNA_TYPES = [_DNA_S1_TYPE, _DNA_S2_TYPE]

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

def generate_naked_memory_files(memory_dir = 'DEFAULT', memory_filename = 'naked.mem', memory_range = 4):
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
        mo.write('[Target]\nnaked DNA\n[Memories]\n')
        mo.write('filename template_res memory_range\n')

    # Read the template dataframe
    memory_list = read_memory_template(template_list_filename = 'data/naked_database_list')

    # Generate memory fragments from the memory templates
    with open(memory_location, 'a') as mo:
        for _, row in memory_list.iterrows():
            fragment = (row['n_bp'])//2
            if fragment + memory_range < row['n_bp']:
                mo.write(f"{row['filename']} {fragment} {memory_range}\n")

    print("done")

def get_naked_memory_template_bonds(memory_list_dir = 'DEFAULT', memory_list_filename = 'naked.mem', database_location = 'data/memory_template/'):
    """Extract the positions of template associative memory fragments

    Args:
        memory_list_dir (str, optional): directory of the memory list. Defaults to current workin directory.
        memory_list_filename (str, optional): filename of the memory list. Defaults to 'naked.mem'.
        database_location (str, optional): location of the database. Defaults to 'data/memory_template/'.

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
        fragment_filename = pkg_resources.resource_filename(__name__, database_location + row['filename'])

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

def generate_nucleo_memory_files(memory_dir = 'DEFAULT', memory_filename = 'nucleo.mem', database_location = 'data/memory_template/', neighbor_radius = 100.0):
    """generate the file to store associative memory fragments for nucleosomes

    Args:
        memory_dir (str, optional): Defaults to current workin directory.
        memory_filename (str, optional): Defaults to 'nucleo.mem'.
        database_location (str, optional): Defaults to 'data/memory_template/'.
        neighbor_radius (float, optional): the radius to select neighbor contacts. in Angstrom for MDAnalysis. Defaults to 100.0
    """
    print("Preparing the nucleosome associative memory files......", end = ' ')
    if memory_dir == 'DEFAULT':
        memory_dir = os.getcwd()
    
    # Read the template dataframe
    memory_list = read_memory_template(template_list_filename = 'data/nucleo_database_list')

    # Save the memory fragments in desired file
    memory_location = os.path.join(memory_dir, memory_filename)
    with open(memory_location, 'w') as mo:
        mo.write('[Target]\nnucleosome\n[Memories]\n')
        mo.write('filename nucleo_center type1 res1 type2 res2_list\n')

    file_query = defaultdict(dict)
    file_center = {}
    # Generate memory fragments from the memory templates
    with open(memory_location, 'a') as mo:
        for _, row in memory_list.iterrows():
            filename = row['filename'] 
            nucleo_center = row['nucleo_center']
            contacts = [int(x) for x in row['contacts'][1:-1].split(",")]
            dna_type = row['type']

            # generate center memory (between hitone and contact DNA)
            type1 = _PRO_TYPE
            res1 = row['histone_res']
            type2 = dna_type
            res2_list = contacts            
            mo.write(f"{filename} {nucleo_center} {type1} {res1} {type2} [" + ",".join([str(int(x)) for x in res2_list]) + "]\n")

            # prepare for neighbor memory
            file_center[filename] = nucleo_center

            # prepare query to select contact atoms from the file
            file_query[filename][dna_type] = f"name {dna_type} and ("
            for res_id in contacts:
                file_query[filename][dna_type] += f"resid {res_id + 1} or " #openmm index starts from 0, mda from 1
            file_query[filename][dna_type] = file_query[filename][dna_type][:-4] + ")"
    
    for filename in file_query:
        # go through memory files one by one
        fragment_filename = pkg_resources.resource_filename(__name__, database_location + filename)
        frag = Universe(PDBxFile(fragment_filename))

        # select atoms for the two strands
        dna_atoms = {}
        for dna_type in _DNA_TYPES:
            dna_atoms[dna_type] = frag.select_atoms(file_query[filename][dna_type])

        for res1_type in _DNA_TYPES:
            for res2_type in _DNA_TYPES:
                res1_res2_set = defaultdict(set)
                res1_atoms = dna_atoms[res1_type]
                res2_atoms= dna_atoms[res2_type]
                for atom1 in res1_atoms:
                    for atom2 in res2_atoms:
                        res1_id = atom1.resid-1 #openmm index from 0, mdanalysis index from 1
                        res2_id = atom2.resid-1
                        if res2_id < res1_id:
                            continue
                        if res2_id == res1_id and res1_type == res2_type:
                            continue
                        pos1 = atom1.position
                        pos2 = atom2.position
                        rm = np.linalg.norm(pos1 - pos2)
                
                        if rm < neighbor_radius:
                            res1_res2_set[res1_id].add(res2_id)
                for res1_id in res1_res2_set:
                    res1_res2_set[res1_id] = sorted(res1_res2_set[res1_id])
                    res2_list = "["+",".join([str(int(x)) for x in res1_res2_set[res1_id]])+"]"
                    with open(memory_location, 'a') as mo:
                        mo.write(f"{filename} {file_center[filename]} {res1_type} {res1_id} {res2_type} {res2_list}\n")
    print("done")

def get_nucleo_memory_template_bonds(memory_list_dir = 'DEFAULT', memory_list_filename = 'nucleo.mem', database_location = 'data/memory_template/'):
    """Extract the positions of template associative memory fragments

    Args:
        memory_list_dir (str, optional): directory of the memory list. Defaults to current workin directory.
        memory_list_filename (str, optional): filename of the memory list. Defaults to 'nucleo.mem'.
        database_location (str, optional): location of the database. Defaults to 'data/memory_template/'.

    Returns:
        dict, 2 x N: lengths of center templates bonds, 1st layer is type, 2nd layer is relative res_id
        dict, 2 x 2 x N1 x N2: lengths of neighbor templates bonds between the pair of (type, res_id) atoms
    """
    if memory_list_dir == 'DEFAULT':
        memory_list_dir = os.getcwd()
    # read the associative memory list files
    memory_list_location = os.path.join(memory_list_dir, memory_list_filename)
    memory_list = pd.read_csv(memory_list_location, skiprows=4, sep='\s+', names=['filename', 'nucleo_center', 'type1', 'res1', 'type2', 'res2_list'])
    
    # bonds dict initialization
    centerBonds = {}
    for dna_type in _DNA_TYPES:
        centerBonds[dna_type] = {}
    neighborBonds = {}
    for dna_type1 in _DNA_TYPES:
        neighborBonds[dna_type1] = {}
        for dna_type2 in _DNA_TYPES:
            neighborBonds[dna_type1][dna_type2] = defaultdict(dict)

    # read each associative memory fragment
    old_filename = None
    for index, row in memory_list.iterrows():
        if old_filename != row['filename']:
            fragment_filename = pkg_resources.resource_filename(__name__, database_location + row['filename'])
            # read the memory file
            frag = Universe(PDBxFile(fragment_filename))
        old_filename = row['filename']

        nucleo_zero = row['nucleo_center']
        type1 = row['type1']
        res1 = row['res1']
        type2 = row['type2']
        res2_list = [int(x) for x in row['res2_list'][1:-1].split(",")]

        

        atom1 = frag.select_atoms(f"name {type1} and resid {res1 + 1}")[0] #openmm index starts from 0, mda from 1
        pos1 = atom1.position / 10 # MDanalysis default units Angstrom, openmm default units nm
        for res2 in res2_list:
            atom2 = frag.select_atoms(f"name {type2} and resid {res2 + 1}")[0] #openmm index starts from 0, mda from 1
            pos2 = atom2.position / 10 # MDanalysis default units Angstrom, openmm default units nm

            rm = np.linalg.norm(pos1 - pos2)
            if type1 == _PRO_TYPE:
                centerBonds[type2][res2 - nucleo_zero] = rm
            else:
                neighborBonds[type1][type2][res1 - nucleo_zero][res2 - nucleo_zero] = rm

    return centerBonds, neighborBonds
