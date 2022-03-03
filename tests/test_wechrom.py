from simtk.openmm import HarmonicBondForce, CustomNonbondedForce
import os
from wechrom import WechromSystem


# test input file
_dir_path = os.path.dirname(os.path.realpath(__file__))
_cif_path = os.path.join(_dir_path, 'naked_2bp.cif')
_test_bp = 2

def test_pass():
    assert True, "dummy sample test"

def test_WechromSystem_init():
	"""A comprehensive test for WechromSystem.__init__()
	"""

	# test dna_res setup
	we = WechromSystem(_cif_path)
	assert we.nBp == _test_bp, "Check num of base pairs"

	# test dna_atom2res setup
	for dnaType in we.dnaTypes:
			assert len(we.dna_Atom2Res[dnaType]) == _test_bp, "Check dna atom 2 res relation"

	# test openmm system setup
	assert we.system.getNumParticles() == 2 * _test_bp, "Check num of system particle"

def test_WechromSystem_addConnetivityForce():
	"""A comprehensive test for WechromSystem.addConnetivityForce()
	"""
	we = WechromSystem(_cif_path)
	we.addConnetivityForce()

	force = we.system.getForce(1)
	assert isinstance(force, HarmonicBondForce), "Force should be HarmonicBondForce"

	nBonds = force.getNumBonds()
	assert nBonds == 2*(_test_bp - 1), "Check num of Harmonic Bonds"

def test_WechromSystem_addExcludVolumeForce():
	"""A comprehensive test for WechromSystem.addExcludVolumeForce()
	"""
	we = WechromSystem(_cif_path)
	we.addExcludVolumeForce()

	force = we.system.getForce(1)
	assert isinstance(force, CustomNonbondedForce), "Force should be CustomNonbondedForce"

	nParticles = force.getNumParticles()
	assert nParticles == 2*_test_bp, "Check num of particles"

	nExcl = force.getNumExclusions()
	assert nExcl == 2*(_test_bp - 1), "Check num of exclustion"

def test_WechromSystem_initializeSimulation():
	we = WechromSystem(_cif_path)
	we.initializeSimulation()

	