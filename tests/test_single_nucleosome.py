from wechrom import SingleNucleoSystem
from simtk.openmm import CustomNonbondedForce
import os

# test input file
_dir_path = os.path.dirname(os.path.realpath(__file__))
_cif_path = os.path.join(_dir_path, 'nucleo_2bp.cif')
_test_bp = 2
_test_pro_res = 1

def test_SingleNucleoSystem_init():
	"""A simple test for SingleNucleoSystem.__init__()
	"""
	# test dnaRes setup
	we = SingleNucleoSystem(_cif_path)
	assert we.nBp == _test_bp, "Check num of base pairs"

	# test proRes setup
	assert len(we.proRes) == _test_pro_res, "Check num of protein residues"

    # test openmm system setup
	assert we.system.getNumParticles() == 2 * _test_bp + _test_pro_res, "Check num of system particle"

def test_addExcludVolumeForce():
    """A simple test for SingleNucleoSystem.addExcludVolumeForce()
	"""

    we = SingleNucleoSystem(_cif_path)
    we.addExcludVolumeForce()

    force = we.system.getForce(1)
    assert isinstance(force, CustomNonbondedForce), "Force should be CustomNonbondedForce"

    nParticles = force.getNumParticles()
    assert nParticles == 2*_test_bp + _test_pro_res, "Check num of particles"

    nExcl = force.getNumExclusions()
    assert nExcl == 2*(_test_bp - 1), "Check num of exclustion"