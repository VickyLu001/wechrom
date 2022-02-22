from wechrom import prepare_memory
import tempfile
import os

def test_pass():
    assert True, "dummy sample test"

def test_read_memory_template():
    # test if we get a non-empty list
    memory_list = prepare_memory.read_memory_template()
    assert 'filename' in memory_list.columns
    assert len(memory_list) > 0

def test_generate_memory_files():
    # test if the memory list file is generated
    outfile_handle, outfile_path = tempfile.mkstemp()
    os.close(outfile_handle)
    dir, filename = os.path.split(outfile_path)
    try:
        prepare_memory.generate_memory_files(memory_dir=dir, memory_filename=filename)
        with open(outfile_path, 'r') as fi:
            contents = fi.read()
    finally:
        os.remove(outfile_path)
    assert contents.startswith('[Target]\nquery\n[Memories]\n')

def test_get_memory_template_bonds():
	"""An output size test for WEChroMSystem.getMemoryTemplateBonds()
	This function depends on prepare_memory.generate_memory_files
	"""
	outfile_handle, outfile_path = tempfile.mkstemp()
	os.close(outfile_handle)
	dir, filename = os.path.split(outfile_path)
	try:
		prepare_memory.generate_memory_files(memory_dir=dir, memory_filename=filename)
		template_positions = prepare_memory.get_memory_template_bonds(memory_list_dir = dir, memory_list_filename = filename)
		n_template = len(template_positions)
		n_type = len(template_positions[0])
		assert n_type == 2, "two dna types"
		for dna_type in template_positions[0]:
			memory_length, n_dim = template_positions[0][dna_type].shape
			assert n_dim == 3, "x, y, z in positions"
			break
	finally:
		os.remove(outfile_path)