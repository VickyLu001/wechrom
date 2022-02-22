Tutorial
=========

Naked DNA
---------
This section will instruct you to build a vanilla wechrom system (DNA only) and run simulations. You can find the example notebook at ``wechrom/examples/naked_75bp/naked_75bp.ipynb``. You can run this example interactively
using the Jupyter notebook.

After installing the wechrom package following the :ref:`Installation` section, you should be able to import the modules, classes and functions from the package. In this example, we utilize the `WechromSystem` class.

.. code-block:: python

    from wechrom import WechromSystem

We have prepared a pdbx to set up the system. It's a  corase grained 75-bp naked DNA molecule. You can coarse grain an atomistic DNA pdb file into a wechrom cif file with our utility function.

.. code-block:: python

    from wechrom import coarse_grain_atompdb_2_wechromcif

    # a cif file will be generated at your working directory
    coarse_grain_atompdb_2_wechromcif(YOUR_PDB_FILE, out_dir = os.getcwd(), out_file_prefix = 'coarse_grain')

Initialize the wechrom system with a pdbx/mmcif file

.. code-block:: python

    naked_75bp = WechromSystem("naked_75bp.cif", verbose = True)

Apply default forces designed for wechrom

.. code-block:: python

    naked_75bp.addDefaultForces()

Output with verbose:

.. code-block::

    Building connectivity terms...... done
    Building excluded volume force...... done
    Preparing the associative memory files...... done
    Building intra-strand associative memory force...... done
    Building inter-strand associative memory force...... done

Initialize the simulation with default integrator and wechrom topology

.. code-block:: python

    naked_75bp.initializeSimulation()

Output with verbose:

.. code-block::

    Langevin integrator and simulation initialized

Run simulation with trajectory and energy reported.

.. code-block:: python

    naked_75bp.runSteps(steps = 1000, reportFreq = 100)

Output with verbose:

.. code-block::

    Simulation will take 1000 steps and get reported every 100 steps
    ----------------Simulation Starts----------------
    100%|██████████| 10/10 [00:02<00:00,  4.95it/s]
    Simulation done.
    Please check your trajectory file movie.dcd, energy file energy.txt at your output directory d:\wechrom\examples\naked_75bp