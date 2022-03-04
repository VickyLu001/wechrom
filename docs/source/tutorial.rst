Tutorial
=========

Usage
---------
We provide a description for general usage here. It might be helpful to learn the usage by the two examples in the followed section :ref:`Example 1: Naked DNA`.

After installing the wechrom package following the :ref:`Installation` section, you should be able to import the modules by

.. code-block:: python

    import wechrom

We provide two classes, ``WechromSystem`` for DNA-only systems and ``SingleNucleoSystem`` for single nucleosomes. To initiate the system, you need to provide a pdbx/mmcif file or a pdb file to set up the topology.

.. code-block:: python

    your_system = WechromSystem(your_cif_file) # vanilla DNA system

or

.. code-block:: python

    your_system = SingleNucleoSystem(your_cif_file) # single nucleosome system

The default forces designed in WEChroM can be applied by

.. code-block:: python

    your_system.addDefaultForces()

You can also add any custom forces supported by openmm

.. code-block:: python

    your_system.addForce(your_custom_force, "name_of_your_force")

At this stage, you can make use of ``your_system.topology``, ``your_system.system`` to set up your openmm simultions, or you can also initialize the simulation with default langevin integrator and wechrom topology with our method

.. code-block:: python

    your_system.initializeSimulation(platform = 'CPU') # platform can be 'CPU', 'CUDA' or 'OpenCL' depending on your openmm installation

After simulation initialization, you can make use of ``your_system.simulation``, ``your_system.integrator`` to run simulations, or you can use our default dcd reporter and energy reporter to run the simulations

.. code-block:: python

    your_system.runSteps(steps = your_steps, reportFreq = your_reportFreq)
    

Example 1: Naked DNA
---------------------
This section will instruct you to build a vanilla wechrom system (DNA only) and run simulations. You can find the example notebook at ``wechrom/examples/naked_75bp/naked_75bp.ipynb``. You can run this example interactively using the Jupyter notebook.

After installing the wechrom package following the :ref:`Installation` section, you should be able to import the modules, classes and functions from the package. In this example, we utilize the ``WechromSystem`` class.

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
    Please check your trajectory file movie.dcd, energy file energy.txt at your output directory your_path\examples\naked_75bp

    
Example 2: Single nucleosome
---------------------------------
This section will instruct you to build a single nucleosome wechrom system and run simulations. You can find the example notebook at ``wechrom/examples/nucleosome_223bp/nucleosome_223bp.ipynb``. You can run this example interactively using the Jupyter notebook.

After installing the wechrom package following the :ref:`Installation` section, you should be able to import the modules, classes and functions from the package. In this example, we utilize the ``SingleNucleoSystem`` class.

.. code-block:: python

    import wechrom

We have prepared a pdbx to set up the system.Initialize the nucleosome system with a pdbx file we prepared with 147 bp wrapped DNA, 38 bp linker DNA on each end and a histone core particle. This cif file also includes two virtual sites at the two ends for external force illustration. 

.. code-block:: python

    nuc_223bp = wechrom.SingleNucleoSystem("singleN_L38_endvs.cif", verbose = True)

Apply default forces designed for wechrom

.. code-block:: python

    nuc_223bp.addDefaultForces()

Output with verbose:

.. code-block::

    Building connectivity terms...... done
    Building excluded volume force...... done
    Preparing the associative memory files...... done
    Building DNA intra-strand associative memory force...... done
    Building DNA inter-strand associative memory force...... done
    Preparing the nucleosome associative memory files...... done
    Building nucleosome center associative memory force...... done
    Building nucleosome neighbor associative memory force...... done

Here we provide an example of applying external force to the system. The desired external force is to pull the two ends apart in x direction. You may use any force supported by openmm

.. code-block:: python

    from openmm import CustomExternalForce
    def stretch_term(we, k_stretching=2.0, forceDirect="x", appliedSite=0):
        stretching = CustomExternalForce(f"({k_stretching})*({forceDirect})")
        stretching.addParticle(we.virtualSites[appliedSite])
        return stretching

    nuc_223bp.addForce(stretch_term(nuc_223bp, forceDirect="x", appliedSite=0), "pull_h")
    nuc_223bp.addForce(stretch_term(nuc_223bp, forceDirect="-x", appliedSite=-1), "pull_t")

Initialize the simulation with default integrator and wechrom topology

.. code-block:: python

    nuc_223bp.initializeSimulation(platform='CUDA')

Output with verbose:

.. code-block::

    Langevin integrator and simulation initialized

Here we provide an example of adding a pdb reporter to the simulation
.. code-block:: python

    from openmm.app import PDBReporter # openmm version >= 7.6
    # from simtk.openmm.app import PDBReporter # openmm version < 7.6
    steps = 1e3
    reportFreq = 1e1
    nuc_223bp.simulation.reporters.append(PDBReporter("./movie.pdb", reportFreq))

Run simulation with trajectory and energy reported.

.. code-block:: python

    nuc_223bp.runSteps(steps = steps, reportFreq = reportFreq)

Output with verbose:

.. code-block::

    Simulation will take 1000 steps and get reported every 100 steps
    ----------------Simulation Starts----------------
    100%|██████████| 10/10 [00:02<00:00, 38.37it/s]
    Simulation done.
    Please check your trajectory file movie.dcd, energy file energy.txt at your output directory your_path\examples\nucleosome_223bp