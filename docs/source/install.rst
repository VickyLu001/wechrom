Installation
============

Requirements
------------

wechrom requires the following python libraries:

* **openmm**
* **numpy**
* **pandas**
* **MDAnalysis**
* **setuptools**
* **tqdm**
  
To install **openmm**, you may refer to the `openmm installation guide <http://docs.openmm.org/7.5.0/userguide/application.html#installing-openmm>`_.

Install with pip
-------------------

.. code-block:: bash

    $ pip install wechrom

Install from source
-------------------
Download the source code from

.. code-block:: bash
    
    $ git clone https://github.com/VickyLu001/wechrom.git

Then you can install the package by

.. code-block:: bash

    $ python ./wechrom/setup.py install

or install it in develop mode

.. code-block:: bash

    $ pip install -e ./wechrom