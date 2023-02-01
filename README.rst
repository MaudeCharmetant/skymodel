ccatp_sky_model
===============

**ccatp_sky_model** 
is a tool that generates high-resolution maps of the microwave sky based on numerical 
simulations of extragalactic emission and the python sky model (PySM) library, which generates 
maps of galactic foreground emission. The sky model allows generating the most relevant microwave
sky components and return them in any desired combination. Maps can be returned at arbitrary 
frequencies between 30 and 860 GHz, a variety of commonly used units, custom angular resulution, 
and including instrumental and/or atmospheric noise components. It was build to be easy to use 
and versatile for anyone working on the microwave sky. 

.. image:: https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white
.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org/
.. image:: https://img.shields.io/badge/license-MIT-red.svg?style=flat
    :target: https://github.com/MaudeCharmetant/CCATp_sky_model/blob/master/LICENSE    
.. image:: https://img.shields.io/badge/Doc-readthedoc-green.svg
    :target: https://skymodel.readthedocs.io/en/latest/index.html 
.. image:: https://img.shields.io/badge/Doc-PDF-green.svg
    :target: https://github.com/MaudeCharmetant/CCATp_sky_model/blob/master/Documentation.pdf


How to install the CCATp sky model
==================================

ccatp_sky_model is a pure Python module.
You'll need the `NumPy <https://numpy.org/>`_ , `Astropy <https://www.astropy.org/>`_, `PySM <https://github.com/bthorne93/PySM_public>`_ , `tqdm <https://github.com/tqdm/tqdm>`_, and `healpy <https://github.com/healpy/healpy>`_ packages. The required template maps can be dowloaded `here <https://uni-bonn.sciebo.de/s/zgPsb7qvXTnNsrO>`_ with the password: pw4referee, then simply change the data path at the beginning of the code file to indicate where you saved the data. The original extra-galactic simulation templates can be downloaded from `WebSky <https://mocks.cita.utoronto.ca/data/websky/v0.0/>`_ and `SO&Sehgal <https://lambda.gsfc.nasa.gov/simulation/tb_sim_ov.cfm>`_ .

.. _source:

From source
-----------

ccatp_sky_model is developed on `GitHub <https://github.com/MaudeCharmetant/CCATp_sky_model>`_ and can be 
installed by cloning the source repository and install from there:

.. code-block:: bash

    git clone https://github.com/MaudeCharmetant/CCATp_sky_model.git
    cd CCATp_sky_model
    python setup.py install


Test the installation
=====================

To make sure that the installation was successful and to familiarise yourself with 
ccatp_sky_model, we recommend running the provided jupyter notebooks that can be found in
the /examples directory. 


Read the Doc page
==================
A webpage containing all the functions, inputs, outputs, and a short description of their 
roles, can be found at:  https://skymodel.readthedocs.io/en/latest/index.html 


Community guidelines
====================

Contributions are welcome and appreciated please submit an issue or a pull request. 


Copyright 2020 Maude Charmetant, Jens Erler, and contributors.

The ccatp_sky_model is free software made available under the MIT License. For details see
the LICENSE file.
