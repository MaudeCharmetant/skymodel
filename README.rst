ccatp_sky_model
===============

**ccatp_sky_model** 
generate high-resolution maps of the microwave sky using maps 
from numerical simulations of extragalactic emissions and PySM 
code to generate maps of the galactic foregrounds. 
The sky model allows generating any component or combination of 
components at a desired frequency, unit, beam size, including 
instrumental and atmospheric noises. 
It was thought to be easy and versatile to use for anyone working
on the microwave sky. 

.. image:: https://img.shields.io/badge/Doc-readthedoc-green.svg
    :target: https://skymodel.readthedocs.io/en/latest/index.html 
.. image:: https://img.shields.io/badge/Doc-PDF-green.svg
    :target: https://github.com/MaudeCharmetant/CCATp_sky_model/blob/master/Documentation.pdf
    


Install the CCATp sky model
===========================

ccatp_sky_model is a pure Python module.
You'll need is numpy, astropy, `PySM <https://github.com/bthorne93/PySM_public>`_ , tqdm, and `healpy <https://github.com/healpy/healpy>`_ .
The extra-galactic simulations templates are not provided but can be downloaded from `WebSky <https://mocks.cita.utoronto.ca/data/websky/v0.0/>`_ , `SO&Sehgal <https://lambda.gsfc.nasa.gov/simulation/tb_sim_ov.cfm>`_ then simply change the data path at the beginning of the code file to indicate where you saved the data.

.. _source:

From source
-----------

ccatp_sky_model is developed on `GitHub <https://github.com/MaudeCharmetant/CCATp_sky_model>`_ and can be 
installed by cloning the source repository and install from there

.. code-block:: bash

    git clone https://github.com/MaudeCharmetant/CCATp_sky_model.git
    cd CCATp_sky_model
    python setup.py install


Test the installation
=====================

To make sure that the installation went alright and to familiarise yourself with 
ccatp_sky_model, we recommend running the provided jupyter notebooks that can be found in
the /examples directory. 


Read the Doc page
==================
A webpage containing all the functions, inputs, outputs, and a short description of their 
roles, can be found at:  https://skymodel.readthedocs.io/en/latest/index.html 


Community guidelines
====================

Contributions are welcome and appreciated please submit an Issue or a Pull request. 


Copyright 2020 Maude Charmetant, Jens Erler, and contributors.

The ccatp_sky_model is free software made available under the MIT License. For details see
the LICENSE file.
