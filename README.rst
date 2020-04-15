ccatp_sky_model
===============

**ccatp_sky_model** 
provides a model of the microwave sky that will be used for 
science forecasts for the CCAT-prime survey. The sky model combines maps from
numerical simulations of extragalactic sources with the most recent Galactic
foregrounds maps based on Planck data.


Install the CCATp sky model
===========================

ccatp_sky_model is a pure Python module and should therefore be pretty easy to install.
All you'll need is numpy, astropy, PySM <https://github.com/bthorne93/PySM_public> , tqdm, and healpy <https://github.com/healpy/healpy>.

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


Copyright 2020 Maude Charmetant, Jens Erler, and contributors.

The ccatp_sky_model is free software made available under the MIT License. For details see
the LICENSE file.
