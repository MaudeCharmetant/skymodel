.. CCATp sky model documentation master file, created by
   sphinx-quickstart on Thu Jun  3 13:50:11 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the skymodel's documentation!
===========================================

skymodel provides a model of the microwave sky based on templates maps from numerical simulations of extragalactic sources with the most recent Galactic foregrounds maps from PySM. The sky model generate full-sky maps of the microwave sky containing the components of you choice at the frequencies of your choice between ~[27GHz,853GHz].


Citing
------

For now, no paper is attached to the code, mentioning the Github page is sufficient. 

Installation
------------

.. toctree::
   :maxdepth: 2
   :caption: Contents:

skymodel is a pure Python module. You'll need is numpy, astropy, PySM , tqdm, and healpy . The extra-galactic simulations templates are not provided but can be downloaded from WebSky , SO&Sehgal then simply change the data path at the beginning of the code file to indicate where you saved the data.
From source

ccatp_sky_model is developed on GitHub and can be installed by cloning the source repository and install from there
 
1. git clone https://github.com/MaudeCharmetant/skymodel.git 
2. cd skymodel 
3. python setup.py install

   
 
Reference
---------

.. toctree::
   :maxdepth: 2
   
   code


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
