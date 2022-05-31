---
title: 'skymodel: Python package for versatile high-resolution microwave sky simulation'
tags:
  - Python
  - astronomy
  - CMB
  - galactic science
  - extra-galactic science
  - microwave sky
authors:
  - name: Jens Erler
    orcid : 0000-0002-5486-6745
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Maude Charmetant # note this makes a footnote saying 'co-first author'
    orcid : 0000-0002-3137-1328 
    affiliation: 2
  - name: Frank Bertoldi 
    affiliation: 2
affiliations:
 - name: Deutsches Zentrum für Luft- und Raumfahrt (DLR)
   index: 1
 - name: Argelander-Institut für Astronomie, University of Bonn, Auf dem Hügel 71, 53121 Bonn, Germany
   index: 2

date: 1 Mach 2022
bibliography: paper.bib

# Summary
%must be for non specialists and less than 1000 words long. 

The study of microwave sky emissions is rich in information for Cosmology 
and Astrophysics. Observation of extra-galactic components such as the Cosmic Microwave Background (CMB) allowed 
to put unpreceded constraints on cosmological parameters and confirmed the big bang 
theory [@COBE_1994][@Bennett_2003][@Planck_2014]. The Cosmic Infrared Background (CIB) is a probe of the Large Scale Structure (LSS) [@Lenz_2019]
of the universe. The thermal Sunyaev-Zeldovich (tSZ) effect probe the integrated pressure, 
the hot baryons in the Universe and was vastly used to detect galaxy clusters [@Planck_2011]. 
Another flavor of it, the kinematic SZ (kSZ) effect probe the bulk motion of those 
clusters and the LSS [@Hern_2015]. All the galactic foregrounds, galactic dust, synchrotron, free-free emission, 
and Anomalous Microwave Emission (AME) are often treated as contaminants but can also 
be studied to learn more about our galaxy properties. Therefore, having a code 
that allows to easily generate those components, together, or separately at any frequency, 
resolution or unit is essential to many different fields. 


# Statement of need

`skymodel` is a Python package that generate high-resolution ($p_{size}\approx 0.86'$) full-sky map in HEALPy [@Healpy_2005]
format for the various components of the Microwave sky at the desired frequency
$\nu \in [27,860]$GHz. It makes use of extra-galactic template maps from three existing 
simulations WebSky [@CITA_2020], Sehgal [@Sehgal_2010] and SO [@SO_2019], and of the PySM [@Ben_2017] code to simulate galactic foregrounds.
Switching components on and off, switching between simulations, units conversion
$\mathrm{K}_{\mathrm{CMB}}$, $\mathrm{K}_{\mathrm{RJ}}$, MJy/sr$^{-1}$, beam, instrumental white noise and atmospheric red noise handling 
is made easy thanks to one single function. The skymodel deals with all the calculation 
and ensure consistency between the different data for the user. To function it requires 
having HEALPy, PySM, numpy, astropy[@astropy], and the extra-galactic templates maps that are publicly 
available and can be downloaded from WebSky\footnote{\url{https://mocks.cita.utoronto.ca/data/websky/v0.0/}}
SO&Sehgal\footnote{\url{https://lambda.gsfc.nasa.gov/simulation/tb_sim_ov.cfm}}.

`skymodel` is designed so that only one short easy-to-use function could generate
a full-sky map of the microwave sky with the desired parameters so that it can 
be used by researchers in Astrophysics that are not experts on microwave sky and 
by students. It was also designed for the (Cerro Chajnantor Atacama Telescope - prime) 
CCAT-prime collaboration whose upcoming Fred Young Submillimeter Telescope (FYST) first 
light is scheduled for 2023. Having a code that allows modeling the sky as FYST would 
see it is essential to make predictions. The skymodel was already used in the first collaboration 
paper to make predictions [@CCAT_2021]. Note, that a more advanced version of the PySM software, 
PySM3[@PySM_2021] extends the PySM code not only to simulate galactic foregrounds but also extra-galactic foregrounds
thus simulation the microwave sky exactly like the skymodel. However, differences exists, PySM3 
does not only simulate Intensity maps of the microwave sky but also Polarization, which the skymodel 
does not do. The skymodel offers more variety in terms of the base simulations to generate the 
extra-galactic components, some like the Sehgal, going to very high-resolution ($p_{size}\approx 0.43'$).

# Acknowledgements

What to put here ? 

# References