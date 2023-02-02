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
    affiliation: 2, 1 # (Multiple affiliations must be quoted)
  - name: Maude Charmetant # note this makes a footnote saying 'co-first author'
    orcid : 0000-0002-3137-1328 
    affiliation: 2
affiliations:
 - name: Deutsches Zentrum für Luft- und Raumfahrt e.V. (DLR) Projektträger, Joseph-Beuys-Allee 4, 53113 Bonn, Germany 
   index: 1
 - name: Argelander-Institut für Astronomie, University of Bonn, Auf dem Hügel 71, 53121 Bonn, Germany
   index: 2

date: 1 Mach 2022
bibliography: paper.bib

# Statement of need

The study of the sky at microwave wavelengths is rich in information for cosmology 
and astrophysics. In particular, the study of the Cosmic Microwave Background (CMB) has led to tremendous advances in our understanding of the Universe, its age, its history, its composition, and the properties of its components. The CMB is a relic radiation that carries information from the hot and dense Universe roughly 380.000 years after the big bang. The study of emissions from distant galaxies through the Cosmic Infrared Background allows us to probe the Large-Scale Structure (LSS) of the cosmos. Emissions, such as the Thermal and Kinematic Sunyaev-Zeldovich effect, respectively allow us to detect galaxy clusters and measure the LSS velocity field. Besides these extragalactic sources, our own Milky Way Galaxy emits microwaves through various mechanisms like synchrotron and bremsstrahlung emission or thermal radiation from interstellar dust grains, all of which allow to probe complex physics, such as star formation or astrochemistry. A major challenge in studies of the microwave sky is the separation of these different astrophysical components and the reduction of instrumental effects like noise. Sophisticated component separation and data analysis software are required and usually custom-made for each experiment. To develop and test such software and create forecasts for future instruments, realistic mock data sets are needed. With the recent push towards ground-based wide-field multifrequency surveys, a tool to produce mock images of the microwave based on the most recent observations and high-resolution simulations is needed. We provide an easy-to-use and open-source tool to create high-resolution maps of the microwave sky in the frequency range between $27$ GHz and $860$ GHz.


# Summary
%must be for non-specialists and less than 1000 words long. 

`skymodel` is a Python package that generates high-resolution ($p_{size}\approx 0.86'$) 
full-sky maps in the HEALPix [@Healpy_2005] format that include various components of 
the microwave sky at frequencies $\nu \in [27,860]$ GHz. It makes use of extra-galactic 
template maps from three existing simulations: WebSky [@CITA_2020], Sehgal [@Sehgal_2010] 
and SO [@SO_2019]. Galactic foregrounds are simulated based on maps generated with the 
Python Sky Model (PySM) [@Ben_2017]. To function it requires having HEALPy, PySM, numpy, 
astropy[@astropy], and full-sky extra-galactic template maps that can be downloaded [here](https://uni-bonn.sciebo.de/s/zgPsb7qvXTnNsrO/authenticate) 
and are derived by spectral modeling and subsequent processing of the WebSky\footnote{\url{https://mocks.cita.utoronto.ca/data/websky/v0.0/}}, 
SO and Sehgal extragalactic maps \footnote{\url{https://lambda.gsfc.nasa.gov/simulation/tb_sim_ov.cfm}}.

`skymodel` is designed around one easy-to-use function that generates
a full-sky map of the microwave sky with the desired parameters so that it can 
be easily and intuitively used by researchers in Astrophysics. The package was initially designed for the (Cerro Chajnantor Atacama Telescope - prime) 
CCAT-prime collaboration, whose upcoming Fred Young Submillimeter Telescope (FYST) first 
light is scheduled for 2023. Having a code that allows modeling the sky as FYST would 
see it is essential to make predictions. The skymodel was already used in the first collaboration 
paper to make predictions [@CCAT_2021]. 
Even though PySM3[@PySM_2021] extends the PySM code, not only to simulate galactic foregrounds but also extra-galactic foregrounds, both in intensity and polarisation, the skymodel aim is different. The skymodel offers more flexibility in terms of the base simulations to generate the 
extra-galactic components, some like the Sehgal, going to very high-resolution ($p_{size}\approx 0.43'$). Choosing components, switching between simulations, unit conversion between $\mathrm{K}_{\mathrm{CMB}}$, $\mathrm{K}_{\mathrm{RJ}}$, and MJy/sr$^{-1}$, choice of spatial resolution are made easy thanks to a slim user interface through a single function. The skymodel deals with all the calculation and ensure consistency between the different data for the user. The skymodel can simulate any instrumental white noise and atmospheric red noise tailored for the upcoming FYST telescope. 


# Acknowledgements

The authors would like to thank Ben Thorne, Kaustuv Basu, Frank Bertoldi, Steve Choi, and the CCAT-p collaboration members for insightful discussions. MC and JE acknowledge support from Frank Bertoldi, by the Bonn-Cologne Graduate School of Physics and Astronomy (BCGS) and the International Max Planck Research School (IMPRS). The galactic foregrounds are generated using PySM [@Ben_2017]. The simulations of the extra-galactic components used in this paper were developed by the WebSky Extragalactic CMB Mocks team, with the continuous support of the Canadian Institute for Theoretical Astrophysics (CITA), the Canadian Institute for Advanced Research (CIFAR), and the Natural Sciences and Engineering Council of Canada (NSERC), and were generated on the GPC supercomputer at the SciNet HPC Consortium. SciNet is funded by the Canada Foundation for Innovation under the auspices of Compute Canada, the Government of Ontario, the Ontario Research Fund – Research Excellence, and the University of Toronto. This work made use of Astropy:\footnote{http://www.astropy.org} a community-developed core Python package and an ecosystem of tools and resources for astronomy.


# References
