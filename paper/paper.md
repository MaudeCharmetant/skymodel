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

The study of the microwave sky emissions is rich of informations for Cosmology 
and Astrophysics. Obervation of extra-galactic components such as the Cosmic Microwave Background (CMB) allowed 
to put unpreceed constrains on cosmological parameters and confirmed the big bang 
theory. The Cosmic Infrared Background (CIB) is a probe of the large sclae structure (LSS)
of the universe. The thermal Sunyaev-Zeldovich (tSZ) effect probe the integrated pressure, 
the hot baryons in the Universe, this was wastly used to detect galaxy clusters. 
Another flavour of it, the kinematic SZ (kSZ) effect probe the bulk motion of those 
clusters and of the LSS. All the galactic foregrounds, galactid dust, syncrothron, free-free emission
and Anomalous Microwave Emission (AME) are often treated as contaminants but can also 
be studied to learn more about our galaxy properties. Therefore, having a code 
that allow to easily generate those components, together, or separately at any frequency, 
resolution or units is essential to many differents fields. 


# Statement of need

`skymodel` is a Python package that generate high-resolution ($p_{size}\approx 0.86'$) full-sky map in HEALPy [@healpy]
format for the various components of the Microwave sky at a desired frequency
$\nu \in [27,860]$GHz. It make use of extra-galactic templates maps from three existing 
simulations WebSky, Sehgal and SO and of the PySM code to simulate galactic foregrounds.
Switching components on and off, switching between simulations, units conversion
$\mathrm{K}_{\mathrm{CMB}}$, $\mathrm{K}_{\mathrm{RJ}}$, MJy/sr$^{-1}$, beam, instrumental white noise and atmospheric red noise handeling 
is made easy thanks to one single function. The skymodel deals with all the calculation 
and ensure consistency between the different datas for the user. To function it requires 
have HEALPy, PySM, numpy, astropy and the extra-galactic templates maps that are publicly 
available and can be downloaded from WebSky\footnote{\url{https://mocks.cita.utoronto.ca/data/websky/v0.0/}}
SO&Sehgal\footnote{\url{https://lambda.gsfc.nasa.gov/simulation/tb_sim_ov.cfm}}.

`skymodel` is designed so that only this one function could be run by anyone wishing to 
simulate de Microwave sky. 

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
