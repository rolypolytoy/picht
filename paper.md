---
title: 'Picht: Electron Optics in Python using the Finite Difference Method'
tags:
  - Python
  - electron optics
  - electrodynamics
  - relativistic physics
authors:
  - name: Rishiit Sharma
    orcid: 0009-0006-4895-2686
    equal-contrib: true
    affiliation: 1
affiliations:
  - index: 1
    name: None
date: 17 May 2025
bibliography: paper.bib
---

# Summary

Electron optics is the field of physics that deals with the study of electron and ion beams under the influence of electric and magnetic 'lenses'. Electromagnetic lenses are called 'lenses' because they make charged particle beams converge, diverge and reflect, in a way analogous to how optical lenses affect beams of light.  Electron optics is an important field of modern science, and is an integral part of the operations of many machines including electron microscopes, ion implantation and milling machines, and mass spectrometers. Picht provides researchers and engineers with the ability to design and simulate such systems, to reduce the barrier of entry to perform meaningful research in this field. 

# Statement of Need

`Picht` is a Python package to simulate electron optics systems. Using Python enables an easy and intuitive API, without sacrificing on performance when compared to writing the package in a compiled language (like C++). This is because it uses libraries like Numba [lam2015numba], which provide performance in some cases comparable to CUDA while maintaining the standards of usability expected of libraries in the Python ecosystem  [@9092407]. `Picht` exists to provide an open-source alternative to commercial electrodynamics tools like SIMION and COMSOL, without sacrifices in user-friendliness or performance. It solves magnetic fields using an optimized multigrid solver by PyAMG [pyamg2023], and uses joblib [joblib2020] for parallelizing trajectory calculations. 

`Picht` is architecturally focused on parameterization rather than scripting. The classes are constrained, and creating a beam is simple, however adjusting parameters massively affects trajectories and outcomes. The intended workflow is creating a system, observing the user interface, and tweaking voltage, magnetomotive force, permeability, or geometric features, until the system fits its intended purpose. It thus enables research or education in electron optics for amateurs and experts alike, due to its intentionally constrained syntax, but computationally efficient physics engine.

# Computational Physics



# References
