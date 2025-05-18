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

`Picht` is architecturally focused on parameterization rather than scripting. The classes are constrained, and creating a beam is simple, however adjusting parameters massively affects trajectories and outcomes. The ideal workflow is creating a system, observing the user interface, and tweaking voltage, magnetomotive force, permeability, or geometric features, until the system fits its intended purpose. It thus enables research or education in electron optics for amateurs and experts alike, due to its intentionally constrained syntax, but computationally efficient physics engine.

# Computational Physics

The main equation `Picht` uses to solve for its magnetic field is Poisson's equation:

$$\nabla^2 A = -\mu_0 \mu_r J$$

Where $\nabla^2$ is the Laplacian operator, A is the vector potential, $\mu_0$ is the permeability of free space, $\mu_r$ is the dimensionless value for relative permeability of the material, and J is the magnetomotive force in ampere-turns [jackson1999classical]. We use Dirichlet masking to specify magnetomotive force values, as well as relative permeability values, where the user specifies a magnetic lens is. Then, PyAMG solves for vector potential based on the resolution of the grid the user specified.

From this, we can use the equation: 

$$B = \nabla \times A$$

To compute the magnetic field 'B' from the curl of the vector potential A. We use the axisymmetric formulation for curl, by making the assumption that the system is cylindrically symmetric around its optical axis, which is a reasonable assumption to make in electron optics [hawkes1972electron].

To compute the electric field, we solve Laplace's equation first:

$$\nabla^2 V = 0$$

Where V is electric potential  [jackson1999classical]. We use Dirichlet masking to set the exact voltage amount the user specified at the points they specified, and then compute solutions. From the electric potential field, the electric field can be calculated through the following equation: 

$$E = -\nabla V$$

Given the electric and magnetic field components at every point, particle trajectories can be computed. For electric effects, we use the Lorentz equation for electrostatics:

$$F = qE$$

Where q is charge and E is the electric field component. However, for the magnetic field, we don't use the Lorentz equation, because the full Lorentz treatment magnetic focusing is prohibitively expensive to simulate. To provide accuracy and efficiency, we instead use the paraxial ray equation for magnetic lenses, which is a reasonable approximation at small angles from the axis [hawkes1972electron]. The form we use in our code is:

$$F_r =  -\frac{q^2 B_z^2}{4m} r$$

Where $F_r$ is the radial component of force, q is the charge of the particle, $B_z$ is the axial component of the magnetic field, m is the mass of the particle, and r is the distance of the particle from the central axis.

We use second-order (5-stencil) finite difference methods for electric field calculations, and fourth-order (9-stencil) finite difference methods for magnetic field calculations, to provide greater accuracy for magnetic fields to offset some deviations from reality due to using the paraxial ray equation rather than full Lorentz treatment. Combining this treatment of electric and magnetic components of force enables us to simulate systems with both kinds of lenses in an accurate and performant manner.

# References
