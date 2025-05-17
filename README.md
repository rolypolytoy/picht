# Picht
Electron optics in electric and magnetic lenses. Uses the axisymmetric view, computes static E- and B-fields, enables the parameterization of ion and electron beams, and allows you to view their trajectories in the vicinity of various kinds of electromagnetic lenses, with physically accurate and relativistically corrected dynamics. Allows for custom mesh sizes using the finite difference method, and is incredibly performant due to integrations with PyAMG, Numba, and Joblib. Works on all operating systems, is unit-tested, and works locally and on Jupyter Notebook with all dependencies handled by PyPi.

## Installation
```bash
pip install picht
```
[![PyPI version](https://img.shields.io/pypi/v/picht.svg)](https://pypi.org/project/picht/) ![tests](https://github.com/rolypolytoy/picht/actions/workflows/tests.yml/badge.svg) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15376240.svg)](https://doi.org/10.5281/zenodo.15376240)

## Gallery

API documentation, computational physics methods, and tutorials can be found at the official website: https://rolypolytoy.github.io/picht/. 

All of the examples in the gallery were scripted with under 100 lines of code and generated in under a minute of real time.

### Electric Lens
Focusing electrons with a cylindrical electrostatic lens. Reference implementation [here](https://rolypolytoy.github.io/picht/auto_examples/01_example_cylindrical_lens.html)
![cylindricallens](https://github.com/user-attachments/assets/3e3e00a7-009d-4cf6-a026-860ee545307e)

### Magnetic Lens
Focusing electrons with a magnetic lens. Reference implementation [here](https://rolypolytoy.github.io/picht/auto_examples/05_example_magnetic.html)
![ezgif-61b10d1bda4d4b](https://github.com/user-attachments/assets/333a8a55-f8b3-47ad-8dd7-268c9f8edb46)

### Einzel Lens
Focusing electrons with three electrodes in a unipotential lens arrangement. Reference implementation [here](https://rolypolytoy.github.io/picht/auto_examples/02_example_einzel_focusing.html).
![einzel](https://github.com/user-attachments/assets/44ba059c-c857-42fe-8d60-5b0ea997467f)

### Scanning Electron Microscope
Controlling electrons with a Wehnelt cylinder, cathode/anode, a condenser (einzel) lens and an objective (einzel) lens. Reference implementation [here](https://rolypolytoy.github.io/picht/auto_examples/04_example_sem.html).

![sem5](https://github.com/user-attachments/assets/4d6d943f-faa6-4b7e-8e78-7df707f885b3)
