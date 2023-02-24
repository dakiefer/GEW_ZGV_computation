# GEW ZGV computation [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7537441.svg)](https://doi.org/10.5281/zenodo.7537441)

**Compute zero-group-velocity (ZGV) points of guided elastic waves (GEWs).** 

Three different computational techniques to locate ZGV points on dispersion curves are implemented. They are all based on the discretized waveguide problem.

- **Newton-type iteration:** super fast but needs initial guesses.
- **Method of fixed relative distance (MFRD):** scans a wavenumber interval without initial guesses and is likely to locate all ZGV points but is substantially slower. It refines computed approximations with the Newton-type iteration.
- **Direct method:** does not need initial guesses and guarantees to find all ZGV points. It is slow and can, therefore, only be used with rather small matrices.

<img src="https://user-images.githubusercontent.com/3725269/210227577-cfc9a367-2dd5-4c1c-95b0-f1c1a3cbfdce.png"  alt="ZGV points on dispersion curves" width="400px" />

Code repository: [<img src="https://www.svgrepo.com/show/35001/github.svg" alt="GitHub" width="27px" />](https://github.com/dakiefer/gew_zgv_computation) [https://github.com/dakiefer/gew_zgv_computation](https://github.com/dakiefer/gew_zgv_computation)

The methods have been presented in:

> D. A. Kiefer, B. Plestenjak, H. Gravenkamp, and C. Prada, “Computing zero-group-velocity points in anisotropic elastic waveguides: Globally and locally convergent methods,” The Journal of the Acoustical Society of America, vol. 153, no. 2, pp. 1386–1398, Feb. 2023, doi: [10.1121/10.0017252](http://doi.org/10.1121/10.0017252)

## How to use

1. Change into the `GEW_ZGV_computation` folder or add it to the Matlab path.
2. Execute `example.m` . Enjoy!

## Dependencies

The *direct method* is based on the solver for singular two-parameter eigenvalue problems implemented by Bor Plestenjak and Andrej Muhič in `MultiParEig`: 
> Bor Plestenjak (2023). MultiParEig (https://www.mathworks.com/matlabcentral/fileexchange/47844-multipareig), MATLAB Central File Exchange. Retrieved January 14, 2023.

## Authors

Code created 2022–2023 by

Bor Plestenjak, Faculty of Mathematics and Physics, University of Ljubljana, Slovenia<br/>
[bor.plestenjak@fmf.uni-lj.si](bor.plestenjak@fmf.uni-lj.si) &nbsp; ● &nbsp; Follow me on [ResearchGate](https://www.researchgate.net/profile/Bor-Plestenjak)!

Daniel A. Kiefer, Institut Langevin, ESPCI Paris | PSL, France<br/>
[daniel.kiefer@espci.fr](mailto:daniel.kiefer@espci.fr) &nbsp; ● &nbsp; [dakiefer.net](https://dakiefer.net) &nbsp; ● &nbsp; Follow me on [ResearchGate](https://www.researchgate.net/profile/Daniel-Kiefer-5)!

If this code is useful to you, please cite it as:

> B. Plestenjak and D. A. Kiefer, GEW ZGV computation [Computer software], 2023. doi: [10.5281/zenodo.7537441](http://doi.org/10.5281/zenodo.7537441)

and also the related publication:

> D. A. Kiefer, B. Plestenjak, H. Gravenkamp, and C. Prada, “Computing zero-group-velocity points in anisotropic elastic waveguides: Globally and locally convergent methods,” The Journal of the Acoustical Society of America, vol. 153, no. 2, pp. 1386–1398, Feb. 2023, doi: [10.1121/10.0017252](http://doi.org/10.1121/10.0017252)

[<img src="https://user-images.githubusercontent.com/3725269/210226492-a2a56855-f1ea-4e36-96bf-0c66e28e6d7b.svg" alt="Logo FMF" width="190px" />](https://www.fmf.uni-lj.si/en/) &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [<img src="https://user-images.githubusercontent.com/3725269/210226450-69833a58-8e9e-4dca-8401-2ab5119cc1da.svg" alt="Logo Institut Langevin" width="210px" />](https://www.institut-langevin.espci.fr) &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [<img src="https://user-images.githubusercontent.com/3725269/210226449-caa43ed7-0385-4546-a9ad-5d3c99339b8c.svg" alt="Logo ESPCI Paris" width="190px" />](https://www.espci.psl.eu/en/)
