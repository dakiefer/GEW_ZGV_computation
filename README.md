# GEW ZGV computation

**Compute zero-group-velocity (ZGV) points of guided elastic waves (GEWs).** 

Three different computational techniques to locate ZGV points are implemented. They are all based on the discretized guided wave problem.

- **Newton-type iteration:** super fast but needs initial guesses.
- **Method of fixed relative distance (MFRD):** scans a wavenumber interval without initial guesses and is likely to locate all ZGV points but is substantially slower. It refines computed approximations with the Newton-type iteration.
- **Direct method:** does not need initial guesses and guarantees to find all ZGV points. It is slow and can, therefore, only be used with rather small matrices.

![dispersion_fig](https://user-images.githubusercontent.com/3725269/185589376-c991b579-2550-40e8-8d7b-d90de9edec77.png)

Code repository: [<img src="https://www.svgrepo.com/show/35001/github.svg" alt="GitHub" width="27px" />](https://github.com/dakiefer/gew_zgv_computation) [https://github.com/dakiefer/gew_zgv_computation](https://github.com/dakiefer/gew_zgv_computation)

The methods have been presented in:

> D. A. Kiefer, B. Plestenjak, H. Gravenkamp, and C. Prada, "Computing zero-group-velocity points in anisotropic elastic waveguides: globally and locally convergent methods." arXiv, Nov. 2022. doi: [10.48550/arXiv.2211.01995](http://doi.org/10.48550/arXiv.2211.01995)

## How to use

1. Change into the `GEW_ZGV_computation` folder or add it to the Matlab path.
2. Execute `example.m` . Enjoy!

## Dependencies

The *direct method* is based on the solver for singular two-parameter eigenvalue problems implemented by Bor Plestenjak and Andrej Muhič in `MultiParEig`: 
> Bor Plestenjak (2023). MultiParEig (https://www.mathworks.com/matlabcentral/fileexchange/47844-multipareig), MATLAB Central File Exchange. Retrieved January 2, 2023.

## Authors

Code created 2022–2023 by

Bor Plestenjak, Faculty of Mathematics and Physics, University of Ljubljana, Slovenia<br/>
[bor.plestenjak@fmf.uni-lj.si](bor.plestenjak@fmf.uni-lj.si) &nbsp; ● &nbsp; Follow me on [ResearchGate](https://www.researchgate.net/profile/Bor-Plestenjak)!

Daniel A. Kiefer, Institut Langevin, ESPCI Paris | PSL, France<br/>
[daniel.kiefer@espci.fr](mailto:daniel.kiefer@espci.fr) &nbsp; ● &nbsp; [dakiefer.net](https://dakiefer.net) &nbsp; ● &nbsp; Follow me on [ResearchGate](https://www.researchgate.net/profile/Daniel-Kiefer-5)!

If this code is useful to you, please cite it as:

> TODO

and also the related publication:

> D. A. Kiefer, B. Plestenjak, H. Gravenkamp, and C. Prada, "Computing zero-group-velocity points in anisotropic elastic waveguides: globally and locally convergent methods." arXiv, Nov. 2022. doi: [10.48550/arXiv.2211.01995](http://doi.org/10.48550/arXiv.2211.01995)

[<img src="https://user-images.githubusercontent.com/3725269/185571121-f5fcd518-32de-40b2-b4b1-f4ef0610ccd1.svg" alt="Logo Institut Langevin" width="250px" />](https://www.institut-langevin.espci.fr) &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [<img src="https://user-images.githubusercontent.com/3725269/185570398-ca2796ab-2bd3-4171-a7a6-af1f74014504.svg" alt="Logo ESPCI Paris" width="260px" />](https://www.espci.psl.eu/en/)
