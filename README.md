# ConservationLawsParticles.jl

<!-- ![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg) -->
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://FedericoStra.github.io/ConservationLawsParticles.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://FedericoStra.github.io/ConservationLawsParticles.jl/dev)
[![Build Status](https://github.com/FedericoStra/ConservationLawsParticles.jl/workflows/CI/badge.svg)](https://github.com/FedericoStra/ConservationLawsParticles.jl/actions)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

Particle methods for 1D conservation laws.

This package implements the deterministic particle schemes described in the article

> **E. Radici, F. Stra**:
> *Entropy solutions of mildly singular non-local scalar conservation laws with congestion via deterministic particle method*.
> SIAM Journal on Mathematical Analysis 55.3 (2023), pp. 2001-2041.
> DOI: [https://doi.org/10.1137/21M1462994](https://doi.org/10.1137/21M1462994).
> arXiv: [https://arxiv.org/abs/2107.10760](https://arxiv.org/abs/2107.10760).

You can cite the article as

```bibtex
@article{doi:10.1137/21M1462994,
    author = {Radici, Emanuela and Stra, Federico},
    title = {Entropy solutions of mildly singular non-local scalar conservation
             laws with congestion via deterministic particle method},
    journal = {SIAM Journal on Mathematical Analysis},
    volume = {55},
    number = {3},
    pages = {2001-2041},
    year = {2023},
    doi = {10.1137/21M1462994},
    URL = {https://doi.org/10.1137/21M1462994},
    eprint = {2107.10760},
    eprinttype = {arxiv},
    eprintclass = {math.AP},
}
```

The convergence rate of the scheme is studied in the follow-up article (among other results)

```bibtex
@online{arxiv/2211.02450,
    author = {Marconi, Elio and Radici, Emanuela, and Stra, Federico},
    title = {Stability of quasi-entropy solutions of non-local scalar conservation laws},
    year = {2022},
    eprint = {2211.02450},
    eprinttype = {arxiv},
    eprintclass = {math.AP},
    url = {https://arxiv.org/abs/2211.02450},
}
```
