# A uniform lower bound for the Faltings height of Jacobians of Shimura curves X_0^D(M)

[![verify](https://github.com/Li-Hongmin/shimura-jacobian-faltings-height/actions/workflows/verify.yml/badge.svg)](https://github.com/Li-Hongmin/shimura-jacobian-faltings-height/actions/workflows/verify.yml) [![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.24054785-blue)](https://doi.org/10.5281/zenodo.24054785)

**Hongmin Li** (Institute of Science Tokyo; The University of Tokyo) · ORCID [0000-0003-0228-0600](https://orcid.org/0000-0003-0228-0600)

Main result of the manuscript: for squarefree N = DM prime to 6 and genus g ≥ 2, the Faltings height of the Jacobian of the Shimura curve X₀^D(M) satisfies a uniform lower bound of order (1/3 − o(1)) g log N.

- Manuscript on rStation: <https://rstation.ai/problems/math/papers/4/v5> (version 5)
- Archived on Zenodo: concept DOI [10.5281/zenodo.24054785](https://doi.org/10.5281/zenodo.24054785) (all versions); this repository mirrors version 5, DOI [10.5281/zenodo.24054831](https://doi.org/10.5281/zenodo.24054831)

## Contents

- `paper/paper.md`, `paper/shimura-jacobian-faltings-height.pdf` — the manuscript (version 5).
- `verification/` — the supplement of that version, unpacked unchanged; its own README describes each file.
- `.github/workflows/verify.yml` — reruns every check below on each push.

## What the programs check

| Check | Command (run inside `verification/`) |
|---|---|
| Finite arithmetic checks | `python checks.py` |

The programs check the finite and exact computations used in the manuscript; the mathematical arguments are those written in the paper.

## Reproducing

```sh
cd verification
python checks.py
```

Requirements: Python 3.9+. Run the commands from a copy if you want to keep the recorded outputs untouched.

## Cite

```bibtex
@misc{li2026shimurajacobianfaltingsh,
  author    = {Li, Hongmin},
  title     = {A uniform lower bound for the Faltings height of Jacobians of Shimura curves X_0^D(M)},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.24054785},
  url       = {https://doi.org/10.5281/zenodo.24054785}
}
```

## License

Manuscript and data: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Code: MIT (see `LICENSE`).
