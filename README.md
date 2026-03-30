# Unified Dilation Theory — Replication Repository

**Authors:** Charles Rotter (corresponding, c.rotter@udtphysics.com) and Anthony Watts

**Paper:** *Cosmological distance relations from a screened scalar modification of the Schwarzschild metric*

**Zenodo DOI:** [10.5281/zenodo.19324497](https://doi.org/10.5281/zenodo.19324497)
**Data DOI:** [10.5281/zenodo.19323253](https://doi.org/10.5281/zenodo.19323253)

**License:** Academic use only (see [LICENSE](LICENSE)). Commercial rights reserved under pending patents.

## Overview

This repository contains everything needed to reproduce every numerical result in the manuscript. One metric, three quantum numbers, zero fitted physics parameters.

The UDT metric:

$$ds^2 = -e^{-2\phi(r)}c^2 dt^2 + e^{2\phi(r)}dr^2 + r^2 d\Omega^2$$

where $\phi(r)$ satisfies a nonlinear screened Klein-Gordon equation, produces **34 derived quantities** (all 24 Standard Model free parameters plus 10 QCD observables) from three quantum numbers $(j, \ell, |\kappa_\text{max}|) = (1/2, 1, 3)$ and the electron mass:

- **Fine structure constant** $1/\alpha = 137.036$ (+0.0004%, three-level chain)
- **Particle masses**: $m_p/m_e = 6\pi^5$ (-0.002%), $m_\tau$ from Koide $Z_3$ (+0.007%), 15+ hadrons sub-3%
- **PMNS mixing**: all 3 angles + $\delta_\text{CP}$ within $1\sigma$ of NuFIT (0 free parameters)
- **CKM mixing**: $\sin\theta_C = 9/40$ (-0.13%), $\delta_\text{CKM} = \arctan(9/4)$ (0.1$\sigma$), $\theta_\text{QCD} = 0$
- **Higgs sector**: $v = 504\pi^6 m_e$ (+0.6%), $m_H = 126.7$ GeV (+1.1%)
- **Nuclear force**: $g_A = 4/\pi$ (-0.23%), $g^2/(4\pi) = \pi^7/225$ (-0.6%), su(3) closure at machine precision
- **Neutrino mass** $m_\nu = \alpha^3 m_e/4 = 0.049$ eV (-0.8%)
- **CMB three-spectrum** (TT/TE/EE): EE at 9.3% outperforms $\Lambda$CDM Planck EE (12.6%), zero parameters vs six
- **SNe Ia** at 0.166 mag RMS, **BAO** at 3.8% RMS
- **Universe**: $r_* = 9.164$ Gpc, $c^2 = 2GM/r_*$ (Machian closure, 0.0000%)

External inputs: $c$, $G$, $m_e$, $T_\text{CMB}$, $T_\text{starlight}$, $\hbar$ (conversion only).

## Quick Start

```bash
# Clone and setup
git clone <repo-url>
cd udt-repo
pip install -r requirements.txt

# Download observational data from Zenodo
bash download_data.sh

# Run full replication (~60 min GPU, ~5 hrs CPU)
bash run_all.sh

# Or run individual analyses
python3 analysis/01_vacuum_profile.py
python3 analysis/02_eigenvalues.py
# ... etc.
```

## Requirements

- Python 3.10+
- NumPy 1.26+
- SciPy 1.13+
- Matplotlib 3.9+
- PyTorch 2.3+ (optional, for GPU acceleration)

### GPU Acceleration

Scripts 02, 08, and 10 use GPU acceleration for eigenvalue scanning if a CUDA-capable GPU is available. All include automatic CPU fallback:

```python
import torch
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
```

### Runtime Estimates

| Component | GPU (V100) | CPU (single core) |
|-----------|-----------|-------------------|
| Vacuum profile (01) | 5 sec | 30 sec |
| All eigenvalues (02) | 2 min | 30 min |
| Source integrals (03) | 3 min | 40 min |
| Cosmological fits (12-16) | 1 min | 5 min |
| New analyses (18-23) | 1 min | 1 min |
| Figure generation | 2 min | 2 min |
| **Total** | **~15 min** | **~2 hrs** |

## Repository Structure

```
├── README.md                    # This file
├── requirements.txt             # Pinned dependencies
├── LICENSE                      # Academic-use license
├── run_all.sh                   # Master replication script (23 steps)
├── download_data.sh             # Fetch datasets from Zenodo
├── manuscript/
│   ├── manuscript.tex           # LaTeX primary (20 pages)
│   ├── manuscript.pdf           # Compiled PDF
│   ├── figures/                 # 13 figures (.pdf, .png)
│   └── figure_sources.md        # Figure → script cross-reference
├── data/
│   ├── external/                # Observational datasets (via download_data.sh)
│   └── generated/               # All computed outputs (JSON)
├── lib/
│   ├── constants.py             # Locked parameters, physical constants
│   ├── vacuum_phi.py            # Vacuum φ ODE solver
│   ├── dirac_formT.py           # Dirac Form-T eigenvalue solver (GPU)
│   ├── angular_integrals.py     # S² coupling integrals, mass formulas
│   └── utils.py                 # I/O, comparison, plotting utilities
├── analysis/                    # 23 analysis scripts
│   ├── 01_vacuum_profile.py     # Solve vacuum φ, compute I₂
│   ├── 02_eigenvalues.py        # Dirac eigenvalues for κ = ±1,±2,±3
│   ├── 03_sources.py            # Source integrals
│   ├── 04_bridge_identities.py  # Bridge identity verification
│   ├── 05_alpha_em.py           # α_EM three-level derivation chain
│   ├── 06_particle_masses.py    # Full mass table
│   ├── 07_angular_sector.py     # Lepton masses, pion, tau
│   ├── 08_neutrino_mass.py      # m_ν = α³mₑ/4
│   ├── 09_pmns_angles.py        # PMNS mixing angles
│   ├── 10_wedge_product.py      # Exterior algebra structure
│   ├── 11_alpha3_factorization.py # α³ factorization verification
│   ├── 12_sne_fit.py            # Type Ia supernovae
│   ├── 13_bao_fit.py            # Baryon acoustic oscillations
│   ├── 14_cmb_peaks.py          # CMB peak positions
│   ├── 15_cmb_spectrum.py       # Full CMB power spectrum (TT/TE/EE)
│   ├── 16_sparc_remnants.py     # SPARC rotation curves
│   ├── 17_nuclear_potential.py  # Nuclear potential and phase shifts
│   ├── 18_weinberg_couplings.py # Weinberg angle and coupling hierarchy
│   ├── 19_ckm_mixing.py         # CKM matrix and quark mixing
│   ├── 20_quark_masses.py       # Quark charges, color, mass ratios
│   ├── 21_higgs_sector.py       # Higgs VEV, self-coupling, boson mass
│   ├── 22_su3_closure.py        # su(3) algebra closure verification
│   └── 23_nuclear_couplings.py  # g_A, f_π, g_πNN, proton radius
├── figures/                     # 12 figure generation scripts
└── supplements/                 # 17 supplementary documents (LaTeX)
    ├── S1–S12                   # Original supplements
    ├── S13_ckm_derivation.tex   # CKM from vertex partition
    ├── S14_quark_charges.tex    # Color, charges, mass ratios
    ├── S15_higgs_mechanism.tex  # VEV as angular capacity
    ├── S16_qcd_nuclear.tex      # su(3) closure, nuclear force
    └── S17_cmb_three_spectrum.tex # Kaleidoscope Principle, 5 experiments
```

## Locked Parameters

All three geometric parameters are derived from the metric. The algebraic structure exists ONLY at this triple.

| Parameter | Value | Origin |
|-----------|-------|--------|
| $\phi_0$ | $-\cos(\pi/5) = -0.80902$ | Geometric: $\langle G \rangle = 2/\pi$ |
| $\mu^2$ | $\pi/3$ | Geometric: $\pi \times \langle\cos^2\theta\rangle$ on $S^2$ |
| $r_*$ | $7 - 1/80 = 6.9875$ | Derived: $(2|\kappa_\max|+1) - \text{source}^2/(2|\kappa_\max|-1)$ |
| $C$ | $4\pi^2 m_e r_* = 140.96$ MeV | Electron anchor |

## Consistency Gates

Every computation passes internal consistency checks before results are used:

| Gate | Threshold | Checked in |
|------|-----------|------------|
| $\Box_g\phi - \mu^2\phi = 0$ (vacuum) | $< 10^{-8}$ | `01_vacuum_profile.py` |
| $G'(r_*) = 0$ (Neumann BC) | $< 10^{-10}$ | `02_eigenvalues.py` |
| SL normalization | $< 10^{-8}$ | `03_sources.py` |
| Bridge identity | $< 0.3\%$ | `04_bridge_identities.py` |
| su(3) closure | $< 10^{-14}$ | `22_su3_closure.py` |
| Casimir $C_2 = 4/3$ | $< 10^{-14}$ | `22_su3_closure.py` |

## Citation

If you use this code or results, please cite:

```bibtex
@article{rotter2026udt,
  title={Cosmological distance relations from a screened scalar
         modification of the {Schwarzschild} metric},
  author={Rotter, Charles and Watts, Anthony},
  year={2026},
  doi={10.5281/zenodo.19324497}
}
```

## Contact

Charles Rotter: c.rotter@udtphysics.com
