# Unified Dilation Theory — Replication Repository

**Authors:** Charles Rotter and Anthony Watts

**Paper:** *Geodesic Implications of a Screened Scalar Field in General Relativity: Mass Emergence, Coupling Constants, and Cosmological Predictions from a Single Metric*

**License:** Academic use only (see [LICENSE](LICENSE)). Commercial rights reserved under pending patents.

## Overview

This repository contains everything needed to reproduce every numerical result in the manuscript. One metric, three quantum numbers, zero fitted physics parameters.

The UDT metric:

$$ds^2 = -e^{-2\phi(r)}c^2 dt^2 + e^{2\phi(r)}dr^2 + r^2 d\Omega^2$$

where $\phi(r)$ satisfies a nonlinear screened Klein-Gordon equation, produces:
- **18 particle masses** at sub-percent accuracy (4 parameter-free from angular sector)
- **Fine structure constant** $1/\alpha = 36\pi/I_2 = 137.4$ (+0.29%)
- **Neutrino mass** $m_\nu = \alpha^3 m_e/4 = 0.049$ eV (-0.6%)
- **PMNS mixing angles** all within 1$\sigma$ of NuFIT (0 free parameters)
- **CMB peak positions** at 1.32% RMS
- **SNe Ia Hubble diagram** at 0.166 mag RMS
- **BAO** at 3.8% RMS

External inputs: $c$, $G$, $m_e$, $T_\text{CMB}$, $T_\text{starlight}$, $\hbar$ (conversion only).

## Quick Start

```bash
# Clone and setup
git clone <repo-url>
cd udt-repo
pip install -r requirements.txt

# Run full replication (45 min GPU, 4 hrs CPU)
bash run_all.sh

# Or run individual analyses
python analysis/01_vacuum_profile.py
python analysis/02_eigenvalues.py
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
| Figure generation | 2 min | 2 min |
| **Total** | **~15 min** | **~2 hrs** |

## Repository Structure

```
├── README.md                    # This file
├── requirements.txt             # Pinned dependencies
├── LICENSE                      # Academic-use license
├── run_all.sh                   # Master replication script
├── manuscript/
│   ├── manuscript.tex           # LaTeX primary
│   ├── manuscript.md            # Markdown companion
│   ├── manuscript.pdf           # Compiled PDF
│   ├── figures/                 # Generated figures (.pdf, .png)
│   └── figure_sources.md        # Figure → script cross-reference
├── data/
│   ├── external/                # External datasets with provenance
│   │   ├── pantheon_plus/       # Pantheon+ SNe data
│   │   ├── planck_2018/         # Planck CMB data
│   │   ├── boss_bao/            # BOSS BAO measurements
│   │   ├── desi_bao/            # DESI BAO measurements
│   │   ├── sparc/               # SPARC rotation curves
│   │   └── nufit/               # NuFIT neutrino data
│   └── generated/               # All computed outputs (JSON)
├── lib/
│   ├── constants.py             # Locked parameters, physical constants
│   ├── vacuum_phi.py            # Vacuum φ ODE solver
│   ├── dirac_formT.py           # Dirac Form-T eigenvalue solver (GPU)
│   ├── scalar_perturbation.py   # Scalar boson solver
│   ├── vector_perturbation.py   # EM vector mode solver
│   ├── angular_integrals.py     # S² coupling integrals
│   └── utils.py                 # I/O, comparison, plotting utilities
├── analysis/
│   ├── 01_vacuum_profile.py     # Solve vacuum φ, compute I₂
│   ├── 02_eigenvalues.py        # Dirac eigenvalues for κ = ±1,±2,±3
│   ├── 03_sources.py            # Source integrals
│   ├── 04_bridge_identities.py  # Bridge identity verification
│   ├── 05_alpha_em.py           # α_EM derivation
│   ├── 06_particle_masses.py    # Full mass table
│   ├── 07_angular_sector.py     # Lepton masses, pion
│   ├── 08_neutrino_mass.py      # m_ν = α³mₑ/4
│   ├── 09_pmns_angles.py        # PMNS mixing angles
│   ├── 10_wedge_product.py      # Exterior algebra structure
│   ├── 11_alpha3_factorization.py # α³ factorization verification
│   ├── 12_sne_fit.py            # Type Ia supernovae
│   ├── 13_bao_fit.py            # Baryon acoustic oscillations
│   ├── 14_cmb_peaks.py          # CMB peak positions
│   ├── 15_cmb_spectrum.py       # Full CMB power spectrum
│   ├── 16_sparc_remnants.py     # SPARC rotation curves
│   └── 17_nuclear_potential.py  # Nuclear potential
├── figures/
│   ├── fig01_sne_hubble.py      # Hubble diagram
│   ├── fig02_bao.py             # BAO D_V/r_d
│   ├── fig03_mass_ratios.py     # Lepton mass ratios
│   ├── fig04_angular_sector.py  # Angular sector schematic
│   ├── fig05_mass_ladder.py     # Mass ladder
│   ├── fig06_neutrino_scale.py  # Neutrino mass scale
│   ├── fig07_pmns_contours.py   # PMNS contour plots
│   ├── fig08_exterior_algebra.py # Exterior algebra diagram
│   ├── fig09_connection_map.py  # Unified connection map
│   ├── fig10_cmb_spectrum.py    # CMB power spectrum
│   ├── fig11_sparc_rotation.py  # SPARC rotation curve
│   └── fig12_falsification_roadmap.py # Falsification timeline
└── supplements/
    ├── S1_metric_derivation.tex
    ├── S2_dirac_formT.tex
    ├── S3_angular_sector.tex
    ├── S4_pion_diophantine.tex
    ├── S5_alpha_em.tex
    ├── S6_neutrino_derivation.tex
    ├── S7_pmns_derivation.tex
    ├── S8_bridge_identities.tex
    ├── S9_orbit_matching.tex
    ├── S10_cosmological_profile.tex
    ├── S11_unified_table.tex
    └── S12_predictions.tex
```

## Locked Parameters

All three geometric parameters are derived from the metric. The algebraic structure exists ONLY at this triple.

| Parameter | Value | Origin |
|-----------|-------|--------|
| $\phi_0$ | $-\cos(\pi/5) = -0.80902$ | Geometric: $\langle G \rangle = 2/\pi$ |
| $\mu^2$ | $\pi/3$ | Geometric: $\pi \times \langle\cos^2\theta\rangle$ on $S^2$ |
| $r_*$ | $7 - 1/80 = 6.9875$ | Derived: $(2|\kappa_\max|+1) - \text{source}^2/(2|\kappa_\max|-1)$ |
| $C$ | $4\pi^2 m_e r_* = 140.95$ MeV | Electron anchor |

## Consistency Gates

Every computation passes internal consistency checks before results are used:

| Gate | Threshold | Checked in |
|------|-----------|------------|
| $\Box_g\phi - \mu^2\phi = 0$ (vacuum) | $< 10^{-8}$ | `01_vacuum_profile.py` |
| $G'(r_*) = 0$ (Neumann BC) | $< 10^{-10}$ | `02_eigenvalues.py` |
| SL normalization | $< 10^{-8}$ | `03_sources.py` |
| Bridge identity | $< 0.1\%$ | `04_bridge_identities.py` |

## Citation

If you use this code or results, please cite:

```bibtex
@article{rotter2026udt,
  title={Geodesic Implications of a Screened Scalar Field in General Relativity},
  author={Rotter, Charles and Watts, Anthony},
  year={2026},
  note={Zenodo DOI: [pending]}
}
```

## Contact

For questions about the theory, code, or commercial licensing: [contact information pending]
