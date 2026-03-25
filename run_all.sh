#!/bin/bash
# UDT Full Replication — reproduces every number in the manuscript
# Expected runtime: ~45 min on GPU (V100), ~4 hrs on CPU
# Requirements: Python 3.10+, NumPy, SciPy, PyTorch (CUDA optional), Matplotlib
#
# Authors: Charles Rotter and Anthony Watts
# License: Academic use only (see LICENSE)

set -e  # Exit on any error

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== UDT Full Replication ==="
echo "Start: $(date)"
echo "Working directory: $SCRIPT_DIR"
echo ""

# Check Python
python3 -c "import numpy, scipy; print(f'NumPy {numpy.__version__}, SciPy {scipy.__version__}')" || {
    echo "ERROR: NumPy/SciPy not found. Install: pip install -r requirements.txt"
    exit 1
}

# Check GPU (optional)
python3 -c "import torch; print(f'PyTorch {torch.__version__}, CUDA: {torch.cuda.is_available()}')" 2>/dev/null || {
    echo "NOTE: PyTorch not found or no CUDA. Running in CPU mode (slower)."
}

echo ""
echo "=== Analysis Pipeline ==="
echo ""

echo "[1/17] Vacuum phi profile and I_2..."
python3 analysis/01_vacuum_profile.py

echo "[2/17] Dirac eigenvalues (all kappa)..."
python3 analysis/02_eigenvalues.py

echo "[3/17] Source integrals..."
python3 analysis/03_sources.py

echo "[4/17] Bridge identities..."
python3 analysis/04_bridge_identities.py

echo "[5/17] alpha_EM derivation..."
python3 analysis/05_alpha_em.py

echo "[6/17] Particle masses..."
python3 analysis/06_particle_masses.py

echo "[7/17] Angular sector (leptons, pion)..."
python3 analysis/07_angular_sector.py

echo "[8/17] Neutrino mass..."
python3 analysis/08_neutrino_mass.py

echo "[9/17] PMNS mixing angles..."
python3 analysis/09_pmns_angles.py

echo "[10/17] Wedge product (exterior algebra)..."
python3 analysis/10_wedge_product.py

echo "[11/17] alpha^3 factorization..."
python3 analysis/11_alpha3_factorization.py

echo "[12/17] SNe fit..."
python3 analysis/12_sne_fit.py

echo "[13/17] BAO fit..."
python3 analysis/13_bao_fit.py

echo "[14/17] CMB peaks..."
python3 analysis/14_cmb_peaks.py

echo "[15/17] CMB spectrum..."
python3 analysis/15_cmb_spectrum.py

echo "[16/17] SPARC rotation curves..."
python3 analysis/16_sparc_remnants.py

echo "[17/17] Nuclear potential..."
python3 analysis/17_nuclear_potential.py

echo ""
echo "=== Generating Figures ==="
echo ""

for fig in figures/fig*.py; do
    echo "  $(basename $fig)..."
    python3 "$fig" || echo "  WARNING: $fig failed (may need data from analysis)"
done

echo ""
echo "=== Replication Complete ==="
echo "End: $(date)"
echo ""
echo "All numerical outputs:  data/generated/"
echo "All figures:            manuscript/figures/"
echo ""
echo "To verify: compare data/generated/*.json values against manuscript claims."
