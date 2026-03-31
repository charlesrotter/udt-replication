#!/bin/bash
# UDT Full Replication — reproduces every number in the manuscript
# Expected runtime: ~10 min on GPU (V100), ~60 hrs on CPU
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

# Check data
if [ ! -d "data/external" ] || [ -z "$(ls data/external/ 2>/dev/null)" ]; then
    echo "NOTE: No data in data/external/. Run 'bash download_data.sh' first"
    echo "      to fetch observational datasets from Zenodo."
    echo "      (Scripts 12-16 require external data; others will still run.)"
    echo ""
fi

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

echo "[1/23] Vacuum phi profile and I_2..."
python3 analysis/01_vacuum_profile.py

echo "[2/23] Dirac eigenvalues (all kappa)..."
python3 analysis/02_eigenvalues.py

echo "[3/23] Source integrals..."
python3 analysis/03_sources.py

echo "[4/23] Bridge identities..."
python3 analysis/04_bridge_identities.py

echo "[5/23] alpha_EM derivation..."
python3 analysis/05_alpha_em.py

echo "[6/23] Particle masses..."
python3 analysis/06_particle_masses.py

echo "[7/23] Angular sector (leptons, pion)..."
python3 analysis/07_angular_sector.py

echo "[8/23] Neutrino mass..."
python3 analysis/08_neutrino_mass.py

echo "[9/23] PMNS mixing angles..."
python3 analysis/09_pmns_angles.py

echo "[10/23] Wedge product (exterior algebra)..."
python3 analysis/10_wedge_product.py

echo "[11/23] alpha^3 factorization..."
python3 analysis/11_alpha3_factorization.py

echo "[12/23] SNe fit..."
python3 analysis/12_sne_fit.py

echo "[13/23] BAO fit..."
python3 analysis/13_bao_fit.py

echo "[14/23] CMB peaks..."
python3 analysis/14_cmb_peaks.py

echo "[15/23] CMB spectrum..."
python3 analysis/15_cmb_spectrum.py

echo "[16/23] SPARC rotation curves..."
python3 analysis/16_sparc_remnants.py

echo "[17/23] Nuclear potential..."
python3 analysis/17_nuclear_potential.py

echo "[18/23] Weinberg angle and couplings..."
python3 analysis/18_weinberg_couplings.py

echo "[19/23] CKM matrix..."
python3 analysis/19_ckm_mixing.py

echo "[20/23] Quark masses and charges..."
python3 analysis/20_quark_masses.py

echo "[21/23] Higgs sector..."
python3 analysis/21_higgs_sector.py

echo "[22/23] su(3) closure verification..."
python3 analysis/22_su3_closure.py

echo "[23/23] Nuclear couplings..."
python3 analysis/23_nuclear_couplings.py

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
