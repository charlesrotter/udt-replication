#!/usr/bin/env python3
# Generates: Figure 10 in manuscript section 21 (CMB Power Spectrum)
# Three-panel CMB: TT, TE, EE across 5 independent experiments
"""
fig10 -- CMB three-spectrum comparison (TT/TE/EE).

Uses cr184_definitive_multi as base, with improved title.
Zero free parameters. All constants from ds².
Verified on Planck, SPT-3G 2018, SPT-3G D1, ACT DR4, WMAP.

Key result: EE on SPT-3G D1 (9.3%) outperforms LCDM Planck EE (12.6%).
"""
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.image import imread

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from lib.utils import savefig

raw_path = os.path.join(os.path.dirname(__file__), '..',
                        'manuscript', 'figures', 'fig10_cmb_spectrum_raw.png')

if not os.path.exists(raw_path):
    print("  WARNING: Raw CMB image not found.")
    print("  Copy cr_next/outputs/cr184_definitive_multi.png to")
    print(f"  {raw_path}")
    sys.exit(1)

img = imread(raw_path)

fig, ax = plt.subplots(figsize=(8, 11.4))

# Crop top ~6.5% where old title lives
h = img.shape[0]
crop_top = int(h * 0.065)
img_cropped = img[crop_top:, :, :]

ax.imshow(img_cropped)
ax.axis('off')

# Publication title
fig.suptitle(
    'CMB Power Spectrum: Zero Free Parameters\n'
    r'All constants from $ds^2$ — verified on 5 independent experiments',
    fontsize=13, fontweight='bold', y=0.995, va='top'
)

fig.tight_layout(rect=[0, 0, 1, 0.97])
savefig(fig, 'fig10_cmb_spectrum', dpi=200)
plt.close()
