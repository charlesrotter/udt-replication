# Pantheon+ Supernova Data

## Source
- **Paper:** Scolnic et al. (2022), "The Pantheon+ Analysis: The Full Data Set and Light-curve Release"
- **URL:** https://github.com/PantheonPlusSH0ES/DataRelease
- **Download date:** [To be filled at data download time]
- **Version:** Pantheon+ (1701 SNe Ia)

## Files
- `Pantheon+SH0ES.dat` — Full dataset (z, distance modulus, uncertainties)

## Notes
- Data is model-agnostic: distance moduli from light curve fitting only
- No LCDM assumptions in the data compilation
- If this file is not present, `analysis/12_sne_fit.py` uses representative bin averages

## Citation
```bibtex
@article{Scolnic2022,
  author = {Scolnic, D. and others},
  title = {The Pantheon+ Analysis},
  journal = {ApJ},
  year = {2022},
  volume = {938},
  pages = {113}
}
```
