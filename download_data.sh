#!/bin/bash
# download_data.sh — Fetch all datasets from Zenodo for UDT replication
#
# Source: https://zenodo.org/records/19323253
# DOI: 10.5281/zenodo.19323253
# Authors: Charles Rotter, Anthony Watts
#
# Total download: ~101 MB (6 zip files)
# Datasets: Planck, Pantheon+SH0ES, ACT, SPT-3G, WMAP, SPARC

set -e

ZENODO_BASE="https://zenodo.org/records/19323253/files"
DATA_DIR="data/external"

mkdir -p "$DATA_DIR"

echo "=== UDT Replication Data Download ==="
echo "Source: https://zenodo.org/records/19323253"
echo "Target: $DATA_DIR/"
echo ""

download_and_extract() {
    local filename="$1"
    local desc="$2"
    local url="${ZENODO_BASE}/${filename}?download=1"
    local target="${DATA_DIR}/${filename}"

    echo -n "[$desc] Downloading ${filename}... "
    if [ -f "$target" ]; then
        echo "already exists, skipping download."
    else
        curl -sL -o "$target" "$url"
        echo "done."
    fi

    echo -n "  Extracting... "
    unzip -qo "$target" -d "$DATA_DIR/"
    echo "done."
}

download_and_extract "planck.zip"          "1/6  Planck 2018 CMB (60.3 MB)"
download_and_extract "Pantheon+SH0ES.zip"  "2/6  Pantheon+ SNe (10.6 MB)"
download_and_extract "act.zip"             "3/6  ACT DR4 CMB (29.7 MB)"
download_and_extract "spt3g.zip"           "4/6  SPT-3G CMB (79.1 kB)"
download_and_extract "wmap.zip"            "5/6  WMAP CMB (30.9 kB)"
download_and_extract "SPARC.zip"           "6/6  SPARC rotation curves (221.8 kB)"

echo ""
echo "=== Download complete ==="
echo "Contents of $DATA_DIR/:"
ls -la "$DATA_DIR/"
echo ""
echo "To verify checksums, run: cd $DATA_DIR && sha256sum *.zip"
