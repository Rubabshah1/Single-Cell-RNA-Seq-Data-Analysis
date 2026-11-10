"""
02_fetch_cell_lines.py

Downloader for cell-line scRNA-seq datasets, registry-style (same pattern as
the organ fetcher): a small mapping from a short dataset name to a known
accession, plus a generic downloader per source.

Two sources, since cell-line panels aren't consistently on Census or 3CA:

  1. Broad Single Cell Portal (SCP) -- e.g. Kinker et al. 2020 (SCP542,
     ~198 CCLE cell lines, pan-cancer heterogeneity). SCP requires a logged-in
     auth code to bulk-download; there's no anonymous path. See --help-scp-auth
     for how to get one.

  2. GEO -- e.g. Gambardella et al. 2022 (GSE173634, 32 breast cancer cell
     lines, Drop-seq). Same generic tarball downloader as the organ script.

Registry entries were checked against the source papers as of writing
(see CELL_LINE_REGISTRY comments). Verify an accession yourself before
depending on it -- GEO/SCP records occasionally get superseded or split.

Usage:
    python 02_fetch_cell_lines.py --list
    python 02_fetch_cell_lines.py --dataset breast_gambardella2022
    python 02_fetch_cell_lines.py --dataset pancancer_kinker2020 --scp-auth-code XXXXXX
    python 02_fetch_cell_lines.py --source geo --accession GSE12345
    python 02_fetch_cell_lines.py --help-scp-auth

Requirements:
    pip install requests --break-system-packages
"""

import argparse
import os
import subprocess
import tarfile
import requests

TIMEOUT = (10, 60)

SCP_API_BASE = "https://singlecell.broadinstitute.org/single_cell/api/v1"

SCP_AUTH_HELP = """
SCP bulk downloads require a short-lived auth code from your Broad Single
Cell Portal account -- there's no anonymous bulk-download path.

  1. Log in at https://singlecell.broadinstitute.org (create a free account
     if needed) and open the study's page, e.g.:
     https://singlecell.broadinstitute.org/single_cell/study/SCP542
  2. Click "Download" on the study page. SCP will show a "Bulk Download"
     tab with your personal auth code (a short alphanumeric string), and/or
     a ready-made curl command containing it as ?auth_code=...
  3. Pass that code with --scp-auth-code. It's valid for roughly 30 minutes,
     so fetch it right before running this script.

Some studies (e.g. those derived from CCLE) also require accepting a terms
of use / data use agreement in the portal UI before download links resolve --
if you get 401s even with a fresh code, check the study page for that.
"""

# ---------------------------------------------------------------------------
# DATASET REGISTRY
# ---------------------------------------------------------------------------
# source: "scp" or "geo"
# accession: SCP study ID or GEO series accession

CELL_LINE_REGISTRY = {
    "pancancer_kinker2020": {
        "source": "scp",
        "accession": "SCP542",
        "description": (
            "Kinker et al. 2020, Nature Genetics -- multiplexed scRNA-seq of "
            "198 CCLE cancer cell lines across 22 cancer types."
        ),
    },
    "breast_gambardella2022": {
        "source": "geo",
        "accession": "GSE173634",
        "description": (
            "Gambardella et al. 2022, Nature Communications -- Drop-seq atlas "
            "of 35,276 cells across 32 breast cancer cell lines."
        ),
    },
    # Add more as you find them, e.g.:
    # "cellline_name": {"source": "geo"|"scp", "accession": "...", "description": "..."},
}


def list_registry():
    print(f"{'name':<26} {'source':<6} {'accession':<10} description")
    print("-" * 100)
    for name, entry in CELL_LINE_REGISTRY.items():
        print(f"{name:<26} {entry['source']:<6} {entry['accession']:<10} {entry['description']}")


# ---------------------------------------------------------------------------
# GEO
# ---------------------------------------------------------------------------

def fetch_from_geo(accession: str, outdir: str = "data") -> str:
    """
    Generic downloader for any GEO series (GSExxxxxx). Pulls the supplementary
    files tarball and extracts it.
    """
    os.makedirs(outdir, exist_ok=True)
    url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={accession}&format=file"
    tar_path = os.path.join(outdir, f"{accession}_RAW.tar")

    print(f"Downloading {accession} ...")
    r = requests.get(url, stream=True, timeout=TIMEOUT)
    r.raise_for_status()
    with open(tar_path, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)

    extract_dir = os.path.join(outdir, accession)
    os.makedirs(extract_dir, exist_ok=True)
    with tarfile.open(tar_path) as tar:
        tar.extractall(extract_dir, filter="data")

    print(f"Extracted to {extract_dir}")
    return extract_dir


# ---------------------------------------------------------------------------
# Broad Single Cell Portal
# ---------------------------------------------------------------------------

def fetch_from_scp(accession: str, auth_code: str, outdir: str = "data") -> str:
    """
    Bulk-download every file for an SCP study using the portal's documented
    curl-config endpoint: it returns a curl config (a list of signed,
    short-lived URLs + output paths) which we save and hand to `curl -K`.

    auth_code: short-lived code from your logged-in SCP account -- see
    SCP_AUTH_HELP / --help-scp-auth. Not the same as an API key; SCP has no
    anonymous bulk-download path.
    """
    os.makedirs(outdir, exist_ok=True)
    study_dir = os.path.join(outdir, accession)
    os.makedirs(study_dir, exist_ok=True)

    config_url = (
        f"{SCP_API_BASE}/bulk_download/generate_curl_config"
        f"?accessions={accession}&auth_code={auth_code}&directory=all"
    )

    print(f"Requesting download manifest for {accession} ...")
    r = requests.get(config_url, timeout=TIMEOUT)
    if r.status_code == 401:
        raise RuntimeError(
            "SCP rejected the auth code (401) -- it's likely expired (~30 min "
            "lifetime) or the study requires accepting a data use agreement "
            "in the portal UI first. Run with --help-scp-auth for how to get "
            "a fresh code."
        )
    r.raise_for_status()

    config_path = os.path.join(study_dir, "curl_config.txt")
    with open(config_path, "w") as f:
        f.write(r.text)

    if not r.text.strip():
        raise RuntimeError(
            f"SCP returned an empty download manifest for {accession}. "
            "Check the accession is correct and that you accepted any "
            "required data use agreement on the study page."
        )

    print(f"Downloading files for {accession} via curl ...")
    result = subprocess.run(
        ["curl", "-K", config_path],
        cwd=study_dir,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise RuntimeError(f"curl exited with status {result.returncode}")

    print(f"Downloaded {accession} -> {study_dir}")
    return study_dir


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

def fetch_cell_line_dataset(name: str = None, source: str = None,
                             accession: str = None, scp_auth_code: str = None,
                             outdir: str = "data") -> str:
    if name:
        if name not in CELL_LINE_REGISTRY:
            raise ValueError(
                f"'{name}' not in registry. Known datasets: "
                f"{', '.join(CELL_LINE_REGISTRY)}. "
                f"Or pass --source/--accession directly for anything else."
            )
        entry = CELL_LINE_REGISTRY[name]
        source = entry["source"]
        accession = entry["accession"]

    if not source or not accession:
        raise ValueError("Need either --dataset, or both --source and --accession.")

    if source == "geo":
        return fetch_from_geo(accession, outdir=outdir)
    elif source == "scp":
        if not scp_auth_code:
            raise ValueError(
                "SCP downloads need --scp-auth-code. Run with --help-scp-auth "
                "for how to get one."
            )
        return fetch_from_scp(accession, scp_auth_code, outdir=outdir)
    else:
        raise ValueError("source must be 'geo' or 'scp'")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cell-line scRNA-seq dataset downloader.")
    parser.add_argument("--dataset", default=None,
                        help=f"registry key, one of: {', '.join(CELL_LINE_REGISTRY)}")
    parser.add_argument("--source", default=None, choices=["scp", "geo"],
                        help="used instead of --dataset for an accession not in the registry")
    parser.add_argument("--accession", default=None,
                        help="SCP study ID or GEO series accession, used with --source")
    parser.add_argument("--scp-auth-code", default=None,
                        help="required for --source scp / registry entries with source=scp")
    parser.add_argument("--outdir", default="data")
    parser.add_argument("--list", action="store_true", help="print the registry and exit")
    parser.add_argument("--help-scp-auth", action="store_true",
                        help="print instructions for getting an SCP auth code and exit")
    args = parser.parse_args()

    if args.list:
        list_registry()
    elif args.help_scp_auth:
        print(SCP_AUTH_HELP)
    else:
        fetch_cell_line_dataset(
            name=args.dataset,
            source=args.source,
            accession=args.accession,
            scp_auth_code=args.scp_auth_code,
            outdir=args.outdir,
        )