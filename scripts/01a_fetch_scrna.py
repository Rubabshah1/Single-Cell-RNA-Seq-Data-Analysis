"""
01_fetch_dataset.py

Downloader for 3CA (Curated Cancer Cell Atlas) scRNA-seq studies.

3CA has no public API. It's a set of per-cancer-type category pages at a fixed
URL pattern, each listing studies in an HTML table with per-study
"Download Data" / "Download Meta-data" links. This scrapes that table.

NOTE ON RELIABILITY: the category page HTML structure was confirmed live
(see e.g. https://www.weizmann.ac.il/sites/3CA/prostate), but the exact href
pattern behind each "Download Data" link was NOT independently verified. The
parser below looks for <a> tags whose visible text contains "Download Data
about" / "Download Meta-data about" -- inspect one category page's HTML source
yourself (browser devtools -> Elements/Network) before relying on this in
production, and adjust the selector if it doesn't match. Run with
--list-3ca --debug-3ca-links first if data_url/meta_url come back None.

Usage:
    python 01_fetch_dataset.py --organ lung --list-3ca
    python 01_fetch_dataset.py --organ lung --list-3ca --debug-3ca-links
    python 01_fetch_dataset.py --organ lung
    python 01_fetch_dataset.py --organ breast --study "Chen et al. 2021"

Requirements:
    pip install pandas requests beautifulsoup4 --break-system-packages
"""

import argparse
import os
import re
import requests
import pandas as pd

THREECA_BASE_URL = "https://www.weizmann.ac.il/sites/3CA"

# Network timeouts (connect, read) in seconds -- without these a hung server
# stalls the script indefinitely.
TIMEOUT = (10, 60)


# ---------------------------------------------------------------------------
# ORGAN -> 3CA CATEGORY SLUG
# ---------------------------------------------------------------------------
# CONFIRMED slugs (checked live): "prostate", "hematologic", "brain".
# The rest are inferred from the site's known cancer-type categories and NEED
# VERIFICATION against the live site (visit weizmann.ac.il/sites/3CA and check
# each category link) before you rely on them -- a wrong slug will just 404,
# which fetch_3ca_catalog() will surface as an empty result.

ORGAN_REGISTRY = {
    "bladder":      "bladder",
    "blood":        "hematologic",      # confirmed
    "brain":        "brain",            # confirmed
    "breast":       "breast",
    "colorectal":   "colorectal",
    "esophagus":    "esophageal",
    "head_and_neck": "head_and_neck",
    "kidney":       "kidney",
    "liver":        "liver",
    "lung":         "lung",
    "ovary":        "ovarian",
    "pancreas":     "pancreatic",
    "prostate":     "prostate",         # confirmed
    "rectum":       "colorectal",       # likely folded into colorectal
    "skin":         "skin",
    "soft_tissue":  "sarcoma",
    "stomach":      "gastric",
    "uterus":       "endometrial",
}


def fetch_3ca_catalog(category_slug: str, debug: bool = False) -> "pd.DataFrame":
    """
    Scrape a 3CA category page into a table of studies: title, disease,
    technology, #samples, #cells, and the (best-guess) download links.
    Use this FIRST to see what's actually in a category before downloading.

    debug=True: also stores every <a> tag found in each row (text + href + any
    onclick/data-* attributes) in a 'raw_links' column, so you can see exactly
    what's on the page when the normal link-detection comes up empty (e.g.
    JS-triggered downloads, hrefs on a wrapping element, etc.).
    """
    from bs4 import BeautifulSoup
    import pandas as pd

    url = f"{THREECA_BASE_URL}/{category_slug}"
    resp = requests.get(url, timeout=TIMEOUT)
    resp.raise_for_status()
    soup = BeautifulSoup(resp.text, "html.parser")

    rows = []
    for row in soup.select("table tr"):
        cells = row.find_all("td")
        if not cells:
            continue
        row_text = [c.get_text(strip=True) for c in cells]

        # Use substring match against get_text() (handles nested <span> etc.)
        # rather than an exact `string=` match, which silently fails on any
        # markup nested inside the <a> tag.
        anchors = row.find_all("a")
        data_link = next((a for a in anchors if "Download Data" in a.get_text()), None)
        meta_link = next((a for a in anchors if "Download Meta-data" in a.get_text()), None)

        entry = {
            "raw_row": row_text,
            "data_url": data_link.get("href") if data_link else None,
            "meta_url": meta_link.get("href") if meta_link else None,
        }
        if debug:
            entry["raw_links"] = [
                {"text": a.get_text(strip=True), "href": a.get("href"),
                 "attrs": {k: v for k, v in a.attrs.items() if k != "href"}}
                for a in anchors
            ]
        rows.append(entry)

    df = pd.DataFrame(rows)
    if df.empty:
        print(f"No rows parsed for category '{category_slug}' -- check the slug is correct: {url}")
    return df


def _resolve_url(row, key: str, link_text: str):
    """Get a download URL from the parsed column, falling back to raw_links."""
    url = row.get(key)
    if not url and "raw_links" in row and row["raw_links"]:
        for link in row["raw_links"]:
            if link_text in link["text"]:
                url = link["href"]
                break
    if url and not url.startswith("http"):
        url = THREECA_BASE_URL + url
    return url


def _download(url: str, dest_dir: str, stem: str) -> str:
    """Stream a file to disk, preserving the server's extension where possible."""
    ext = os.path.splitext(url.split("?")[0])[1] or ".tar.gz"
    path = os.path.join(dest_dir, stem + ext)

    if os.path.exists(path):
        print(f"Already exists: {path}")
        return path

    print(f"Downloading: {url}")
    r = requests.get(url, stream=True, timeout=TIMEOUT)
    r.raise_for_status()
    with open(path, "wb") as f:
        for chunk in r.iter_content(chunk_size=1024 * 1024):
            if chunk:
                f.write(chunk)

    print(f"Saved: {path}")
    return path


def fetch_from_3ca(category_slug: str, study_match: str = None,
                   outdir: str = "data") -> str:
    """
    Download one or all studies from a 3CA category.

    For every study, downloads BOTH expression/count data and metadata:

        outdir/
            Study_Name/
                data.tar.gz
                metadata.tar.gz
    """
    os.makedirs(outdir, exist_ok=True)

    catalog = fetch_3ca_catalog(category_slug, debug=True)

    if study_match and not catalog.empty:
        catalog = catalog[
            catalog["raw_row"].apply(lambda r: any(study_match in c for c in r))
        ]

    if catalog.empty:
        print(f"No studies found for '{category_slug}'.")
        return outdir

    saved = []

    for _, row in catalog.iterrows():
        study_name = row["raw_row"][0] if row["raw_row"] else "unknown_study"
        safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", study_name).strip("_")

        study_dir = os.path.join(outdir, safe_name)
        os.makedirs(study_dir, exist_ok=True)

        print("\n" + "=" * 70)
        print(f"Study: {study_name}")
        print("=" * 70)

        for key, label, stem in [
            ("data_url", "Download Data", "data"),
            ("meta_url", "Download Meta-data", "metadata"),
        ]:
            url = _resolve_url(row, key, label)
            if url:
                saved.append(_download(url, study_dir, stem))
            else:
                print(f"WARNING: No {stem} URL found.")

    print("\n" + "=" * 70)
    print(f"Finished. Saved {len(saved)} file(s) from 3CA category '{category_slug}'.")
    print("=" * 70)

    return outdir


def fetch_organ(organ: str, study_match: str = None, outdir: str = "data") -> str:
    if organ not in ORGAN_REGISTRY:
        raise ValueError(
            f"Organ '{organ}' not in registry. Add it to ORGAN_REGISTRY with its "
            f"3CA category slug. Known organs: {', '.join(sorted(ORGAN_REGISTRY))}"
        )
    return fetch_from_3ca(ORGAN_REGISTRY[organ], study_match=study_match, outdir=outdir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="3CA scRNA-seq dataset downloader.")
    parser.add_argument("--organ", required=True,
                        help=f"one of: {', '.join(sorted(ORGAN_REGISTRY))}")
    parser.add_argument("--study", default=None,
                        help="substring to match one study within a category, e.g. "
                             "'Chen et al. 2021'; omit to download every study")
    parser.add_argument("--list-3ca", action="store_true",
                        help="just print the 3CA catalog table for this organ and exit")
    parser.add_argument("--debug-3ca-links", action="store_true",
                        help="with --list-3ca: dump every <a> tag (text/href/attrs) per row, "
                             "to diagnose why data_url/meta_url came back None")
    parser.add_argument("--outdir", default="data")
    args = parser.parse_args()

    if args.organ not in ORGAN_REGISTRY:
        raise SystemExit(
            f"Organ '{args.organ}' not in registry. "
            f"Known organs: {', '.join(sorted(ORGAN_REGISTRY))}"
        )

    if args.list_3ca:
        cat = fetch_3ca_catalog(ORGAN_REGISTRY[args.organ], debug=args.debug_3ca_links)
        if args.debug_3ca_links and not cat.empty:
            for i, row in cat.iterrows():
                print(f"\n--- Row {i}: {row['raw_row'][0] if row['raw_row'] else '?'} ---")
                for link in row["raw_links"]:
                    print(" ", link)
        else:
            print(cat)
    else:
        fetch_organ(args.organ, study_match=args.study, outdir=args.outdir)