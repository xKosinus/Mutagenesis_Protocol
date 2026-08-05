"""
primer_io.py
============
All file I/O for the Site-Directed Mutagenesis Primer Designer.

Reads and writes:
  - primer_list.csv / primer_temp.csv      (flat per-primer rows)
  - primer_list.txt                        (human-readable table)
  - primer_list.json                       (structured, dual-section format)
  - slice_primers.csv                      (SLiCE assembly primers)
  - variant_databank.json                  (variant → mutation-string map)
  - wildtype_sequences.json                (stored WT sequences)
  - protocol_data.json                     (protocol steps for GUI display)

No tkinter / GUI imports.  Callers pass plain Python objects; this module
handles all path logic, encoding, and format details.
"""

from __future__ import annotations

import copy
import csv
import json
import os
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from primer_core import SLiCEPrimer


# ---------------------------------------------------------------------------
# Default file/directory names
# ---------------------------------------------------------------------------
DEFAULT_DIR          = Path.home() / "Mutagenesis"
PRIMER_LIST_CSV      = "primer_list.csv"
PRIMER_TEMP_CSV      = "primer_temp.csv"
PRIMER_LIST_TXT      = "primer_list.txt"
PRIMER_LIST_JSON     = "primer_list.json"
SLICE_PRIMERS_CSV    = "slice_primers.csv"
DATABANK_JSON        = "variant_databank.json"
WT_SEQUENCES_JSON    = "wildtype_sequences.json"
PROTOCOL_DATA_JSON   = "protocol_data.json"

# CSV column names (canonical — used for both reading and writing)
CSV_COLUMNS = [
    "Primer Name",
    "Primer Sequence",
    "Length",
    "Tm (C)",
    "GC Content (%)",
    "Overlap Length",
    "Overlap Tm",
    "Overlap GC (%)",
    "Mutations",
]


# ===========================================================================
# Primer CSV helpers
# ===========================================================================

def _primer_set_to_rows(entry: dict) -> List[dict]:
    """
    Convert one PrimerSet dict (pair-level) into two flat CSV row dicts
    (one forward, one reverse).
    """
    all_muts   = ",".join(entry["all_covered_mutations"])
    base_name  = "_".join(entry["all_covered_mutations"])
    overlap_gc = entry.get("overlap_gc_content", 0.0)

    return [
        {
            "Primer Name":    f"{base_name}_for",
            "Primer Sequence": entry["forward_primer"],
            "Length":          entry["forward_length"],
            "Tm (C)":          entry["forward_tm"],
            "GC Content (%)":  entry.get("forward_gc_content", 0.0),
            "Overlap Length":  entry["overlap_length"],
            "Overlap Tm":      entry["overlap_tm"],
            "Overlap GC (%)":  overlap_gc,
            "Mutations":       all_muts,
        },
        {
            "Primer Name":    f"{base_name}_rev",
            "Primer Sequence": entry["reverse_primer"],
            "Length":          entry["reverse_length"],
            "Tm (C)":          entry["reverse_tm"],
            "GC Content (%)":  entry.get("reverse_gc_content", 0.0),
            "Overlap Length":  entry["overlap_length"],
            "Overlap Tm":      entry["overlap_tm"],
            "Overlap GC (%)":  overlap_gc,
            "Mutations":       all_muts,
        },
    ]


def save_primer_csv(
    primer_data: List[dict],
    output_dir: Optional[Path] = None,
    filename: str = PRIMER_LIST_CSV,
) -> Path:
    """
    Write primer_list.csv from a list of PrimerSet dicts.
    Returns the path that was written.
    """
    out_dir = Path(output_dir) if output_dir else DEFAULT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / filename

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for entry in primer_data:
            for row in _primer_set_to_rows(entry):
                writer.writerow(row)

    return path


def load_primer_csv(path: Path) -> List[dict]:
    """
    Read a primer CSV file (primer_list.csv or primer_temp.csv).
    Returns a list of dicts with canonical CSV_COLUMNS keys.
    Missing columns get an empty-string default.
    Tries utf-8 first, falls back to latin-1 for files saved on older Windows.
    """
    for encoding in ("utf-8", "latin-1"):
        try:
            with open(path, newline="", encoding=encoding) as f:
                reader = csv.DictReader(f)
                rows = list(reader)
            # Normalise: ensure all canonical columns exist
            for row in rows:
                for col in CSV_COLUMNS:
                    row.setdefault(col, "")
            return rows
        except UnicodeDecodeError:
            continue
    raise ValueError(f"Cannot decode {path} as utf-8 or latin-1")


def save_primer_temp_csv(primer_rows: List[dict], output_dir: Path) -> Path:
    """
    Write unsaved editor changes to primer_temp.csv.
    primer_rows is already a flat list of per-primer dicts (from PrimerEditorWindow).
    """
    path = Path(output_dir) / PRIMER_TEMP_CSV
    fieldnames = list(primer_rows[0].keys()) if primer_rows else CSV_COLUMNS

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(primer_rows)

    return path


def promote_temp_to_main(output_dir: Path) -> Optional[Path]:
    """
    Copy primer_temp.csv → primer_list.csv (with backup), then delete temp.
    Returns the main path, or None if there was no temp file.
    """
    out_dir   = Path(output_dir)
    temp_path = out_dir / PRIMER_TEMP_CSV
    main_path = out_dir / PRIMER_LIST_CSV

    if not temp_path.exists():
        return None

    if main_path.exists():
        shutil.copy(main_path, out_dir / "primer_list.csv.backup")

    shutil.copy(temp_path, main_path)
    temp_path.unlink()
    return main_path


def active_primer_csv_path(output_dir: Path) -> Optional[Path]:
    """
    Return the path of the CSV that should be displayed:
    primer_temp.csv if it exists, otherwise primer_list.csv.
    Returns None if neither exists.
    """
    out_dir   = Path(output_dir)
    temp_path = out_dir / PRIMER_TEMP_CSV
    main_path = out_dir / PRIMER_LIST_CSV

    if temp_path.exists():
        return temp_path
    if main_path.exists():
        return main_path
    return None


# ===========================================================================
# Primer TXT (human-readable table)
# ===========================================================================

def save_primer_txt(
    primer_data: List[dict],
    output_dir: Optional[Path] = None,
    filename: str = PRIMER_LIST_TXT,
) -> Path:
    """
    Write primer_list.txt — the wide fixed-width human-readable table.
    Returns the path that was written.
    """
    out_dir = Path(output_dir) if output_dir else DEFAULT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / filename

    line_width = 165

    with open(path, "w", encoding="utf-8") as f:
        f.write("Site-Directed Mutagenesis Primer List\n")
        f.write("=" * line_width + "\n\n")

        header = (
            f"{'Primer Name':<20}{'Primer Sequence':<60}{'Len':>6}"
            f"{'Tm (C)':>9}{'GC%':>6}{'OverlapLen':>11}"
            f"{'OverlapTm':>10}{'OverlapGC%':>11}{'Notes':>30}\n"
        )
        f.write(header)
        f.write("-" * line_width + "\n")

        for entry in primer_data:
            base_name  = "_".join(entry["all_covered_mutations"])
            notes      = ",".join(entry["all_covered_mutations"])
            overlap_gc = entry.get("overlap_gc_content", 0.0)

            for direction in ("for", "rev"):
                is_fwd = direction == "for"
                seq    = entry["forward_primer"]    if is_fwd else entry["reverse_primer"]
                length = entry["forward_length"]    if is_fwd else entry["reverse_length"]
                tm     = entry["forward_tm"]        if is_fwd else entry["reverse_tm"]
                gc     = entry.get("forward_gc_content", 0.0) if is_fwd else entry.get("reverse_gc_content", 0.0)

                f.write(f"{base_name + '_' + direction:<20}")
                f.write(f"{seq:<60}")
                f.write(f"{length:>6}")
                f.write(f"{tm:>9.1f}")
                f.write(f"{gc:>6.1f}")
                f.write(f"{entry['overlap_length']:>11}")
                f.write(f"{entry['overlap_tm']:>10.1f}")
                f.write(f"{overlap_gc:>11.1f}")
                f.write(f"{notes:>30}\n")

            f.write("-" * line_width + "\n")

        f.write(f"\nSummary:\nTotal primer pairs: {len(primer_data)}\n")

    return path


# ===========================================================================
# Primer JSON  (dual-section format)
# ===========================================================================

def save_primer_json(
    primer_data: List[dict],
    output_dir: Optional[Path] = None,
    filename: str = PRIMER_LIST_JSON,
) -> Path:
    """
    Write primer_list.json with two sections:
      "primer_sets"  — original pair-level dicts (for protocol generation)
      "primers"      — flat per-primer rows matching CSV_COLUMNS (for display / tools)

    Both old callers (reading a plain list) and new callers (reading the dict)
    are handled by load_primer_json().
    """
    out_dir = Path(output_dir) if output_dir else DEFAULT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / filename

    flat_primers: List[dict] = []
    for entry in primer_data:
        for row in _primer_set_to_rows(entry):
            flat_primers.append(row)

    payload = {
        "primer_sets": primer_data,
        "primers":     flat_primers,
    }

    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    return path


def load_primer_json(path: Path) -> Dict:
    """
    Load primer_list.json.

    Returns a dict with keys:
      "primer_sets"  — list of pair-level dicts
      "primers"      — flat per-primer list

    Handles the old format (plain list of pair dicts) transparently.
    """
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)

    if isinstance(raw, dict):
        return {
            "primer_sets": raw.get("primer_sets", []),
            "primers":     raw.get("primers", []),
        }

    # Old format: plain list → reconstruct flat list on the fly
    flat: List[dict] = []
    for entry in raw:
        for row in _primer_set_to_rows(entry):
            flat.append(row)

    return {"primer_sets": raw, "primers": flat}


def save_primer_files(
    primer_data: List[dict],
    output_dir: Optional[Path] = None,
) -> Dict[str, Path]:
    """
    Convenience wrapper: write all three primer output formats at once.
    Returns a dict of {format: path}.
    """
    out_dir = Path(output_dir) if output_dir else DEFAULT_DIR
    paths = {
        "txt":  save_primer_txt(primer_data,  out_dir),
        "csv":  save_primer_csv(primer_data,  out_dir),
        "json": save_primer_json(primer_data, out_dir),
    }
    print("Primer files saved:")
    for fmt, p in paths.items():
        print(f"  [{fmt}] {p}")
    return paths


# ===========================================================================
# SLiCE primers CSV
# ===========================================================================

_SLICE_CSV_COLUMNS = [
    "Primer Name",
    "Sequence",
    "Length (bp)",
    "Tm (°C)",
    "Fragment",
    "Direction",
    "Overlap Length",
    "Overlap Tm",
]


def save_slice_primer_csv(primers: List[SLiCEPrimer], path: Path) -> Path:
    """Write SLiCE primers to a CSV file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=_SLICE_CSV_COLUMNS)
        writer.writeheader()
        for p in primers:
            writer.writerow({
                "Primer Name":   p.name,
                "Sequence":      p.sequence,
                "Length (bp)":   p.length,
                "Tm (°C)":       p.tm,
                "Fragment":      p.fragment_name,
                "Direction":     "Forward" if p.is_forward else "Reverse",
                "Overlap Length": p.overlap_length,
                "Overlap Tm":    p.overlap_tm,
            })

    return path


# ===========================================================================
# Variant databank JSON
# ===========================================================================

def load_databank(protocols_dir: Path) -> Dict[str, str]:
    """
    Load variant_databank.json from protocols_dir.
    Returns {} if the file does not exist or is malformed.
    """
    path = Path(protocols_dir) / DATABANK_JSON
    if not path.exists():
        return {}
    try:
        with open(path, "r", encoding="utf-8") as f:
            content = f.read().strip()
        return json.loads(content) if content else {}
    except json.JSONDecodeError:
        print(f"Warning: {path} is malformed — starting with empty databank.")
        return {}


def save_databank(databank: Dict[str, str], protocols_dir: Path) -> Path:
    """Write variant_databank.json."""
    path = Path(protocols_dir) / DATABANK_JSON
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(databank, f, indent=2)
    return path


def save_databank_undo(
    databank: Dict[str, str], protocols_dir: Path, undo_stack: list
) -> None:
    """
    Push current on-disk databank state onto undo_stack, then save new state.
    """
    path = Path(protocols_dir) / DATABANK_JSON
    if path.exists():
        try:
            with open(path, "r", encoding="utf-8") as f:
                current = json.load(f)
        except json.JSONDecodeError:
            current = {}
        undo_stack.append(copy.deepcopy(current))
    else:
        undo_stack.append({})

    save_databank(databank, protocols_dir)


def undo_databank(protocols_dir: Path, undo_stack: list) -> bool:
    """
    Pop the last saved databank state from undo_stack and restore it.
    Returns True on success, False if stack is empty.
    """
    if not undo_stack:
        return False
    last_state = undo_stack.pop()
    save_databank(last_state, protocols_dir)
    return True


# ===========================================================================
# Wildtype sequences JSON
# ===========================================================================

def wt_sequences_path(output_dir: Path) -> Path:
    return Path(output_dir) / WT_SEQUENCES_JSON


def load_wt_sequences(output_dir: Path) -> Dict[str, str]:
    """Load the wildtype_sequences.json store.  Returns {} if missing."""
    path = wt_sequences_path(output_dir)
    if not path.exists():
        return {}
    for encoding in ("utf-8", "latin-1"):
        try:
            with open(path, "r", encoding=encoding) as f:
                return json.load(f)
        except (UnicodeDecodeError, json.JSONDecodeError):
            continue
    print(f"Warning: Could not read {path}")
    return {}


def save_wt_sequences(sequences: Dict[str, str], output_dir: Path) -> Path:
    """Write the wildtype_sequences.json store."""
    path = wt_sequences_path(output_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(sequences, f, indent=2)
    return path


def save_wt_sequences_undo(
    sequences: Dict[str, str], output_dir: Path, undo_stack: list
) -> None:
    """Push current on-disk WT sequences onto undo_stack, then save new state."""
    path = wt_sequences_path(output_dir)
    if path.exists():
        try:
            with open(path, "r", encoding="utf-8") as f:
                current = json.load(f)
        except (json.JSONDecodeError, UnicodeDecodeError):
            current = {}
        undo_stack.append(copy.deepcopy(current))
    else:
        undo_stack.append({})
    save_wt_sequences(sequences, output_dir)


# ===========================================================================
# Protocol data JSON  (read-only from GUI perspective)
# ===========================================================================

def load_protocol_data(protocols_dir: Path) -> Optional[Dict]:
    """
    Load protocol_data.json written by MutagenesisProtocol.run().
    Returns None if the file does not exist.
    """
    path = Path(protocols_dir) / PROTOCOL_DATA_JSON
    if not path.exists():
        return None
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError) as exc:
        print(f"Warning: Could not read protocol_data.json: {exc}")
        return None


# ===========================================================================
# Cleanup helpers  (used by undo / reset operations in the GUI)
# ===========================================================================

def remove_primer_files(output_dir: Path) -> List[Path]:
    """
    Delete all three primer output files if they exist.
    Returns list of paths that were actually removed.
    """
    removed: List[Path] = []
    out_dir = Path(output_dir)
    for name in (PRIMER_LIST_TXT, PRIMER_LIST_CSV, PRIMER_LIST_JSON):
        p = out_dir / name
        if p.exists():
            p.unlink()
            removed.append(p)
    return removed


def remove_protocol_files(protocols_dir: Path) -> List[Path]:
    """
    Delete protocol output files (PDF, JSON, labels) if they exist.
    Returns list of paths that were actually removed.
    """
    removed: List[Path] = []
    pdir = Path(protocols_dir)
    for name in (
        "mutagenesis_protocol.pdf",
        PROTOCOL_DATA_JSON,
        "selected_variant_labels.pdf",
    ):
        p = pdir / name
        if p.exists():
            p.unlink()
            removed.append(p)
    return removed


# ===========================================================================
# FASTA parsing  (moved from MutationExtractorPage in GUI)
# ===========================================================================

def parse_fasta(path: Path) -> List[tuple]:
    """
    Parse a FASTA file and return a list of (header, sequence) tuples.
    Sequences are uppercased and concatenated across wrapped lines.
    Tries utf-8 first, falls back to latin-1 for files from older tools.
    """
    for encoding in ("utf-8", "latin-1"):
        try:
            seqs: List[tuple] = []
            with open(path, 'r', encoding=encoding) as f:
                header: Optional[str] = None
                parts:  List[str]     = []
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if header is not None:
                            seqs.append((header, ''.join(parts)))
                        header, parts = line[1:], []
                    elif line and header is not None:
                        parts.append(line.upper())
                if header is not None:
                    seqs.append((header, ''.join(parts)))
            return seqs
        except UnicodeDecodeError:
            continue
    raise ValueError(f"Cannot decode {path} as utf-8 or latin-1")
