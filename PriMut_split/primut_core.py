"""
primer_core.py
==============
Pure business logic for the Site-Directed Mutagenesis Primer Designer.
No GUI imports, no tkinter, no file I/O — only calculations and data structures.

Modules in this package
------------------------
primer_core.py   ← you are here  (logic / calculations)
primer_io.py     ← file I/O: CSV, JSON, PDF, TXT
primer_gui.py    ← all customtkinter / tkinter windows
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Try to import primer3.  On Windows the C extension sometimes loads but
# returns 0 or raises silently — we validate it once at import time with a
# known sequence.  _PRIMER3_OK is False if primer3 is missing or broken.
# The GUI reads this flag to show the user a clear warning banner.
# ---------------------------------------------------------------------------
_PRIMER3_OK = False
try:
    import primer3 as _primer3
    _test_tm = _primer3.calc_tm("ATGCGATCGATCGATCGATCGATCGATCGA")
    _PRIMER3_OK = (_test_tm is not None and _test_tm > 10.0)
except Exception:
    pass

if not _PRIMER3_OK:
    print(
        "WARNING: primer3 is not working correctly on this system.\n"
        "         Tm values will be reported as 0.0.\n"
        "         To fix: pip install primer3-py --force-reinstall"
    )


# ===========================================================================
# Constants
# ===========================================================================

GENETIC_CODE: Dict[str, str] = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

AMINO_ACID_MAP: Dict[str, str] = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
    'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu',
    'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
    'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
    '*': 'Ter', 'X': 'Xxx',
}




# ===========================================================================
# Utility functions (no GUI dependency)
# ===========================================================================

def natural_sort_key(text: str) -> list:
    """Natural (alphanumeric) sort key — sorts 'A9' before 'A10'."""
    def _atoi(t: str):
        return int(t) if t.isdigit() else t.lower()
    return [_atoi(c) for c in re.split(r'(\d+)', text)]


def translate_dna_to_protein(dna_sequence: str) -> str:
    """Translate a DNA sequence to a protein sequence (stops at first stop codon)."""
    protein: List[str] = []
    for i in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i + 3]
        protein.append(GENETIC_CODE.get(codon, 'X'))
    return ''.join(protein)


def get_three_letter_code(one_letter: str) -> str:
    """One-letter → three-letter amino acid code."""
    return AMINO_ACID_MAP.get(one_letter, 'Xxx')


def parse_mutation(mut_str: str) -> Tuple[str, int, str]:
    """Parse 'A123B' → ('A', 123, 'B').  Raises ValueError on bad format."""
    match = re.match(r'([A-Z])(\d+)([A-Z])', mut_str.strip())
    if not match:
        raise ValueError(f"Invalid mutation format: {mut_str!r}")
    return match.group(1), int(match.group(2)), match.group(3)


def mutation_position(mutation: str) -> int:
    """Extract numeric position from a mutation string like 'A123B'."""
    matches = re.findall(r'\d+', mutation)
    return int(matches[0]) if matches else 0


# ===========================================================================
# Tm calculation — nearest-neighbour with primer3 primary / NN fallback
# ===========================================================================

def calculate_tm(sequence: str) -> float:
    """
    Calculate Tm using primer3.

    Returns 0.0 if primer3 is not installed or not working correctly on this
    platform (see _PRIMER3_OK).  The GUI shows a warning banner in that case
    so the user knows the values are invalid and how to fix the installation.
    """
    if not _PRIMER3_OK or not sequence:
        return 0.0
    seq = sequence.upper().strip()
    if not seq:
        return 0.0
    try:
        return round(_primer3.calc_tm(seq), 2)
    except Exception:
        return 0.0


def calculate_gc_content(sequence: str) -> float:
    """GC content as a percentage (0–100), rounded to 1 decimal place."""
    if not sequence:
        return 0.0
    seq = sequence.upper()
    gc = seq.count('G') + seq.count('C')
    return round(gc / len(seq) * 100, 1)


def check_gc_content(sequence: str) -> Tuple[float, str, str]:
    """
    Returns (gc_percent, status, warning_message).
    status is one of: 'optimal', 'acceptable', 'warning'.
    """
    gc = calculate_gc_content(sequence)
    if 40 <= gc <= 60:
        return gc, 'optimal', ''
    elif 60 < gc <= 70:
        return gc, 'acceptable', 'GC content slightly high'
    elif gc > 70:
        return gc, 'warning', 'GC content too high (>70%)'
    else:
        return gc, 'warning', 'GC content too low (<40%)'


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence (handles upper/lower)."""
    comp = str.maketrans('ATGCatgcNn', 'TACGtacgNn')
    return sequence.translate(comp)[::-1].upper()


# ===========================================================================
# Data classes — shared between core, io, and GUI layers
# ===========================================================================

@dataclass
class SLiCEFragment:
    """A fragment to be amplified for SLiCE/Gibson assembly."""
    name: str
    start: int           # 0-based index into the full plasmid
    end: int             # 0-based index, exclusive
    sequence: str
    is_insert: bool = True  # True = insert, False = backbone


@dataclass
class SLiCEPrimer:
    """A primer designed for SLiCE/Gibson assembly."""
    name: str
    sequence: str
    tm: float
    length: int
    fragment_name: str
    is_forward: bool
    overlap_sequence: str = ""
    overlap_tm: float = 0.0
    overlap_length: int = 0


# PrimerSet is a plain dict in the original code — we keep that for now so
# existing callers don't need changing.  Type alias for documentation only.
PrimerSet = Dict


# ===========================================================================
# PrimerGenerator — SDM primer design logic
# ===========================================================================

class PrimerGenerator:
    """
    Site-directed mutagenesis primer generator.

    All methods are pure functions of sequence + mutation data.
    No file I/O, no GUI calls.
    """

    def __init__(self, flank_size_bases: int = 30) -> None:
        self.stop_codons = {"TAA", "TAG", "TGA"}

        # Reverse codon table: amino acid → list of coding triplets
        self.aa_to_codons: Dict[str, List[str]] = defaultdict(list)
        for codon, aa in GENETIC_CODE.items():
            self.aa_to_codons[aa].append(codon)

        self.flank_size_bases  = flank_size_bases
        self.flank_size_codons = flank_size_bases // 3
        self.MIN_PRIMER_LEN    = 30
        self.MIN_EXTENSION_LEN = 15
        self.TARGET_TM         = 68.0
        self.MIN_TM            = 55.0
        self.MAX_TM_DIFF       = 3.0
        self.MAX_OVERLAP_LEN   = 15

    # ------------------------------------------------------------------
    # Sequence utilities
    # ------------------------------------------------------------------

    @staticmethod
    def reverse_complement(sequence: str) -> str:
        return reverse_complement(sequence)

    @staticmethod
    def calculate_tm(sequence: str) -> float:
        return calculate_tm(sequence)

    @staticmethod
    def calculate_gc_content(sequence: str) -> float:
        return calculate_gc_content(sequence)

    @staticmethod
    def check_gc_content(sequence: str) -> Tuple[float, str, str]:
        return check_gc_content(sequence)

    # ------------------------------------------------------------------
    # Sequence validation
    # ------------------------------------------------------------------

    def check_sequence_validity(
        self, dna_sequence: str, flank_size: int = 30
    ) -> Tuple[bool, Optional[Tuple[str, str]]]:
        """
        Validate DNA sequence structure and reading frame.
        Returns (is_valid, (start_codon, first_aa)) or (False, None).
        """
        if len(dna_sequence) % 3 != 0:
            print(f"Error: Sequence length {len(dna_sequence)} is not a multiple of 3.")
            return False, None

        if len(dna_sequence) < (flank_size * 2 + 6):
            print("Error: Sequence too short for flanking regions.")
            return False, None

        start_codon = dna_sequence[flank_size:flank_size + 3]
        first_aa = GENETIC_CODE.get(start_codon, 'X')

        stop_codon_pos = dna_sequence[-(flank_size + 3):-flank_size]
        if stop_codon_pos not in self.stop_codons:
            print(
                f"Error: No valid stop codon before last {flank_size} bases, "
                f"found {stop_codon_pos!r}."
            )
            return False, None

        print(f"Sequence check passed. First codon: {start_codon} ({first_aa})")
        return True, (start_codon, first_aa)

    def validate_mutations_against_sequence(
        self, dna_sequence: str, variant_list: List[List[str]]
    ) -> None:
        """
        Confirm every mutation matches the wildtype residue in dna_sequence.
        Raises ValueError on mismatch.
        """
        for variant in variant_list:
            for mut in variant:
                original_aa, pos, _ = parse_mutation(mut)
                adj_pos_codons = pos + self.flank_size_codons
                codon_start = (adj_pos_codons - 1) * 3
                codon_seq = dna_sequence[codon_start:codon_start + 3]
                wt_aa = GENETIC_CODE.get(codon_seq, 'X')
                if wt_aa != original_aa:
                    raise ValueError(
                        f"Mutation {mut} does not match WT residue {wt_aa!r} "
                        f"at coding AA position {pos} (codon {codon_seq})"
                    )

    # ------------------------------------------------------------------
    # Codon helpers
    # ------------------------------------------------------------------

    def get_best_codon(self, amino_acid: str, original_codon: str) -> str:
        """Return the original codon if it codes for amino_acid, else the first available."""
        codons = self.aa_to_codons.get(amino_acid, ["NNN"])
        return original_codon if original_codon in codons else codons[0]

    # ------------------------------------------------------------------
    # Mutation grouping
    # ------------------------------------------------------------------

    @staticmethod
    def flatten_mutations(mutation_groups) -> list:
        """Accept mixed list of strings / tuples and return a flat list."""
        flat = []
        for group in mutation_groups:
            if (
                isinstance(group, (list, tuple))
                and all(isinstance(i, tuple) and len(i) == 3 for i in group)
            ):
                flat.extend(group)
            else:
                flat.append(group)
        return flat

    def merge_close_mutations(
        self, mutations: list, max_distance_nt: int = 21
    ) -> List[list]:
        """Group mutations that are within max_distance_nt of each other."""
        sorted_muts = sorted(mutations, key=lambda x: x[1])
        if not sorted_muts:
            return []

        merged: List[list] = []
        current_group = [sorted_muts[0]]
        for mut in sorted_muts[1:]:
            if (mut[1] - current_group[-1][1]) * 3 <= max_distance_nt:
                current_group.append(mut)
            else:
                merged.append(current_group)
                current_group = [mut]
        merged.append(current_group)
        return merged

    def count_nt_changes(self, wt_seq: str, mutations: list) -> int:
        """Count nucleotide changes introduced by mutations relative to wt_seq."""
        mutated = self._create_mutated_sequence(wt_seq, mutations)
        return sum(1 for a, b in zip(wt_seq, mutated) if a != b)

    def split_mutations_by_max_nt_changes(
        self, wt_seq: str, mutation_group: list, max_nt_changes: int = 6
    ) -> List[list]:
        """Split a mutation group so each sub-group ≤ max_nt_changes."""
        sorted_muts = sorted(mutation_group, key=lambda m: m[1])
        split_groups: List[list] = []
        current: list = []

        for mut in sorted_muts:
            test = current + [mut]
            if self.count_nt_changes(wt_seq, test) > max_nt_changes:
                if current:
                    split_groups.append(current)
                current = [mut]
            else:
                current = test

        if current:
            split_groups.append(current)
        return split_groups

    # ------------------------------------------------------------------
    # Sequence manipulation
    # ------------------------------------------------------------------

    def _create_mutated_sequence(self, dna_sequence: str, mutations: list) -> str:
        mutated = list(dna_sequence)
        for aa, pos, new_aa in mutations:
            idx = ((pos - 1) + self.flank_size_codons) * 3
            original_codon = dna_sequence[idx:idx + 3]
            new_codon = self.get_best_codon(new_aa, original_codon)
            if GENETIC_CODE.get(original_codon, 'X') != aa:
                print(f"Warning: Original amino acid mismatch at position {pos}")
            mutated[idx:idx + 3] = list(new_codon)
        return ''.join(mutated)

    # ------------------------------------------------------------------
    # Overlap optimisation
    # ------------------------------------------------------------------

    def _extension_mutation_span(
        self, mutation_positions: List[int], overlap_start: int, overlap_end: int
    ) -> int:
        outside = [p for p in mutation_positions if p < overlap_start or p >= overlap_end]
        return (max(outside) - min(outside)) if outside else 0

    def _find_optimal_overlap(
        self,
        mutated_seq: str,
        codon_starts: List[int],
        codon_ends: List[int],
    ) -> Tuple[int, int, float]:
        """Return (overlap_start, overlap_end, overlap_tm)."""
        min_inside_flank       = 4
        crit_flank_limit       = 2
        reverse_flank_weight   = 2.0
        forward_flank_weight   = 1.5
        wt_inside_weight       = 0.5
        penalty_per_missing    = 50
        upstream_bias_weight   = 0.2

        mutation_positions = sorted(set(codon_starts))
        seq_len = len(mutated_seq)

        def search_len(curr_len):
            best_score = -1e9
            best = (None, None, None, (None, None), best_score)
            for start_pos in range(0, seq_len - curr_len + 1):
                end_pos = start_pos + curr_len
                muts_inside = [m for m in mutation_positions if start_pos <= m < end_pos]
                if not muts_inside:
                    continue
                first_inside = min(muts_inside)
                last_inside  = max(muts_inside)
                rev_flank    = first_inside - start_pos
                fwd_flank    = end_pos - (last_inside + 3)
                if rev_flank < crit_flank_limit or fwd_flank < crit_flank_limit:
                    continue

                count        = len(muts_inside)
                ext_span     = self._extension_mutation_span(mutation_positions, start_pos, end_pos)
                overlap_seq  = mutated_seq[start_pos:end_pos]
                overlap_tm   = calculate_tm(overlap_seq)
                total_mut_nt = count * 3
                total_wt_nt  = curr_len - total_mut_nt

                flank_penalty = 0
                if fwd_flank < min_inside_flank:
                    flank_penalty += (min_inside_flank - fwd_flank) * penalty_per_missing
                if rev_flank < min_inside_flank:
                    flank_penalty += (min_inside_flank - rev_flank) * penalty_per_missing

                score = (
                    count * 1000
                    + fwd_flank * forward_flank_weight
                    + rev_flank * reverse_flank_weight
                    + total_wt_nt * wt_inside_weight
                    - ext_span * 10
                    - flank_penalty
                    - abs(overlap_tm - self.TARGET_TM)
                    - start_pos * upstream_bias_weight
                )
                if score > best[4]:
                    best = (start_pos, end_pos, overlap_tm, (fwd_flank, rev_flank), score)
            return best

        # Stage 1: exact MAX_OVERLAP_LEN
        s, e, tm, (ff, rf), _ = search_len(self.MAX_OVERLAP_LEN)
        if s is not None and ff >= min_inside_flank and rf >= min_inside_flank:
            return s, e, tm

        # Stage 2: allow +1, +2
        best_candidate = None
        best_score = -1e9
        for extra in (1, 2):
            cs, ce, ctm, _, cscore = search_len(self.MAX_OVERLAP_LEN + extra)
            if cs is not None and cscore > best_score:
                best_candidate = (cs, ce, ctm)
                best_score = cscore
        if best_candidate is not None:
            return best_candidate

        # Fallback: centre between first/last mutation
        first_mut  = min(mutation_positions)
        last_mut   = max(mutation_positions)
        centre     = (first_mut + last_mut) // 2
        best_start = max(0, centre - self.MAX_OVERLAP_LEN // 2)
        best_end   = min(seq_len, best_start + self.MAX_OVERLAP_LEN)
        return best_start, best_end, calculate_tm(mutated_seq[best_start:best_end])

    # ------------------------------------------------------------------
    # Primer extension
    # ------------------------------------------------------------------

    def _extend_primer_one_base(
        self,
        mutated_seq: str,
        current_primer: str,
        overlap_seq: str,
        extension_start: int,
        extension_end: int,
        is_reverse: bool,
    ) -> Tuple[str, float, int, int]:
        MAX_PRIMER_LEN = 45
        if len(current_primer) >= MAX_PRIMER_LEN:
            return current_primer, calculate_tm(current_primer), extension_start, extension_end

        if is_reverse:
            if extension_start == 0:
                return current_primer, calculate_tm(current_primer), extension_start, extension_end
            extension_start = max(0, extension_start - 1)
            ext_template = mutated_seq[extension_start:extension_end]
            extension    = reverse_complement(ext_template)
            primer       = reverse_complement(overlap_seq) + extension
        else:
            if extension_end >= len(mutated_seq):
                return current_primer, calculate_tm(current_primer), extension_start, extension_end
            extension_end = min(len(mutated_seq), extension_end + 1)
            extension     = mutated_seq[extension_start:extension_end]
            primer        = overlap_seq + extension

        tm = calculate_tm(primer)
        return primer, tm, extension_start, extension_end

    # ------------------------------------------------------------------
    # Main entry point — generate primer sets for a construct
    # ------------------------------------------------------------------

    def generate_primers_for_construct(
        self, dna_sequence: str, mutations: list
    ) -> List[PrimerSet]:
        """
        Design forward and reverse primers for all mutations on dna_sequence.
        Returns a list of PrimerSet dicts (one per mutation sub-group).
        """
        mutations = self.flatten_mutations(mutations)
        primer_sets: List[PrimerSet] = []
        mutation_groups = self.merge_close_mutations(mutations)

        MAX_PRIMER_LEN = 45
        MAX_MUT_NT     = 6

        for group in mutation_groups:
            split_groups   = self.split_mutations_by_max_nt_changes(dna_sequence, group, MAX_MUT_NT)
            previous_group: list = []

            for sub_group in split_groups:
                codon_starts = [((m[1] - 1) + self.flank_size_codons) * 3 for m in sub_group]
                codon_ends   = [s + 3 for s in codon_starts]

                all_applied  = previous_group + sub_group
                mutated_seq  = self._create_mutated_sequence(dna_sequence, all_applied)

                overlap_start, overlap_end, _ = self._find_optimal_overlap(
                    mutated_seq, codon_starts, codon_ends
                )
                overlap_seq = mutated_seq[overlap_start:overlap_end]

                # --- Forward primer ---
                fwd_start  = overlap_end
                fwd_end    = overlap_end + self.MIN_EXTENSION_LEN
                fwd_primer = overlap_seq + mutated_seq[fwd_start:fwd_end]
                fwd_tm     = calculate_tm(fwd_primer)

                # --- Reverse primer ---
                rev_end      = overlap_start
                rev_start    = max(0, rev_end - self.MIN_EXTENSION_LEN)
                rev_template = mutated_seq[rev_start:rev_end]
                rev_primer   = reverse_complement(overlap_seq) + reverse_complement(rev_template)
                rev_tm       = calculate_tm(rev_primer)

                # Balance Tm
                tm_diff  = abs(fwd_tm - rev_tm)
                attempts = 0
                while (
                    (tm_diff > self.MAX_TM_DIFF or fwd_tm < self.TARGET_TM or rev_tm < self.TARGET_TM)
                    and attempts < 30
                ):
                    if fwd_tm < rev_tm and len(fwd_primer) < MAX_PRIMER_LEN:
                        fwd_primer, fwd_tm, fwd_start, fwd_end = self._extend_primer_one_base(
                            mutated_seq, fwd_primer, overlap_seq, fwd_start, fwd_end, is_reverse=False
                        )
                    elif len(rev_primer) < MAX_PRIMER_LEN:
                        rev_primer, rev_tm, rev_start, rev_end = self._extend_primer_one_base(
                            mutated_seq, rev_primer, overlap_seq, rev_start, rev_end, is_reverse=True
                        )
                    else:
                        break
                    tm_diff = abs(fwd_tm - rev_tm)
                    attempts += 1

                all_mut_set  = {(m[0], m[1], m[2]) for m in all_applied}
                all_mut_list = sorted(all_mut_set, key=lambda x: x[1])
                set_id       = "_".join(f"{m[0]}{m[1]}{m[2]}" for m in all_mut_list)
                depends_on   = [] if not primer_sets else [primer_sets[-1]['set_id']]

                fwd_gc,  fwd_gc_status,  _  = check_gc_content(fwd_primer)
                rev_gc,  rev_gc_status,  _  = check_gc_content(rev_primer)
                ovlp_gc, ovlp_gc_status, _  = check_gc_content(overlap_seq)

                primer_sets.append({
                    'set_id':               set_id,
                    'depends_on':           depends_on,
                    'mutations':            [f"{m[0]}{m[1]}{m[2]}" for m in sub_group],
                    'all_covered_mutations':[f"{m[0]}{m[1]}{m[2]}" for m in all_mut_list],
                    'overlap_sequence':     overlap_seq,
                    'overlap_tm':           calculate_tm(overlap_seq),
                    'forward_primer':       fwd_primer,
                    'reverse_primer':       rev_primer,
                    'forward_tm':           fwd_tm,
                    'reverse_tm':           rev_tm,
                    'forward_length':       len(fwd_primer),
                    'reverse_length':       len(rev_primer),
                    'overlap_length':       len(overlap_seq),
                    'tm_difference':        tm_diff,
                    'forward_gc_content':   fwd_gc,
                    'forward_gc_status':    fwd_gc_status,
                    'reverse_gc_content':   rev_gc,
                    'reverse_gc_status':    rev_gc_status,
                    'overlap_gc_content':   ovlp_gc,
                    'overlap_gc_status':    ovlp_gc_status,
                })

                previous_group = list(all_applied)

        return primer_sets


# ===========================================================================
# SLiCEPrimerDesigner — assembly primer design logic
# ===========================================================================

class SLiCEPrimerDesigner:
    """
    Design primers for SLiCE/Gibson assembly with configurable parameters.
    No GUI, no file I/O.
    """

    def __init__(
        self,
        overlap_length_range: Tuple[int, int] = (15, 25),
        overlap_tm_range:     Tuple[int, int] = (55, 65),
        primer_tm_target:     float           = 60.0,
        primer_tm_range:      Tuple[float, float] = (58.0, 62.0),
    ) -> None:
        self.overlap_min_len  = overlap_length_range[0]
        self.overlap_max_len  = overlap_length_range[1]
        self.overlap_min_tm   = overlap_tm_range[0]
        self.overlap_max_tm   = overlap_tm_range[1]
        self.primer_tm_target = primer_tm_target
        self.primer_tm_min    = primer_tm_range[0]
        self.primer_tm_max    = primer_tm_range[1]

    # Delegate to module-level functions so there is exactly one implementation
    @staticmethod
    def reverse_complement(sequence: str) -> str:
        return reverse_complement(sequence)

    @staticmethod
    def calculate_tm(sequence: str) -> float:
        return calculate_tm(sequence)

    def find_optimal_overlap(
        self, seq1: str, seq2: str, is_circular: bool = False
    ) -> Tuple[str, int, float]:
        """Find best overlap between the 3'-end of seq1 and the 5'-end of seq2."""
        for length in range(self.overlap_max_len, self.overlap_min_len - 1, -1):
            if len(seq1) < length or len(seq2) < length:
                continue
            overlap = seq1[-length:]
            if seq2[:length].upper() == overlap.upper():
                tm = calculate_tm(overlap)
                if self.overlap_min_tm <= tm <= self.overlap_max_tm:
                    return overlap, length, tm

        # Fallback: return shortest acceptable overlap
        length  = min(self.overlap_min_len, len(seq1), len(seq2))
        overlap = seq1[-length:]
        return overlap, length, calculate_tm(overlap)

    def extend_to_tm(
        self,
        template: str,
        start_pos: int,
        direction: str,
        target_tm: float,
        max_length: int = 30,
    ) -> Tuple[str, float]:
        """Extend a primer from start_pos until it reaches primer_tm_min."""
        current_seq = ""
        if direction == 'forward':
            for length in range(15, min(max_length, len(template) - start_pos)):
                current_seq = template[start_pos:start_pos + length]
                if calculate_tm(current_seq) >= self.primer_tm_min:
                    return current_seq, calculate_tm(current_seq)
        else:
            end_pos = start_pos
            for length in range(15, min(max_length, end_pos)):
                current_seq = template[end_pos - length:end_pos]
                if calculate_tm(current_seq) >= self.primer_tm_min:
                    return current_seq, calculate_tm(current_seq)

        return current_seq, calculate_tm(current_seq)

    def design_primers_for_fragments(
        self, fragments: List[SLiCEFragment], circular: bool = True
    ) -> List[SLiCEPrimer]:
        """Design forward and reverse primers for every fragment in the assembly."""
        primers: List[SLiCEPrimer] = []

        for i, fragment in enumerate(fragments):
            n = len(fragments)
            prev_frag = fragments[i - 1] if i > 0 else (fragments[-1] if circular else None)
            next_frag = fragments[i + 1] if i < n - 1 else (fragments[0] if circular else None)

            # Forward primer
            if prev_frag:
                ovlp_seq, ovlp_len, ovlp_tm = self.find_optimal_overlap(
                    prev_frag.sequence[-self.overlap_max_len:],
                    fragment.sequence[:self.overlap_max_len],
                )
            else:
                ovlp_seq, ovlp_len, ovlp_tm = "", 0, 0.0

            binding_seq, _ = self.extend_to_tm(fragment.sequence, 0, 'forward', self.primer_tm_target)
            fwd_seq = ovlp_seq + binding_seq

            primers.append(SLiCEPrimer(
                name=f"{fragment.name}_F",
                sequence=fwd_seq,
                tm=calculate_tm(fwd_seq),
                length=len(fwd_seq),
                fragment_name=fragment.name,
                is_forward=True,
                overlap_sequence=ovlp_seq,
                overlap_tm=ovlp_tm,
                overlap_length=ovlp_len,
            ))

            # Reverse primer
            if next_frag:
                ovlp_seq, ovlp_len, ovlp_tm = self.find_optimal_overlap(
                    fragment.sequence[-self.overlap_max_len:],
                    next_frag.sequence[:self.overlap_max_len],
                )
                ovlp_seq_rc = reverse_complement(ovlp_seq)
            else:
                ovlp_seq_rc, ovlp_len, ovlp_tm = "", 0, 0.0

            binding_seq, _ = self.extend_to_tm(
                fragment.sequence, len(fragment.sequence), 'reverse', self.primer_tm_target
            )
            binding_rc = reverse_complement(binding_seq)
            rev_seq    = ovlp_seq_rc + binding_rc

            primers.append(SLiCEPrimer(
                name=f"{fragment.name}_R",
                sequence=rev_seq,
                tm=calculate_tm(rev_seq),
                length=len(rev_seq),
                fragment_name=fragment.name,
                is_forward=False,
                overlap_sequence=ovlp_seq_rc,
                overlap_tm=ovlp_tm,
                overlap_length=ovlp_len,
            ))

        return primers

    def validate_assembly(
        self, fragments: List[SLiCEFragment], circular: bool = True
    ) -> Dict:
        """
        Validate that all fragments can be assembled.
        Returns {'valid': bool, 'issues': [str], 'total_fragments': int, 'total_length': int}.
        """
        issues: List[str] = []
        n = len(fragments)

        for frag in fragments:
            if len(frag.sequence) < 50:
                issues.append(f"Fragment {frag.name} is very short ({len(frag.sequence)} bp)")

        for i in range(n):
            next_i = (i + 1) % n if circular else i + 1
            if next_i >= n:
                continue
            f1, f2 = fragments[i], fragments[next_i]
            _, length, tm = self.find_optimal_overlap(
                f1.sequence[-self.overlap_max_len:],
                f2.sequence[:self.overlap_max_len],
            )
            if length < self.overlap_min_len:
                issues.append(f"Insufficient overlap between {f1.name} and {f2.name}")
            elif tm < self.overlap_min_tm:
                issues.append(f"Low overlap Tm between {f1.name} and {f2.name}: {tm}°C")
            elif tm > self.overlap_max_tm:
                issues.append(f"High overlap Tm between {f1.name} and {f2.name}: {tm}°C")

        return {
            'valid':            len(issues) == 0,
            'issues':           issues,
            'total_fragments':  n,
            'total_length':     sum(len(f.sequence) for f in fragments),
        }


# ===========================================================================
# PlasmidAnnotator — sequence feature detection (pure logic)
# ===========================================================================

class PlasmidAnnotator:
    """Detect common plasmid features by regex pattern matching."""

    FEATURES: Dict[str, List[Tuple[str, str]]] = {
        'promoter': [
            (r'TTGACA.{15,19}TATAAT', 'bacterial promoter (-35/-10)'),
            (r'TATA[AT]A[AT]',         'TATA box'),
        ],
        'terminator': [
            (r'AAAAAA', 'poly-A signal'),
            (r'TTTTTT', 'transcription terminator'),
        ],
        'restriction_site': [
            (r'GAATTC', 'EcoRI'),  (r'GGATCC', 'BamHI'),
            (r'AAGCTT', 'HindIII'),(r'CTGCAG', 'PstI'),
            (r'GTCGAC', 'SalI'),   (r'GCTAGC', 'NheI'),
            (r'CATATG', 'NdeI'),   (r'GGTCTC', 'BsaI'),
            (r'CACCTGC', 'AarI'),
        ],
        'ori': [
            (r'TTGAGA[GT]ACAGC', 'ColE1 ori'),
            (r'AATGATACGGCGAC',  'pBR322 ori'),
        ],
        'marker': [
            (r'ATG[ATGC]{300,900}(TAG|TAA|TGA)', 'resistance gene'),
        ],
    }

    @staticmethod
    def annotate_sequence(sequence: str) -> List[Dict]:
        """Return list of feature dicts sorted by position."""
        annotations: List[Dict] = []
        seq_upper = sequence.upper()
        for feature_type, patterns in PlasmidAnnotator.FEATURES.items():
            for pattern, description in patterns:
                for match in re.finditer(pattern, seq_upper):
                    annotations.append({
                        'type':        feature_type,
                        'start':       match.start(),
                        'end':         match.end(),
                        'description': description,
                        'sequence':    match.group(),
                    })
        annotations.sort(key=lambda x: x['start'])
        return annotations


# ===========================================================================
# MutagenesisProtocol — protocol planning logic (no GUI, but does file I/O
# via reportlab / json).  File-writing helpers are kept here so callers
# don't need to import primer_io for the protocol run; the pure data
# transformations (get_base_variant, variant_key, etc.) live here.
# ===========================================================================

class MutagenesisProtocol:
    """
    Plan a mutagenesis protocol — which variants to make in which order.

    Heavy on logic, light on I/O.  The generate_pdf / save_protocol_json
    methods write files but do not interact with the GUI.
    """

    def __init__(
        self,
        variant_input,
        max_mutations_per_step: int,
        variant_prefix: str,
        output_dir_path: str,
        start_col: int = 2,
        start_row: int = 16,
        list_existing_as_steps: bool = False,
        undo_stack: Optional[list] = None,
    ) -> None:
        from pathlib import Path
        import json, copy
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.lib import colors
        from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer,
                                        Table, TableStyle, PageBreak)

        self._Path = Path
        self._json = json
        self._copy = copy
        self._A4   = A4
        self._getSampleStyleSheet = getSampleStyleSheet
        self._ParagraphStyle      = ParagraphStyle
        self._colors              = colors
        self._SimpleDocTemplate   = SimpleDocTemplate
        self._Paragraph           = Paragraph
        self._Spacer              = Spacer
        self._Table               = Table
        self._TableStyle          = TableStyle
        self._PageBreak           = PageBreak

        self.variant_input          = variant_input
        self.max_mutations_per_step = max_mutations_per_step
        self.variant_prefix         = variant_prefix
        self.start_col              = start_col
        self.start_row              = start_row
        self.mutagenesis_dir        = Path(output_dir_path)
        self.primer_json_path       = self.mutagenesis_dir / "primer_list.json"
        self.output_dir             = self.mutagenesis_dir / "protocols"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.pdf_path               = self.output_dir / "mutagenesis_protocol.pdf"
        self.databank_file          = self.output_dir / "variant_databank.json"
        self.list_existing_as_steps = list_existing_as_steps
        self.undo_stack             = undo_stack if undo_stack is not None else []

    # ------------------------------------------------------------------
    # Static helpers
    # ------------------------------------------------------------------

    @staticmethod
    def variant_key(mutations) -> str:
        return ','.join(sorted(mutations, key=mutation_position))

    @staticmethod
    def extract_variant_number(name: str) -> int:
        match = re.search(r'variant (\d+)', name)
        return int(match.group(1)) if match else -1

    @staticmethod
    def normalize_mut(m: str) -> str:
        return m.strip().upper()

    def get_base_variant(self, current, existing_variants: dict):
        best = ()
        for var in existing_variants:
            if set(var).issubset(current) and len(var) > len(best):
                best = var
        return best

    # ------------------------------------------------------------------
    # Protocol generation
    # ------------------------------------------------------------------

    def save_protocol_json(
        self, protocol_by_round, final_variants, existing_input_variants,
        variant_to_label_map, variant_to_final_variant, used_variants
    ) -> None:
        protocol_data: Dict = {
            "existing_input_variants": existing_input_variants,
            "final_variants_overview": {},
            "protocol_steps":          {},
            "all_variants":            {},
        }

        for label in sorted(variant_to_label_map, key=natural_sort_key):
            sorted_mutations = sorted(variant_to_label_map[label], key=mutation_position)
            protocol_data["final_variants_overview"][label] = {
                "final_variant": variant_to_final_variant[label],
                "mutations":     sorted_mutations,
            }

        for round_num in sorted(protocol_by_round):
            steps = sorted(protocol_by_round[round_num],
                           key=lambda r: self.extract_variant_number(r[1]))
            step_data = []
            for row in steps:
                sorted_muts = sorted(row[2].split(', '), key=mutation_position)
                step_data.append({
                    "new_variant":         row[0],
                    "parent_variant":      row[1],
                    "mutations_added":     sorted_muts,
                    "mutations_added_str": ', '.join(sorted_muts),
                })
            protocol_data["protocol_steps"][round_num] = step_data

        for muts, name in used_variants.items():
            protocol_data["all_variants"][name] = ', '.join(muts) if muts else '(none)'

        out_path = self.output_dir / "protocol_data.json"
        with open(out_path, "w") as f:
            self._json.dump(protocol_data, f, indent=2)
        print(f"Protocol data saved: {out_path}")

    def generate_pdf(
        self, protocol_by_round, final_variants, filename: str,
        existing_input_variants, variant_to_label_map,
        variant_to_final_variant, used_variants
    ) -> None:
        doc      = self._SimpleDocTemplate(filename, pagesize=self._A4)
        elements = []
        styles   = self._getSampleStyleSheet()
        cell_style = self._ParagraphStyle('mut_cell', fontSize=9, leading=11)
        P  = self._Paragraph
        T  = self._Table
        TS = self._TableStyle
        c  = self._colors
        PB = self._PageBreak
        Sp = self._Spacer

        if existing_input_variants:
            elements.append(P("<b>Pre-existing Variants in Databank</b>", styles['Heading2']))
            elements.append(Sp(1, 12))
            data = [['Input Label', 'Existing Variant Label']]
            for label, existing_label in sorted(existing_input_variants.items()):
                data.append([label, existing_label or '(unknown)'])
            tbl = T(data, colWidths=[150, 250])
            tbl.setStyle(TS([
                ('BACKGROUND', (0, 0), (-1, 0), c.lightgrey),
                ('GRID',       (0, 0), (-1, -1), 0.5, c.grey),
            ]))
            elements += [tbl, PB()]

        if not protocol_by_round:
            elements.append(P("<b>No new variants generated.</b>", styles['Heading2']))
            elements.append(Sp(1, 12))
            msg = ("All input variants already exist in the databank."
                   if existing_input_variants
                   else "No variants were generated because no input variants were provided.")
            elements.append(P(msg, styles['Normal']))
            doc.build(elements)
            return

        elements.append(P("<b>Final Variants Overview</b>", styles['Heading2']))
        final_table_data = [['Final Variant', 'Mutations']]
        for label in sorted(variant_to_label_map, key=natural_sort_key):
            sorted_muts = sorted(variant_to_label_map[label], key=mutation_position)
            final_table_data.append([
                variant_to_final_variant[label],
                P(', '.join(sorted_muts), cell_style),
            ])
        tbl = T(final_table_data, colWidths=[150, 310], repeatRows=1)
        tbl.setStyle(TS([
            ('BACKGROUND', (0, 0), (-1, 0), c.lightgrey),
            ('GRID',       (0, 0), (-1, -1), 0.5, c.grey),
        ]))
        elements += [tbl, Sp(1, 12)]

        elements.append(P("<b>Step-by-Step Mutagenesis Protocol (Grouped by Round)</b>",
                          styles['Heading2']))
        for i, round_num in enumerate(sorted(protocol_by_round), start=1):
            steps = sorted(protocol_by_round[round_num],
                           key=lambda r: self.extract_variant_number(r[1]))
            elements.append(P(f"<b>Step {i}</b>", styles['Heading3']))
            table_data = [['New Variant', 'Parent Variant', 'Mutations Added']]
            for row in steps:
                sorted_muts = ', '.join(sorted(row[2].split(', '), key=mutation_position))
                table_data.append([row[0], row[1], P(sorted_muts, cell_style)])
            step_tbl = T(table_data, colWidths=[100, 120, 160], repeatRows=1)
            step_tbl.setStyle(TS([
                ('BACKGROUND', (0, 0), (-1, 0), c.lightgrey),
                ('GRID',       (0, 0), (-1, -1), 0.5, c.grey),
            ]))
            elements += [step_tbl, Sp(1, 10)]

        elements.append(PB())
        elements.append(P("<b>All Produced Variants</b>", styles['Heading2']))
        all_data = [['Variant Name', 'Mutations']]
        for muts, name in used_variants.items():
            mut_text = ', '.join(muts) if muts else '(none)'
            all_data.append([name, P(mut_text, cell_style)])
        tbl = T(all_data, colWidths=[140, 320], repeatRows=1)
        tbl.setStyle(TS([
            ('BACKGROUND',  (0, 0), (-1, 0), c.lightgrey),
            ('GRID',        (0, 0), (-1, -1), 0.5, c.grey),
            ('VALIGN',      (0, 0), (-1, -1), 'TOP'),
            ('LEFTPADDING', (0, 0), (-1, -1), 5),
            ('RIGHTPADDING',(0, 0), (-1, -1), 5),
        ]))
        elements.append(tbl)
        doc.build(elements)

    def run(self) -> None:
        """Execute the full protocol planning and write all output files."""
        import copy
        from collections import Counter, defaultdict

        if not self.primer_json_path.exists():
            raise FileNotFoundError(f"Primer JSON not found: {self.primer_json_path}")

        with open(self.primer_json_path, "r") as f:
            raw = self._json.load(f)

        # Support both old format (plain list) and new format (dict with primer_sets)
        primer_data = raw.get("primer_sets", []) if isinstance(raw, dict) else raw

        mutation_to_groups: Dict[str, List[int]] = {}
        primer_groups: Dict[int, list] = {}
        for idx, primer_set in enumerate(primer_data):
            group_muts = [self.normalize_mut(m) for m in primer_set["mutations"]]
            primer_groups[idx] = sorted(group_muts, key=mutation_position)
            for mut in group_muts:
                mutation_to_groups.setdefault(mut, []).append(idx)

        if self.list_existing_as_steps:
            variant_databank  = {}
            existing_variants = {(): 'wildtype'}
        else:
            if self.databank_file.exists():
                try:
                    with open(self.databank_file, "r") as f:
                        content = f.read().strip()
                        variant_databank = self._json.loads(content) if content else {}
                except self._json.JSONDecodeError:
                    print("Warning: Databank file invalid. Starting fresh.")
                    variant_databank = {}
            else:
                variant_databank = {}

            existing_variants: Dict = {(): 'wildtype'}
            for key, mutation_str in variant_databank.items():
                muts = tuple(mutation_str.split(',')) if mutation_str != '(none)' else ()
                existing_variants[muts] = key

        all_mutations  = [mut for v in self.variant_input for mut in v]
        mutation_freq  = Counter(all_mutations)
        variants       = [sorted(v, key=lambda x: -mutation_freq[x]) for v in self.variant_input]

        final_variants:            defaultdict = defaultdict(list)
        variant_to_label_map:      Dict        = {}
        variant_to_final_variant:  Dict        = {}
        protocol_by_round:         defaultdict = defaultdict(list)
        used_variants:             Dict        = {}
        existing_input_variants:   Dict        = {}

        existing_ids = [
            int(v_id[len(self.variant_prefix):])
            for v_id in variant_databank
            if v_id.startswith(self.variant_prefix)
            and v_id[len(self.variant_prefix):].isdigit()
        ]
        next_id_num = max(existing_ids, default=0) + 1

        for muts, name in existing_variants.items():
            mutation_str = ','.join(sorted(muts, key=mutation_position))
            if mutation_str not in variant_databank.values():
                new_id = f"{self.variant_prefix}{next_id_num:02d}"
                variant_databank[new_id] = mutation_str
                next_id_num += 1

        planning_variants = {(): 'wildtype'} if self.list_existing_as_steps else dict(existing_variants)

        for idx, target_mutations in enumerate(variants):
            target           = tuple(target_mutations)
            variant_key_str  = self.variant_key(target)
            variant_label    = f"{self.variant_prefix}{next_id_num:02d}"

            if variant_key_str in variant_databank.values():
                if not self.list_existing_as_steps:
                    for muts, name in planning_variants.items():
                        if variant_key_str == self.variant_key(muts):
                            existing_input_variants[variant_label] = name
                            variant_to_label_map[variant_label]    = target
                            variant_to_final_variant[variant_label] = name
                            used_variants[tuple(sorted(
                                map(self.normalize_mut, target), key=mutation_position
                            ))] = name
                            break
                    continue

            base_variant      = self.get_base_variant(target, planning_variants)
            current_mutations = list(base_variant)
            base_name         = planning_variants[base_variant]

            target_group_ids: List[int] = []
            for mut in target:
                mut_norm    = self.normalize_mut(mut)
                possible    = mutation_to_groups.get(mut_norm, [])
                filtered    = [
                    gid for gid in possible
                    if all(m in map(self.normalize_mut, target) for m in primer_groups[gid])
                    and any(m not in map(self.normalize_mut, current_mutations) for m in primer_groups[gid])
                ]
                if filtered:
                    best_gid = max(
                        filtered,
                        key=lambda gid: (
                            len([m for m in primer_groups[gid] if self.normalize_mut(m) in map(self.normalize_mut, target)]),
                            len(primer_groups[gid]),
                        ),
                    )
                    if best_gid not in target_group_ids:
                        target_group_ids.append(best_gid)

            needed = [
                gid for gid in target_group_ids
                if not all(m in map(self.normalize_mut, current_mutations) for m in primer_groups[gid])
            ]

            while needed:
                batch_ids  = needed[:self.max_mutations_per_step]
                batch_muts: List[str] = []
                for gid in batch_ids:
                    batch_muts.extend(primer_groups[gid])

                test_muts = sorted(
                    set(map(self.normalize_mut, current_mutations)) | set(batch_muts),
                    key=mutation_position,
                )
                test_sorted = tuple(test_muts)

                if test_sorted not in planning_variants:
                    new_name = f"{self.variant_prefix}{next_id_num:02d}"
                    planning_variants[test_sorted]   = new_name
                    used_variants[test_sorted]        = new_name
                    protocol_by_round[len(test_sorted)].append(
                        [new_name, base_name, ', '.join(batch_muts)]
                    )
                    base_name         = new_name
                    current_mutations = list(test_muts)
                    variant_databank[new_name] = ','.join(test_sorted)
                    next_id_num += 1

                for gid in batch_ids:
                    needed.remove(gid)

            final_variants[base_name].append(target)
            variant_to_label_map[variant_label]     = target
            variant_to_final_variant[variant_label] = base_name
            used_variants[
                tuple(sorted(map(self.normalize_mut, target), key=mutation_position))
            ] = base_name

        if self.list_existing_as_steps:
            existing_variants.update({
                muts: name
                for muts, name in planning_variants.items()
                if muts not in existing_variants
            })

        self.generate_pdf(
            protocol_by_round, final_variants, str(self.pdf_path),
            existing_input_variants, variant_to_label_map,
            variant_to_final_variant, used_variants,
        )
        self.save_protocol_json(
            protocol_by_round, final_variants, existing_input_variants,
            variant_to_label_map, variant_to_final_variant, used_variants,
        )

        # --- Databank persistence ---
        databank_file = self.output_dir / "variant_databank.json"
        if databank_file.exists():
            with open(databank_file, 'r') as f:
                try:
                    current_state = self._json.load(f)
                except self._json.JSONDecodeError:
                    current_state = {}
            self.undo_stack.append(copy.deepcopy(current_state))
        else:
            self.undo_stack.append({})

        if self.list_existing_as_steps:
            existing_json: Dict = {}
            if databank_file.exists():
                with open(databank_file, "r") as f:
                    try:
                        existing_json = self._json.load(f)
                    except self._json.JSONDecodeError:
                        pass
            existing_json.update(variant_databank)
            with open(databank_file, "w") as f:
                self._json.dump(existing_json, f, indent=2)
        else:
            with open(databank_file, "w") as f:
                self._json.dump(variant_databank, f, indent=2)

        print(f"PDF generated: {self.pdf_path}")
        print(f"Databank updated: {self.databank_file}")


# ===========================================================================
# LabelPrintingSystem — PDF label generation (pure reportlab, no GUI)
# ===========================================================================

class LabelPrintingSystem:
    """Generate PDF label sheets (Herma 10900 format)."""

    def __init__(self) -> None:
        from reportlab.lib.units import mm
        self.columns            = 7
        self.rows               = 27
        self.label_width        = 25.3 * mm
        self.label_height       = 9.9  * mm
        self.horizontal_spacing = 2.5  * mm
        self.top_margin         = 13.5 * mm
        self.left_margin        = 10.3 * mm
        self.font_size          = 5.5

    def generate_labels_pdf(
        self, selected_variants: List[str], start_row: int, start_col: int, output_path
    ) -> bool:
        from reportlab.lib.pagesizes import A4
        from reportlab.pdfgen import canvas as _canvas

        page_width, page_height = A4
        c = _canvas.Canvas(str(output_path), pagesize=A4)
        c.setFont("Helvetica", self.font_size)

        sorted_variants = sorted(selected_variants)
        variant_index   = 0
        total           = len(sorted_variants)
        current_row     = start_row
        current_col     = start_col

        while variant_index < total:
            x = self.left_margin + current_col * (self.label_width + self.horizontal_spacing)
            y = page_height - self.top_margin - (current_row + 1) * self.label_height

            c.drawCentredString(
                x + self.label_width / 2,
                y + self.label_height / 2 - self.font_size / 2,
                sorted_variants[variant_index],
            )
            variant_index += 1
            current_col += 1
            if current_col >= self.columns:
                current_col = 0
                current_row += 1
                if current_row >= self.rows and variant_index < total:
                    c.showPage()
                    c.setFont("Helvetica", self.font_size)
                    current_row = 0
                    current_col = 0

        c.save()
        return True


# ===========================================================================
# ConfigManager — pure file-system persistence, no GUI
# ===========================================================================

class ConfigManager:
    """Persist application configuration to ~/.mutagenesis/mutagenesis_config.json."""

    def __init__(self, config_file: str = "mutagenesis_config.json") -> None:
        from pathlib import Path
        self.config_file = Path.home() / ".mutagenesis" / config_file
        self.config_file.parent.mkdir(parents=True, exist_ok=True)
        self.config: Dict = self._load()

    def _load(self) -> Dict:
        import json
        if self.config_file.exists():
            try:
                with open(self.config_file, 'r') as f:
                    return json.load(f)
            except Exception:
                pass
        return {}

    def save_config(self) -> None:
        import json
        try:
            with open(self.config_file, 'w') as f:
                json.dump(self.config, f, indent=2)
        except Exception as e:
            print(f"Warning: Could not save config: {e}")

    def get(self, key: str, default=None):
        return self.config.get(key, default)

    def set(self, key: str, value) -> None:
        self.config[key] = value
        self.save_config()


# ===========================================================================
# Protein sequence formatting  (moved from BaseProteinPage mixin in GUI)
# ===========================================================================

def format_protein_sequence(protein: str, fmt: str) -> str:
    """
    Format a one-letter protein sequence for display.

    fmt is one of: '1-letter', '3-letter', 'Both', 'FASTA'.
    """
    if fmt == "1-letter":
        lines = [protein[i:i+10] for i in range(0, len(protein), 10)]
        rows  = [" ".join(lines[i:i+6]) for i in range(0, len(lines), 6)]
        return "\n".join(
            f"{i*60+1:>5}: {row}" for i, row in enumerate(rows)
        )
    elif fmt == "3-letter":
        three  = [get_three_letter_code(aa) for aa in protein]
        blocks = ["-".join(three[i:i+10]) for i in range(0, len(three), 10)]
        return "\n".join(blocks)
    elif fmt == "Both":
        one   = format_protein_sequence(protein, "1-letter")
        three = format_protein_sequence(protein, "3-letter")
        return f"=== 1-Letter ===\n{one}\n\n=== 3-Letter ===\n{three}"
    elif fmt == "FASTA":
        seq = "\n".join(protein[i:i+60] for i in range(0, len(protein), 60))
        return f">protein_sequence\n{seq}"
    return protein


# ===========================================================================
# Sequence alignment  (moved from MutationExtractorPage in GUI)
# ===========================================================================

def needleman_wunsch(seq1: str, seq2: str) -> Tuple[str, str, int]:
    """
    Global pairwise alignment (Needleman-Wunsch).
    Scoring: match +2, mismatch -1, gap -2.
    Returns (aligned_seq1, aligned_seq2, score).
    """
    match_s, mismatch_s, gap_s = 2, -1, -2
    n, m = len(seq1), len(seq2)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        dp[i][0] = i * gap_s
    for j in range(1, m + 1):
        dp[0][j] = j * gap_s
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = dp[i-1][j-1] + (
                match_s if seq1[i-1] == seq2[j-1] else mismatch_s)
            dp[i][j] = max(diag, dp[i-1][j] + gap_s, dp[i][j-1] + gap_s)

    a1: List[str] = []
    a2: List[str] = []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            score = dp[i][j]
            diag  = dp[i-1][j-1] + (
                match_s if seq1[i-1] == seq2[j-1] else mismatch_s)
            if score == diag:
                a1.append(seq1[i-1]); a2.append(seq2[j-1]); i -= 1; j -= 1
            elif score == dp[i-1][j] + gap_s:
                a1.append(seq1[i-1]); a2.append('-'); i -= 1
            else:
                a1.append('-'); a2.append(seq2[j-1]); j -= 1
        elif i > 0:
            a1.append(seq1[i-1]); a2.append('-'); i -= 1
        else:
            a1.append('-'); a2.append(seq2[j-1]); j -= 1

    return ''.join(reversed(a1)), ''.join(reversed(a2)), dp[n][m]


def find_mutations_from_alignment(
    ref_protein: str, variant_protein: str
) -> Tuple[List[str], float, int]:
    """
    Align variant_protein against ref_protein and return mutations.

    Returns:
        mutations  — list of strings like 'K49R'
        identity   — fraction of aligned positions that match (0-1)
        n_aligned  — number of aligned (non-gap) positions
    """
    a1, a2, _ = needleman_wunsch(ref_protein, variant_protein)
    ref_pos   = 0
    mutations: List[str] = []
    matched   = 0

    for r, v in zip(a1, a2):
        if r != '-':
            ref_pos += 1
        if r != '-' and v != '-':
            if r != v:
                mutations.append(f"{r}{ref_pos}{v}")
            else:
                matched += 1

    n_aligned = sum(1 for r, v in zip(a1, a2) if r != '-' and v != '-')
    identity  = matched / n_aligned if n_aligned else 0.0
    return mutations, identity, n_aligned


# ===========================================================================
# Variant hierarchy  (moved from MutationFailureHandler in GUI)
# ===========================================================================

def build_variant_hierarchy(databank: Dict[str, str]) -> Dict[str, str]:
    """
    Build a parent→child map from a variant databank.

    For each variant, finds the largest existing variant whose mutation
    set is a strict subset — i.e. the most-mutated direct ancestor.

    Returns {child_id: parent_id}.
    """
    sets: Dict[str, frozenset] = {
        vid: frozenset(m.strip() for m in muts.split(',') if m.strip())
        for vid, muts in databank.items()
        if muts and muts != "(none)"
    }
    hierarchy: Dict[str, str] = {}
    for child, child_muts in sorted(sets.items(), key=lambda x: len(x[1])):
        best_parent, best_size = None, 0
        for parent, parent_muts in sets.items():
            if (parent != child
                    and parent_muts.issubset(child_muts)
                    and len(parent_muts) > best_size):
                best_parent, best_size = parent, len(parent_muts)
        if best_parent:
            hierarchy[child] = best_parent
    return hierarchy
