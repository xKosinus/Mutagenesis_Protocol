"""
primer_gui.py
=============
All GUI windows and pages for the Site-Directed Mutagenesis Primer Designer.

Depends on:
  primer_core.py  — all calculations, data classes, protocol logic
  primer_io.py    — all file I/O (CSV, JSON, TXT, databank, WT sequences)

No business logic lives here: every calculation delegates to primer_core,
every file read/write delegates to primer_io.
"""

from __future__ import annotations

import copy
import re
import threading
import webbrowser
from pathlib import Path
from tkinter import filedialog, messagebox, simpledialog
from typing import Dict, List, Optional

import customtkinter as ctk
from customtkinter import CTkFont, DrawEngine

DrawEngine.preferred_drawing_method = "polygon_shapes"

# Restore original appearance — these were missing and caused the light-mode bug
ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")

# ---------------------------------------------------------------------------
# Core + IO imports — the ONLY place GUI talks to logic/files
# ---------------------------------------------------------------------------
import primer_core as core
import primer_io   as io

from primer_core import (
    GENETIC_CODE, _PRIMER3_OK,
    SLiCEFragment, SLiCEPrimer,
    PrimerGenerator, SLiCEPrimerDesigner,
    MutagenesisProtocol, LabelPrintingSystem,
    calculate_tm, calculate_gc_content, reverse_complement,
    translate_dna_to_protein, parse_mutation, mutation_position,
    natural_sort_key, PlasmidAnnotator,
    format_protein_sequence,
    needleman_wunsch, find_mutations_from_alignment,
    build_variant_hierarchy,
)


# ===========================================================================
# Shared / utility widgets
# ===========================================================================

class ScrollableFrameWithWheel(ctk.CTkScrollableFrame):
    """CTkScrollableFrame that responds to the mouse wheel on Windows & Linux."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        def _on_wheel(event):
            self._parent_canvas.yview_scroll(int(-1 * event.delta * 0.2), "units")

        def _on_up(event):
            self._parent_canvas.yview_scroll(-1, "units")

        def _on_down(event):
            self._parent_canvas.yview_scroll(1, "units")

        self.bind("<Enter>", lambda e: [
            self.bind_all("<MouseWheel>", _on_wheel),
            self.bind_all("<Button-4>",   _on_up),
            self.bind_all("<Button-5>",   _on_down),
        ])
        self.bind("<Leave>", lambda e: [
            self.unbind_all("<MouseWheel>"),
            self.unbind_all("<Button-4>"),
            self.unbind_all("<Button-5>"),
        ])


class PlaceholderTextbox(ctk.CTkTextbox):
    """CTkTextbox with placeholder text support."""

    def __init__(self, master, placeholder="", placeholder_fg="gray",
                 text_fg="white", **kwargs):
        super().__init__(master, **kwargs)
        self._placeholder      = placeholder
        self._placeholder_fg   = placeholder_fg
        self._text_fg          = text_fg
        self._showing_placeholder = False
        self._set_placeholder()
        self.bind("<FocusIn>",  self._clear_placeholder)
        self.bind("<FocusOut>", self._restore_placeholder)

    def _set_placeholder(self):
        if not self.get("1.0", "end").strip():
            self.insert("1.0", self._placeholder)
            self.configure(text_color=self._placeholder_fg)
            self._showing_placeholder = True

    def _clear_placeholder(self, event=None):
        if self._showing_placeholder:
            self.delete("1.0", "end")
            self.configure(text_color=self._text_fg)
            self._showing_placeholder = False

    def _restore_placeholder(self, event=None):
        if not self.get("1.0", "end").strip():
            self._set_placeholder()


class GithubButton(ctk.CTkButton):
    """Button that opens a GitHub URL in the default browser."""

    def __init__(self, master, url: str, **kwargs):
        super().__init__(master, command=lambda: webbrowser.open(url), **kwargs)


class BaseProteinPage:
    """Mixin marker for protein display pages.

    format_protein_sequence() has moved to primer_core and is imported
    at the top of this file — pages call it as a plain function.
    """


# ===========================================================================
# Primer3 warning banner (shown once at startup when _PRIMER3_OK is False)
# ===========================================================================

class Primer3WarningBanner(ctk.CTkFrame):
    """Orange banner shown when primer3 is not working."""

    def __init__(self, master, **kwargs):
        super().__init__(master, fg_color="#B35900", corner_radius=0, **kwargs)
        msg = (
            "⚠  primer3 is not working correctly on this system — "
            "Tm values will show as 0.0.   "
            "Fix: open a terminal and run:  pip install primer3-py --force-reinstall"
        )
        ctk.CTkLabel(self, text=msg, text_color="white",
                     font=ctk.CTkFont(size=12)).pack(padx=20, pady=6)


# ===========================================================================
# MutagenesisApp  (main window)
# ===========================================================================

class MutagenesisApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("Site-Directed Mutagenesis Primer Designer")
        self.geometry("1050x500")
        self.minsize(1200, 700)

        self.config_manager = core.ConfigManager()

        self.LARGEFONT  = ctk.CTkFont(family="Verdana", size=24, weight="bold")
        self.MEDIUMFONT = ctk.CTkFont(family="Verdana", size=14, weight="bold")
        self.SMALLFONT  = ctk.CTkFont(family="Verdana", size=12)

        last_dir = self.config_manager.get(
            'last_output_dir', str(Path.cwd() / "Mutagenesis")
        )
        self.output_dir = ctk.StringVar(value=last_dir)
        self.output_dir.trace_add('write', lambda *_: self._save_output_dir())

        self.modified_variants: set = set()
        self.repair_variants:   set = set()
        self.undo_stack:        list = []

        self._build_layout()

    # ------------------------------------------------------------------
    # Layout
    # ------------------------------------------------------------------

    def _build_layout(self):
        # Optional primer3 warning banner
        if not _PRIMER3_OK:
            banner = Primer3WarningBanner(self)
            banner.pack(fill="x", side="top")

        # Two-column main layout
        main = ctk.CTkFrame(self)
        main.pack(fill="both", expand=True)
        main.grid_columnconfigure(0, weight=0, minsize=200)
        main.grid_columnconfigure(1, weight=1)
        main.grid_rowconfigure(0, weight=1)

        # Navigation
        nav = ctk.CTkFrame(main)
        nav.grid(row=0, column=0, sticky="nsew", padx=(5, 2), pady=5)
        nav.grid_propagate(False)
        nav.grid_columnconfigure(0, weight=1)

        nav_buttons = [
            ("User Input",         InputPage),
            ("Variant Databank",   DatabankPage),
            ("Wildtype Protein",   WildtypeProteinPage),
            ("Variant Proteins",   VariantProteinPage),
            ("Protocol Results",   ProtocolResultsPage),
            ("Primer",             PrimerPage),
            ("Mutation Extractor", MutationExtractorPage),
            ("Mutation Eraser",    MutationFailureHandler),
            ("SLiCE Designer",     SLiCEDesignerPage),
        ]
        nav.grid_rowconfigure(len(nav_buttons), weight=1)
        for i, (text, page_cls) in enumerate(nav_buttons):
            btn = ctk.CTkButton(
                nav, text=text,
                command=lambda p=page_cls: self.show_frame(p)
            )
            btn.grid(row=i, column=0, sticky="ew",
                     pady=(15 if i == 0 else 5, 5), padx=15)

        GithubButton(
            nav,
            "https://github.com/xKosinus/Mutagenesis_Protocol/tree/windows_gui",
            text="GitHub Source",
        ).grid(row=len(nav_buttons) + 1, column=0, sticky="ew", padx=15, pady=20)

        # Page container
        container = ctk.CTkFrame(main)
        container.grid(row=0, column=1, sticky="nsew", padx=(2, 5), pady=5)
        container.grid_columnconfigure(0, weight=1)
        container.grid_rowconfigure(0, weight=1)

        self.frames: Dict = {}
        for page_cls in [p for _, p in nav_buttons]:
            frame = page_cls(container, self)
            self.frames[page_cls] = frame
            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(InputPage)

    # ------------------------------------------------------------------
    # Navigation
    # ------------------------------------------------------------------

    def show_frame(self, page_class):
        frame = self.frames[page_class]
        frame.tkraise()
        if page_class == WildtypeProteinPage:
            frame.update_protein_display()
        elif page_class == VariantProteinPage:
            frame.update_variant_display()
        elif page_class == PrimerPage:
            frame.load_primers()
        elif page_class == MutationExtractorPage:
            frame.refresh_reference_sequences()

    # ------------------------------------------------------------------
    # Databank helpers (delegate to primer_io)
    # ------------------------------------------------------------------

    def get_databank(self) -> Dict[str, str]:
        return io.load_databank(Path(self.output_dir.get()) / "protocols")

    def save_databank(self, databank: Dict[str, str]) -> None:
        io.save_databank(databank, Path(self.output_dir.get()) / "protocols")

    # ------------------------------------------------------------------
    # Variant highlight helpers
    # ------------------------------------------------------------------

    def mark_variants_modified(self, updated_ids, new_repair_ids):
        self.modified_variants.update(updated_ids)
        self.repair_variants.update(new_repair_ids)

    def is_variant_modified(self, variant_id: str) -> bool:
        return variant_id in self.modified_variants

    def is_variant_repair(self, variant_id: str) -> bool:
        return variant_id in self.repair_variants

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _save_output_dir(self):
        self.config_manager.set('last_output_dir', self.output_dir.get())


# ===========================================================================
# InputPage
# ===========================================================================

class InputPage(ctk.CTkFrame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        self.sequence_undo_stack: list = []

        self.grid_columnconfigure(0, weight=0, minsize=200)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=0)
        self.grid_rowconfigure(2, weight=1)
        self.grid_rowconfigure(4, weight=2)

        self._build_ui()
        self.after(200, self._restore_last_sequence)

    # ------------------------------------------------------------------
    def _build_ui(self):
        # Title
        ctk.CTkLabel(self, text="Mutagenesis Input", corner_radius=10,
                     fg_color="#4a90e2", text_color="white",
                     font=self.controller.LARGEFONT, height=50
                     ).grid(row=0, column=0, padx=10, pady=20,
                            columnspan=3, sticky="ew")

        # Output directory
        ctk.CTkLabel(self, text="Output Directory:",
                     font=self.controller.MEDIUMFONT
                     ).grid(row=1, column=0, sticky="w", padx=10, pady=5)
        ctk.CTkEntry(self, textvariable=self.controller.output_dir
                     ).grid(row=1, column=1, sticky="ew", padx=5, pady=5)
        ctk.CTkButton(self, text="Choose...",
                      command=self._select_output_dir, width=70
                      ).grid(row=1, column=2, sticky="w", padx=5, pady=5)

        # DNA sequence
        ctk.CTkLabel(self, text="Wildtype DNA Sequence:",
                     font=self.controller.MEDIUMFONT
                     ).grid(row=2, column=0, sticky="nw", padx=10, pady=(15, 0))
        seq_frame = ctk.CTkFrame(self)
        seq_frame.grid(row=2, column=1, columnspan=2, pady=5,
                       padx=5, sticky="nsew")
        seq_frame.grid_columnconfigure(0, weight=1)
        seq_frame.grid_rowconfigure(0, weight=1)
        self.seq_text = PlaceholderTextbox(
            seq_frame,
            placeholder="Enter wildtype DNA sequence here...",
            placeholder_fg="gray", text_fg="white", height=150,
        )
        self.seq_text.grid(row=0, column=0, sticky="nsew", padx=(0, 5))
        ctk.CTkButton(seq_frame, text="?", width=30, height=30,
                      command=lambda: messagebox.showinfo(
                          "Info", "Include 30 bases at the flanking sites!")
                      ).grid(row=0, column=1, sticky="ne", pady=5)

        # Sequence management buttons
        btn_frame = ctk.CTkFrame(self)
        btn_frame.grid(row=3, column=1, columnspan=3, sticky="ew", padx=5, pady=5)
        for i, (txt, cmd) in enumerate([
            ("Save Sequence",          self._save_sequence),
            ("Load/Manage Sequences",  self._load_manage_sequences),
            ("Undo Sequence Changes",  self._undo_sequence_changes),
        ]):
            ctk.CTkButton(btn_frame, text=txt, command=cmd
                          ).grid(row=0, column=i, padx=10, pady=5, sticky="w")

        # Mutations
        ctk.CTkLabel(self, text="Mutations per Variant:",
                     font=self.controller.MEDIUMFONT
                     ).grid(row=4, column=0, sticky="nw", padx=10, pady=(10, 0))
        var_frame = ctk.CTkFrame(self)
        var_frame.grid(row=4, column=1, columnspan=2, pady=5,
                       padx=5, sticky="nsew")
        var_frame.grid_columnconfigure(0, weight=1)
        var_frame.grid_rowconfigure(0, weight=1)
        self.var_text = ctk.CTkTextbox(var_frame, height=150,
                                       font=self.controller.SMALLFONT)
        self.var_text.grid(row=0, column=0, sticky="nsew", padx=(0, 5))
        self.var_text.insert("1.0",
                             "S19V, S30V, Q424V, A431E\nK3Q, K36R, D76Q, L167A, Q424V\n")
        ctk.CTkButton(var_frame, text="?", width=30, height=30,
                      command=lambda: messagebox.showinfo(
                          "Info",
                          "Enter mutations separated by commas, one variant per line.\n"
                          "Example: S19V, S30V\nMultiple variants separated by new lines.")
                      ).grid(row=0, column=1, sticky="ne", pady=5)

        # Parameters
        params = ctk.CTkFrame(self)
        params.grid(row=5, column=0, columnspan=3, sticky="ew", padx=10, pady=10)
        params.grid_columnconfigure(1, weight=1)
        params.grid_columnconfigure(3, weight=1)

        ctk.CTkLabel(params, text="Max mutations per step:",
                     font=self.controller.MEDIUMFONT
                     ).grid(row=0, column=0, sticky="w", padx=10, pady=5)
        self.max_mut_var = ctk.IntVar(value=1)
        ctk.CTkComboBox(params, values=["1", "2"],
                        variable=self.max_mut_var, width=80
                        ).grid(row=0, column=1, sticky="w", padx=5, pady=5)

        ctk.CTkLabel(params, text="Variant prefix:",
                     font=self.controller.MEDIUMFONT
                     ).grid(row=0, column=2, sticky="w", padx=10, pady=5)
        self.var_prefix = ctk.StringVar(value="KWE_TA_A")
        ctk.CTkEntry(params, textvariable=self.var_prefix, width=150
                     ).grid(row=0, column=3, sticky="w", padx=5, pady=5)

        self.skip_db = ctk.BooleanVar(value=False)
        ctk.CTkCheckBox(params, text="Skip variant databank",
                        variable=self.skip_db
                        ).grid(row=0, column=4, pady=10, padx=10, sticky="e")

        # Action buttons
        acts = ctk.CTkFrame(self)
        acts.grid(row=6, column=0, columnspan=3, sticky="ew", padx=10, pady=10)
        ctk.CTkButton(acts, text="Run Complete Workflow",
                      font=self.controller.MEDIUMFONT, fg_color="green",
                      command=self._run_workflow
                      ).grid(row=0, column=0, padx=10, pady=10, sticky="w")
        ctk.CTkButton(acts, text="Undo Last Change",
                      font=self.controller.MEDIUMFONT, fg_color="orange",
                      command=self._undo_change
                      ).grid(row=0, column=1, padx=10, pady=10, sticky="w")

        # Status
        self.status_label = ctk.CTkLabel(
            self, text="Ready to run workflow",
            font=CTkFont(family="Verdana", size=14), text_color="gray"
        )
        self.status_label.grid(row=7, column=0, columnspan=3,
                               sticky="w", padx=10, pady=(0, 10))

    # ------------------------------------------------------------------
    # Sequence helpers
    # ------------------------------------------------------------------

    def get_dna_sequence(self) -> str:
        seq = self.seq_text.get("1.0", "end").strip()
        if seq == "Enter wildtype DNA sequence here...":
            return ""
        return seq.replace('\n', '').upper()

    def _output_dir(self) -> Path:
        return Path(self.controller.output_dir.get())

    # ------------------------------------------------------------------
    # Sequence save / load / undo
    # ------------------------------------------------------------------

    def _save_sequence(self):
        current_seq = self.get_dna_sequence()
        if not current_seq:
            messagebox.showwarning("Warning",
                                   "Please enter a DNA sequence before saving.")
            return
        while True:
            name = simpledialog.askstring(
                "Save Sequence", "Enter a unique name for this sequence:")
            if name is None:
                return
            name = name.strip()
            if not name:
                messagebox.showwarning("Warning", "Please enter a non-empty name.")
                continue
            seqs = io.load_wt_sequences(self._output_dir())
            if name in seqs:
                if not messagebox.askyesno(
                    "Name Exists",
                    f"A sequence named '{name}' already exists. Overwrite?"
                ):
                    continue
            seqs[name] = current_seq
            io.save_wt_sequences_undo(seqs, self._output_dir(),
                                      self.sequence_undo_stack)
            if len(self.sequence_undo_stack) > 10:
                self.sequence_undo_stack.pop(0)
            self.status_label.configure(
                text=f"Sequence '{name}' saved successfully", text_color="green")
            break

    def _load_manage_sequences(self):
        seqs = io.load_wt_sequences(self._output_dir())
        if not seqs:
            messagebox.showinfo("Info", "No saved sequences found.")
            return

        win = ctk.CTkToplevel(self)
        win.title("Manage Saved Sequences")
        win.geometry("800x500")
        win.transient(self)
        win.after(100, win.grab_set)
        win.grid_columnconfigure(0, weight=1)
        win.grid_rowconfigure(1, weight=1)

        ctk.CTkLabel(win, text="Saved Sequences",
                     font=self.controller.LARGEFONT
                     ).grid(row=0, column=0, pady=20, sticky="ew")

        list_frame = ScrollableFrameWithWheel(win)
        list_frame.grid(row=1, column=0, sticky="nsew", padx=20, pady=10)
        list_frame.grid_columnconfigure(1, weight=1)

        sequence_vars: Dict = {}
        load_var = ctk.StringVar(value="")

        for col, text in [(0, "Load"), (1, "Delete"), (2, "Name"),
                          (3, "Sequence (first 50 chars)")]:
            ctk.CTkLabel(list_frame, text=text,
                         font=self.controller.MEDIUMFONT
                         ).grid(row=0, column=col, padx=10, pady=5, sticky="w")

        for i, (name, sequence) in enumerate(sorted(seqs.items()), start=1):
            ctk.CTkRadioButton(list_frame, text="",
                               variable=load_var, value=name
                               ).grid(row=i, column=0, padx=10, pady=2, sticky="w")
            dv = ctk.BooleanVar(value=False)
            sequence_vars[name] = dv
            ctk.CTkCheckBox(list_frame, text="",
                            variable=dv
                            ).grid(row=i, column=1, padx=10, pady=2, sticky="w")
            ctk.CTkLabel(list_frame, text=name,
                         font=self.controller.SMALLFONT
                         ).grid(row=i, column=2, padx=10, pady=2, sticky="w")
            preview = sequence[:50] + ("..." if len(sequence) > 50 else "")
            ctk.CTkLabel(list_frame, text=preview,
                         font=ctk.CTkFont(family="Courier", size=10)
                         ).grid(row=i, column=3, padx=10, pady=2, sticky="w")

        btn_frame = ctk.CTkFrame(win)
        btn_frame.grid(row=2, column=0, sticky="ew", padx=20, pady=10)

        def _load_selected():
            sel = load_var.get()
            if not sel:
                messagebox.showwarning("Warning",
                                       "Please select a sequence to load.")
                return
            self.seq_text._clear_placeholder()
            self.seq_text.delete("1.0", "end")
            self.seq_text.insert("1.0", seqs[sel])
            self.status_label.configure(
                text=f"Loaded sequence '{sel}'", text_color="green")
            win.destroy()

        def _delete_selected():
            to_del = [n for n, v in sequence_vars.items() if v.get()]
            if not to_del:
                messagebox.showwarning("Warning",
                                       "Please select sequences to delete.")
                return
            if not messagebox.askyesno(
                "Confirm Deletion",
                f"Delete {len(to_del)} sequence(s)?"
            ):
                return
            updated = {n: s for n, s in seqs.items() if n not in to_del}
            io.save_wt_sequences_undo(updated, self._output_dir(),
                                      self.sequence_undo_stack)
            self.status_label.configure(
                text=f"Deleted {len(to_del)} sequence(s)", text_color="orange")
            win.destroy()

        for col, (txt, cmd) in enumerate([
            ("Load Selected",   _load_selected),
            ("Delete Selected", _delete_selected),
            ("Cancel",          win.destroy),
        ]):
            ctk.CTkButton(btn_frame, text=txt, command=cmd
                          ).grid(row=0, column=col, padx=10, pady=10, sticky="w")

    def _undo_sequence_changes(self):
        if not self.sequence_undo_stack:
            self.status_label.configure(
                text="No sequence changes to undo", text_color="orange")
            return
        last = self.sequence_undo_stack.pop()
        io.save_wt_sequences(last, self._output_dir())
        self.status_label.configure(
            text="Sequence changes undone successfully", text_color="green")

    def _restore_last_sequence(self):
        last = self.controller.config_manager.get('last_dna_sequence', '')
        if last:
            self.seq_text._clear_placeholder()
            self.seq_text.delete("1.0", "end")
            self.seq_text.insert("1.0", last)
            self.seq_text.configure(text_color="white")
            self.status_label.configure(
                text="Last used DNA sequence restored", text_color="gray")

    def _select_output_dir(self):
        new_dir = filedialog.askdirectory(
            title="Select Output Directory",
            initialdir=self.controller.output_dir.get(),
        )
        if new_dir:
            self.controller.output_dir.set(new_dir)
            self.status_label.configure(
                text=f"Output directory changed to: {new_dir}",
                text_color="green")

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def _validate_inputs(self) -> bool:
        dna = self.get_dna_sequence()
        raw = self.var_text.get("1.0", "end").strip()
        lines = [l.strip() for l in raw.splitlines() if l.strip()]
        mutations = [[m.strip() for m in l.split(',') if m.strip()]
                     for l in lines]

        if not dna:
            messagebox.showerror("Input Error",
                                 "DNA sequence is empty.")
            return False
        if not re.fullmatch(r"[ATGC]+", dna):
            messagebox.showerror("Input Error",
                                 "DNA sequence must only contain A, T, G, C.")
            return False
        if len(dna) % 3 != 0:
            messagebox.showerror("Input Error",
                                 "DNA sequence length must be divisible by 3.")
            return False
        if len(dna) < 66:
            messagebox.showerror("Input Error",
                                 "DNA sequence is too short (minimum 66 bases).")
            return False
        stop = dna[-33:-30]
        if stop not in {"TAA", "TAG", "TGA"}:
            messagebox.showerror("Input Error",
                                 f"Expected stop codon before last 30 bases, "
                                 f"found '{stop}'.")
            return False

        first_codon = dna[30:33]
        first_aa    = GENETIC_CODE.get(first_codon, 'X')
        self.status_label.configure(
            text=f"First codon after flanking region: {first_codon} "
                 f"({first_aa}) — please verify",
            text_color="orange",
        )

        if not mutations:
            messagebox.showerror("Input Error",
                                 "No mutations entered.")
            return False

        pattern = re.compile(
            r"^[ACDEFGHIKLMNPQRSTVWY]\d+[ACDEFGHIKLMNPQRSTVWY]$"
        )
        for idx, variant in enumerate(mutations):
            seen: set = set()
            for mut in variant:
                if not pattern.match(mut):
                    messagebox.showerror(
                        "Input Error",
                        f"Invalid mutation format: '{mut}' in line {idx+1}. "
                        "Expected format e.g. K49R.")
                    return False
                if mut in seen:
                    messagebox.showerror(
                        "Input Error",
                        f"Duplicate mutation '{mut}' in line {idx+1}.")
                    return False
                seen.add(mut)

        self.status_label.configure(
            text="Validation successful.", text_color="green")
        return True

    # ------------------------------------------------------------------
    # Workflow
    # ------------------------------------------------------------------

    def _run_workflow(self):
        self.status_label.configure(text="Running workflow...",
                                    text_color="blue")
        if not self._validate_inputs():
            return

        seq      = self.get_dna_sequence()
        self.controller.config_manager.set('last_dna_sequence', seq)
        raw      = self.var_text.get("1.0", "end").strip().split('\n')
        variants = [list(map(str.strip, l.split(',')))
                    for l in raw if l.strip()]
        max_muts = self.max_mut_var.get()
        vp       = self.var_prefix.get().strip()
        skip_db  = self.skip_db.get()
        out_dir  = Path(self.controller.output_dir.get())

        # Primer design
        try:
            pg = PrimerGenerator(flank_size_bases=30)
            valid, codon_info = pg.check_sequence_validity(seq, flank_size=30)
            if not valid:
                self.status_label.configure(
                    text="Sequence failed validation", text_color="red")
                return

            if codon_info:
                start_codon, first_aa = codon_info
                messagebox.showinfo(
                    "Sequence Check",
                    f"First codon: {start_codon} ({first_aa})\n"
                    "Please verify this is the correct reading frame."
                )

            parsed: list = []
            parsed_strings: list = []   # string form for MutagenesisProtocol
            for variant in variants:
                row = []
                row_str = []
                for mut in variant:
                    aa, pos, new_aa = parse_mutation(mut)
                    row.append((aa, pos, new_aa))
                    row_str.append(f"{aa}{pos}{new_aa}")
                parsed.append(row)
                parsed_strings.append(row_str)

            pg.validate_mutations_against_sequence(seq, parsed_strings)

            # IMPORTANT: generate primers per-variant, not all at once.
            # generate_primers_for_construct() flattens its whole input and
            # merges nearby mutations together (merge_close_mutations), so
            # passing every variant in a single call lets mutations from
            # *different* variants merge into the same primer group just
            # because they sit close together in the sequence. Each variant
            # needs its own call, mirroring the original workflow.
            primer_data: list = []
            existing_pairs: set = set()
            for row in parsed:
                for primer_set in pg.generate_primers_for_construct(seq, row):
                    key = (primer_set['forward_primer'], primer_set['reverse_primer'])
                    if key not in existing_pairs:
                        primer_data.append(primer_set)
                        existing_pairs.add(key)

            if primer_data:
                io.save_primer_files(primer_data, out_dir)

        except Exception as e:
            self.status_label.configure(
                text=f"Primer design error: {e}", text_color="red")
            return

        # Protocol
        try:
            protocol = MutagenesisProtocol(
                variant_input          = parsed_strings,
                max_mutations_per_step = max_muts,
                variant_prefix         = vp,
                output_dir_path        = str(out_dir),
                list_existing_as_steps = skip_db,
                undo_stack             = self.controller.undo_stack,
            )
            protocol.run()
            self.status_label.configure(
                text="Workflow completed successfully! Check output directory.",
                text_color="green")
            self.controller.frames[DatabankPage].refresh_variants()

        except Exception as e:
            self.status_label.configure(
                text=f"Error: {e}", text_color="red")

    def _undo_change(self):
        out_dir      = Path(self.controller.output_dir.get())
        protocols_dir = out_dir / "protocols"
        try:
            if self.controller.undo_stack:
                last = self.controller.undo_stack.pop()
                io.save_databank(last, protocols_dir)
                io.remove_primer_files(out_dir)
                io.remove_protocol_files(protocols_dir)
                self.status_label.configure(
                    text="Last change undone successfully",
                    text_color="green")
            else:
                self.status_label.configure(
                    text="No changes to undo", text_color="orange")
        except Exception as e:
            self.status_label.configure(
                text=f"Undo failed: {e}", text_color="red")


# ===========================================================================
# WildtypeProteinPage
# ===========================================================================

class WildtypeProteinPage(ctk.CTkFrame, BaseProteinPage):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        self.current_protein_1letter = ""

        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(3, weight=1)

        # Header
        hf = ctk.CTkFrame(self)
        hf.grid(row=0, column=0, sticky="ew", padx=10, pady=20)
        hf.grid_columnconfigure(0, weight=1)
        ctk.CTkLabel(hf, text="Wildtype Protein Sequence",
                     corner_radius=10, fg_color="#4a90e2",
                     text_color="white", font=controller.LARGEFONT, height=50
                     ).grid(row=0, column=0, padx=10, pady=10, sticky="ew")
        ctk.CTkButton(hf, text="Back to Input",
                      command=lambda: controller.show_frame(InputPage)
                      ).grid(row=0, column=1, padx=10, pady=10, sticky="e")

        # Info
        info = ctk.CTkFrame(self)
        info.grid(row=1, column=0, sticky="ew", padx=10, pady=10)
        info.grid_columnconfigure((0, 1, 2), weight=1)
        self.dna_length_label    = ctk.CTkLabel(info, text="DNA Length: N/A",
                                                font=controller.MEDIUMFONT)
        self.coding_length_label = ctk.CTkLabel(info, text="Coding Length: N/A",
                                                font=controller.MEDIUMFONT)
        self.protein_length_label = ctk.CTkLabel(info, text="Protein Length: N/A",
                                                 font=controller.MEDIUMFONT)
        for col, lbl in enumerate([self.dna_length_label,
                                    self.coding_length_label,
                                    self.protein_length_label]):
            lbl.grid(row=0, column=col, padx=10, pady=10, sticky="w")

        # Format controls
        fmt = ctk.CTkFrame(self)
        fmt.grid(row=2, column=0, sticky="ew", padx=10, pady=5)
        ctk.CTkLabel(fmt, text="Display Format:",
                     font=controller.MEDIUMFONT
                     ).grid(row=0, column=0, padx=10, pady=10, sticky="w")
        self.format_var = ctk.StringVar(value="1-letter")
        ctk.CTkOptionMenu(fmt, values=["1-letter", "3-letter", "Both", "FASTA"],
                          variable=self.format_var,
                          command=self._update_display_format, width=120
                          ).grid(row=0, column=1, padx=5, pady=10, sticky="w")
        ctk.CTkButton(fmt, text="Refresh",
                      command=self.update_protein_display, width=80
                      ).grid(row=0, column=2, padx=10, pady=10, sticky="w")
        ctk.CTkButton(fmt, text="Copy to Clipboard",
                      command=self._copy_to_clipboard, width=120
                      ).grid(row=0, column=3, padx=10, pady=10, sticky="w")

        # Display
        self.protein_text = ctk.CTkTextbox(
            self, font=ctk.CTkFont(family="Courier", size=11))
        self.protein_text.grid(row=3, column=0, pady=10, padx=10, sticky="nsew")
        self.protein_text.configure(state="disabled")

        self.status_label = ctk.CTkLabel(self, text="",
                                         text_color="gray",
                                         font=controller.SMALLFONT)
        self.status_label.grid(row=4, column=0, sticky="w", padx=10, pady=5)

    def update_protein_display(self):
        dna = self.controller.frames[InputPage].get_dna_sequence()
        self.protein_text.configure(state="normal")
        self.protein_text.delete("1.0", "end")
        if not dna or len(dna) <= 60:
            self.protein_text.insert("1.0",
                "No DNA sequence provided or sequence too short.")
            self.protein_text.configure(state="disabled")
            for lbl in [self.dna_length_label,
                        self.coding_length_label,
                        self.protein_length_label]:
                lbl.configure(text=lbl.cget("text").split(":")[0] + ": N/A")
            return

        coding = dna[30:-30]
        self.current_protein_1letter = translate_dna_to_protein(coding)
        self.dna_length_label.configure(text=f"DNA Length: {len(dna)}bp")
        self.coding_length_label.configure(text=f"Coding Length: {len(coding)}bp")
        self.protein_length_label.configure(
            text=f"Protein Length: {len(self.current_protein_1letter)}aa")
        self.protein_text.configure(state="disabled")
        self._update_display_format()
        self.status_label.configure(
            text="✓ Removed 30bp flanking regions from both ends")

    def _update_display_format(self, *_):
        if not self.current_protein_1letter:
            return
        fmt = format_protein_sequence(
            self.current_protein_1letter, self.format_var.get())
        self.protein_text.configure(state="normal")
        self.protein_text.delete("1.0", "end")
        self.protein_text.insert("1.0", fmt)
        self.protein_text.configure(state="disabled")

    def _copy_to_clipboard(self):
        if self.current_protein_1letter:
            self.clipboard_clear()
            self.clipboard_append(self.protein_text.get("1.0", "end-1c"))
            self.status_label.configure(
                text="✓ Copied to clipboard", text_color="green")
        else:
            self.status_label.configure(
                text="No protein sequence to copy", text_color="orange")


# ===========================================================================
# VariantProteinPage
# ===========================================================================

class VariantProteinPage(ctk.CTkFrame, BaseProteinPage):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        self.current_variants: Dict = {}

        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(3, weight=1)

        hf = ctk.CTkFrame(self)
        hf.grid(row=0, column=0, sticky="ew", padx=10, pady=20)
        hf.grid_columnconfigure(0, weight=1)
        ctk.CTkLabel(hf, text="Variant Protein Sequences",
                     corner_radius=10, fg_color="#4a90e2",
                     text_color="white", font=controller.LARGEFONT, height=50
                     ).grid(row=0, column=0, padx=10, pady=10, sticky="ew")
        ctk.CTkButton(hf, text="Back to Databank",
                      command=lambda: controller.show_frame(DatabankPage)
                      ).grid(row=0, column=1, padx=10, pady=10, sticky="e")

        info = ctk.CTkFrame(self)
        info.grid(row=1, column=0, sticky="ew", padx=10, pady=10)
        info.grid_columnconfigure((0, 1, 2), weight=1)
        self.variant_count_label    = ctk.CTkLabel(info, text="Selected Variants: 0",
                                                   font=controller.MEDIUMFONT)
        self.total_sequences_label  = ctk.CTkLabel(info, text="Total Sequences: 0",
                                                   font=controller.MEDIUMFONT)
        self.avg_length_label       = ctk.CTkLabel(info, text="Avg Length: N/A",
                                                   font=controller.MEDIUMFONT)
        for col, lbl in enumerate([self.variant_count_label,
                                    self.total_sequences_label,
                                    self.avg_length_label]):
            lbl.grid(row=0, column=col, padx=10, pady=10, sticky="w")

        ctrl = ctk.CTkFrame(self)
        ctrl.grid(row=2, column=0, sticky="ew", padx=10, pady=5)
        ctk.CTkButton(ctrl, text="Refresh from Databank",
                      command=self.update_variant_display, width=150
                      ).grid(row=0, column=0, padx=10, pady=10, sticky="w")
        ctk.CTkButton(ctrl, text="Copy All to Clipboard",
                      command=self._copy_to_clipboard, width=150
                      ).grid(row=0, column=1, padx=10, pady=10, sticky="w")

        self.protein_text = ctk.CTkTextbox(
            self, font=ctk.CTkFont(family="Courier", size=11))
        self.protein_text.grid(row=3, column=0, pady=10, padx=10, sticky="nsew")
        self.protein_text.configure(state="disabled")

        self.status_label = ctk.CTkLabel(self, text="",
                                         text_color="gray",
                                         font=controller.SMALLFONT)
        self.status_label.grid(row=4, column=0, sticky="w", padx=10, pady=5)

    def _get_selected_variants(self) -> Dict:
        db_page = self.controller.frames.get(DatabankPage)
        if not db_page:
            return {}
        return {vid: muts for vid, cv in db_page.check_vars
                if cv.get() == 1
                for muts in [db_page.all_variants.get(vid, "")]}

    def update_variant_display(self):
        input_page = self.controller.frames[InputPage]
        dna = input_page.get_dna_sequence()
        if not dna or len(dna) <= 60:
            self.status_label.configure(
                text="No DNA sequence available", text_color="orange")
            return

        coding_dna     = dna[30:-30]
        wt_protein     = translate_dna_to_protein(coding_dna)
        selected       = self._get_selected_variants()
        databank       = self.controller.get_databank()
        variants_to_show = selected if selected else databank

        if not variants_to_show:
            self.status_label.configure(
                text="No variants found in databank", text_color="orange")
            return

        fasta_lines: List[str] = []
        self.current_variants = {}
        pg = PrimerGenerator(flank_size_bases=30)

        for var_id, mut_str in sorted(variants_to_show.items(),
                                      key=lambda x: natural_sort_key(x[0])):
            if mut_str and mut_str != "(none)":
                try:
                    muts = [parse_mutation(m.strip())
                            for m in mut_str.split(',') if m.strip()]
                    mut_seq = pg._create_mutated_sequence(dna, muts)
                    coding  = mut_seq[30:-30]
                except Exception:
                    coding = coding_dna
            else:
                coding = coding_dna

            protein = translate_dna_to_protein(coding)
            self.current_variants[var_id] = {
                "protein": protein, "length": len(protein)
            }
            seq60 = "\n".join(protein[i:i+60]
                              for i in range(0, len(protein), 60))
            fasta_lines.append(f">{var_id} | {mut_str or 'WT'}\n{seq60}")

        self.protein_text.configure(state="normal")
        self.protein_text.delete("1.0", "end")
        self.protein_text.insert("1.0", "\n\n".join(fasta_lines))
        self.protein_text.configure(state="disabled")

        cnt     = len(self.current_variants)
        lengths = [d["length"] for d in self.current_variants.values()]
        avg     = sum(lengths) / len(lengths) if lengths else 0
        self.variant_count_label.configure(text=f"Selected Variants: {cnt}")
        self.total_sequences_label.configure(text=f"Total Sequences: {cnt}")
        self.avg_length_label.configure(text=f"Avg Length: {avg:.1f} aa")
        self.status_label.configure(
            text=f"✓ Generated {cnt} variant protein sequences",
            text_color="green")

    def _copy_to_clipboard(self):
        if not self.current_variants:
            self.status_label.configure(
                text="No variant sequences to copy", text_color="orange")
            return
        self.clipboard_clear()
        self.clipboard_append(self.protein_text.get("1.0", "end-1c"))
        self.status_label.configure(
            text="✓ Copied all sequences to clipboard", text_color="green")


# ===========================================================================
# DatabankPage
# ===========================================================================

class DatabankPage(ctk.CTkFrame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        self.last_deleted: Dict = {}
        self.check_vars:   list = []
        self.checkboxes:   list = []
        self.all_variants:      Dict = {}
        self.filtered_variants: Dict = {}

        self.grid_columnconfigure(0, weight=1, minsize=400)
        self.grid_columnconfigure(1, weight=1, minsize=400)
        self.grid_rowconfigure(1, weight=1)

        self._build_header()
        self._build_content()
        self._build_action_buttons()

        self.status_box = ctk.CTkLabel(self, text="", text_color="green")
        self.status_box.grid(row=4, column=0, columnspan=2,
                             sticky="w", padx=10, pady=(0, 10))

        self.search_fields: list = []
        self.search_vars:   list = []
        self._add_search_field(is_first=True)
        self.refresh_variants()

    # ------------------------------------------------------------------
    def _build_header(self):
        hf = ctk.CTkFrame(self)
        hf.grid(row=0, column=0, columnspan=2, sticky="ew", padx=10, pady=20)
        hf.grid_columnconfigure(0, weight=1)
        ctk.CTkLabel(hf, text="Variant Databank", corner_radius=15,
                     fg_color="#4a90e2", text_color="white",
                     font=self.controller.LARGEFONT, height=50
                     ).grid(row=0, column=0, padx=10, pady=10, sticky="ew")
        ctk.CTkButton(hf, text="Back to Input",
                      command=lambda: self.controller.show_frame(InputPage)
                      ).grid(row=0, column=1, sticky="e", padx=10, pady=10)

    def _build_content(self):
        cf = ctk.CTkFrame(self)
        cf.grid(row=1, column=0, columnspan=2, sticky="nsew", padx=10, pady=5)
        cf.grid_columnconfigure(0, weight=1)
        cf.grid_columnconfigure(1, weight=1)
        cf.grid_rowconfigure(0, weight=1)

        # Variant list
        vf = ctk.CTkFrame(cf)
        vf.grid(row=0, column=0, sticky="nsew", padx=(5, 2), pady=5)
        vf.grid_columnconfigure(0, weight=1)
        vf.grid_rowconfigure(1, weight=1)
        ctk.CTkLabel(vf, text="Variants:",
                     font=self.controller.MEDIUMFONT
                     ).grid(row=0, column=0, sticky="w", padx=10, pady=(10, 5))
        self.variant_listbox = ScrollableFrameWithWheel(vf)
        self.variant_listbox.grid(row=1, column=0, sticky="nsew",
                                  padx=10, pady=5)

        # Search
        sf = ctk.CTkFrame(cf)
        sf.grid(row=0, column=1, sticky="nsew", padx=(2, 5), pady=5)
        sf.grid_columnconfigure(0, weight=1)
        sf.grid_rowconfigure(1, weight=1)
        ctk.CTkLabel(sf, text="Search & Filter:",
                     font=self.controller.MEDIUMFONT
                     ).grid(row=0, column=0, sticky="w", padx=10, pady=(10, 5))
        self.search_frame = ScrollableFrameWithWheel(sf)
        self.search_frame.grid(row=1, column=0, sticky="nsew",
                               padx=10, pady=5)
        self.search_frame.grid_columnconfigure(0, weight=1)
        self._build_search_options()

        self.results_label = ctk.CTkLabel(self, text="", text_color="gray")
        self.results_label.grid(row=2, column=0, columnspan=2,
                                sticky="w", padx=10, pady=(5, 0))

    def _build_search_options(self):
        of = ctk.CTkFrame(self.search_frame)
        of.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
        of.grid_columnconfigure(5, weight=1)

        ctk.CTkLabel(of, text="Mutations:",
                     font=self.controller.SMALLFONT
                     ).grid(row=0, column=0, padx=5, pady=5, sticky="w")

        self.min_mutations_var = ctk.StringVar(value="")
        self.max_mutations_var = ctk.StringVar(value="")
        for col, var, ph in [(1, self.min_mutations_var, "Min"),
                             (3, self.max_mutations_var, "Max")]:
            e = ctk.CTkEntry(of, textvariable=var,
                             placeholder_text=ph, width=50)
            e.grid(row=0, column=col, padx=2, pady=5, sticky="w")
            var.trace("w", self._on_search_change)
        ctk.CTkLabel(of, text="to",
                     font=self.controller.SMALLFONT
                     ).grid(row=0, column=2, padx=2, pady=5, sticky="w")
        ctk.CTkButton(of, text="Clear", width=60,
                      command=self._clear_mutation_filter
                      ).grid(row=0, column=4, padx=10, pady=5, sticky="w")

        ctk.CTkLabel(of, text="Logic:",
                     font=self.controller.SMALLFONT
                     ).grid(row=1, column=0, padx=5, pady=5, sticky="w")
        self.search_logic_var = ctk.StringVar(value="AND")
        ctk.CTkOptionMenu(of, values=["AND", "OR"],
                          variable=self.search_logic_var,
                          command=self._on_search_change, width=70
                          ).grid(row=1, column=1, padx=5, pady=5, sticky="w")

        self.add_btn = ctk.CTkButton(of, text="+", width=30, height=30,
                                     command=self._add_search_field)
        self.add_btn.grid(row=1, column=2, padx=2, pady=5, sticky="w")
        self.remove_btn = ctk.CTkButton(of, text="−", width=30, height=30,
                                        command=self._remove_search_field,
                                        state="disabled")
        self.remove_btn.grid(row=1, column=3, padx=2, pady=5, sticky="w")
        ctk.CTkButton(of, text="Clear All", width=80,
                      command=self._clear_all_search
                      ).grid(row=1, column=4, padx=10, pady=5, sticky="w")

    def _build_action_buttons(self):
        bf = ctk.CTkFrame(self)
        bf.grid(row=3, column=0, columnspan=2, sticky="ew", padx=10, pady=10)
        row1 = [("Select All", self.select_all),
                ("Deselect All", self.deselect_all),
                ("Delete Selection", self.delete_variants, "orange"),
                ("Undo Delete", self.undo_delete, "gray")]
        row2 = [("Print Labels", self.print_selected_labels),
                ("Refresh", self.refresh_variants),
                ("View Proteins", lambda: self.controller.show_frame(VariantProteinPage)),
                ("Add Variants", self._add_variants_window),
                ("Edit Variants", self._open_variant_editor)]

        for row_idx, buttons in enumerate([row1, row2]):
            for col, spec in enumerate(buttons):
                txt, cmd = spec[0], spec[1]
                color    = spec[2] if len(spec) > 2 else None
                kwargs: Dict = {"text": txt, "command": cmd, "width": 160,
                                "font": self.controller.MEDIUMFONT}
                if color:
                    kwargs["fg_color"] = color
                ctk.CTkButton(bf, **kwargs).grid(
                    row=row_idx, column=col, sticky="w", padx=10, pady=5)

    # ------------------------------------------------------------------
    # Search helpers
    # ------------------------------------------------------------------

    def _add_search_field(self, is_first=False):
        field_num  = len(self.search_fields)
        search_var = ctk.StringVar()
        search_var.trace("w", self._on_search_change)
        self.search_vars.append(search_var)

        ff = ctk.CTkFrame(self.search_frame)
        ff.grid(row=field_num + 2, column=0, sticky="ew", padx=5, pady=2)
        ff.grid_columnconfigure(1, weight=1)
        ctk.CTkLabel(ff, text="Search:",
                     font=self.controller.SMALLFONT
                     ).grid(row=0, column=0, sticky="w", padx=5, pady=5)
        entry = ctk.CTkEntry(ff, textvariable=search_var,
                             placeholder_text="Search by variant ID or mutations...")
        entry.grid(row=0, column=1, sticky="ew", padx=5, pady=5)

        self.search_fields.append({'frame': ff, 'entry': entry, 'var': search_var})
        self.remove_btn.configure(
            state="normal" if len(self.search_fields) > 1 else "disabled")
        self.filter_variants()

    def _remove_search_field(self):
        if len(self.search_fields) <= 1:
            return
        last = self.search_fields.pop()
        self.search_vars.pop()
        last['frame'].destroy()
        self.remove_btn.configure(
            state="normal" if len(self.search_fields) > 1 else "disabled")
        self.filter_variants()

    def _clear_all_search(self):
        for v in self.search_vars:
            v.set("")

    def _clear_mutation_filter(self):
        self.min_mutations_var.set("")
        self.max_mutations_var.set("")

    def _on_search_change(self, *_):
        self.filter_variants()

    def _get_active_search_terms(self) -> List[str]:
        return [v.get().lower().strip()
                for v in self.search_vars if v.get().strip()]

    def _count_mutations(self, mut_str: str) -> int:
        if not mut_str or mut_str == "(none)":
            return 0
        return len([m for m in mut_str.split(',') if m.strip()])

    def filter_variants(self):
        terms    = self._get_active_search_terms()
        min_str  = self.min_mutations_var.get().strip()
        max_str  = self.max_mutations_var.get().strip()
        min_muts = int(min_str) if min_str.isdigit() else None
        max_muts = int(max_str) if max_str.isdigit() else None

        if not terms and min_muts is None and max_muts is None:
            self.filtered_variants = dict(self.all_variants)
        else:
            logic = self.search_logic_var.get()
            self.filtered_variants = {}
            for vid, muts in self.all_variants.items():
                text  = f"{vid} {muts}".lower()
                count = self._count_mutations(muts)

                text_match = True
                if terms:
                    text_match = (
                        all(t in text for t in terms) if logic == "AND"
                        else any(t in text for t in terms)
                    )
                count_match = True
                if min_muts is not None and count < min_muts:
                    count_match = False
                if max_muts is not None and count > max_muts:
                    count_match = False

                if text_match and count_match:
                    self.filtered_variants[vid] = muts

        self._update_variant_display()

    def _update_variant_display(self):
        for w in self.variant_listbox.winfo_children():
            w.destroy()
        self.check_vars.clear()
        self.checkboxes.clear()

        total    = len(self.all_variants)
        filtered = len(self.filtered_variants)
        terms    = self._get_active_search_terms()
        self.results_label.configure(
            text=(f"Showing all {total} variants" if not terms
                  else f"Showing {filtered} of {total}"))

        for vid, muts in sorted(self.filtered_variants.items(),
                                 key=lambda x: natural_sort_key(x[0])):
            var_text  = f"{vid}: {muts if muts != '(none)' else 'no mutations'}"
            fg_color  = None
            text_color = "white"
            hover_color = None

            if self.controller.is_variant_repair(vid):
                fg_color = "#A5D6A7"; text_color = "black"; hover_color = "#81C784"
            elif self.controller.is_variant_modified(vid):
                fg_color = "#FFF59D"; text_color = "black"; hover_color = "#FFF176"

            cv = ctk.IntVar(value=0)
            if fg_color:
                frame = ctk.CTkFrame(self.variant_listbox,
                                     fg_color=fg_color, corner_radius=5)
                frame.pack(fill="x", padx=10, pady=2, anchor="w")
                cb = ctk.CTkCheckBox(frame, text=var_text, variable=cv,
                                     text_color=text_color,
                                     fg_color=hover_color or fg_color,
                                     hover_color=hover_color or fg_color)
                cb.pack(fill="x", padx=5, pady=2, anchor="w")
            else:
                cb = ctk.CTkCheckBox(self.variant_listbox,
                                     text=var_text, variable=cv)
                cb.pack(fill="x", padx=10, pady=2, anchor="w")

            self.check_vars.append((vid, cv))
            self.checkboxes.append(cb)

    # ------------------------------------------------------------------
    # Actions
    # ------------------------------------------------------------------

    def select_all(self):
        for _, cv in self.check_vars:
            cv.set(1)

    def deselect_all(self):
        for _, cv in self.check_vars:
            cv.set(0)

    def refresh_variants(self):
        self.all_variants = self.controller.get_databank()
        self.filter_variants()

    def delete_variants(self):
        to_del = [vid for vid, cv in self.check_vars if cv.get() == 1]
        if not to_del:
            self.status_box.configure(
                text="No variants selected for deletion.", text_color="red")
            return
        db = self.controller.get_databank()
        self.last_deleted = {vid: db[vid] for vid in to_del if vid in db}
        for vid in to_del:
            db.pop(vid, None)
        self.controller.save_databank(db)
        self.refresh_variants()
        self.status_box.configure(
            text=f"Deleted {len(to_del)} variant(s). You can undo this action.",
            text_color="green")

    def undo_delete(self):
        if not self.last_deleted:
            self.status_box.configure(text="Nothing to undo.", text_color="red")
            return
        db = self.controller.get_databank()
        db.update(self.last_deleted)
        self.controller.save_databank(db)
        cnt = len(self.last_deleted)
        self.last_deleted.clear()
        self.refresh_variants()
        self.status_box.configure(
            text=f"Restored {cnt} variant(s).", text_color="green")

    def print_selected_labels(self):
        selected = [vid for vid, cv in self.check_vars if cv.get() == 1]
        if not selected:
            self.status_box.configure(
                text="Please select variants to print labels for.",
                text_color="red")
            return
        LabelSelectionWindow(self, selected,
                             self.controller.output_dir.get(),
                             self.controller)

    def _add_variants_window(self):
        AddVariantsWindow(self, self.controller.output_dir.get(),
                          self.controller)

    def _open_variant_editor(self):
        VariantEditorWindow(self, self.controller)

    def get_selected_variants(self) -> Dict:
        return {vid: self.all_variants.get(vid, "")
                for vid, cv in self.check_vars if cv.get() == 1}


# ===========================================================================
# AddVariantsWindow
# ===========================================================================

class AddVariantsWindow:
    def __init__(self, parent, output_dir: str, controller):
        self.parent     = parent
        self.output_dir = output_dir
        self.controller = controller
        self.variant_entries: list = []

        self.window = ctk.CTkToplevel(parent)
        self.window.title("Add New Variants to Databank")
        self.window.geometry("900x600")
        self.window.transient(parent)
        self.window.grid_columnconfigure(0, weight=1)
        self.window.grid_rowconfigure(2, weight=1)
        self._setup_ui()
        self.window.after(100, self.window.grab_set)

    def _setup_ui(self):
        ctk.CTkLabel(self.window, text="Add New Variants to Databank",
                     font=self.controller.LARGEFONT
                     ).grid(row=0, column=0, pady=20, padx=20, sticky="ew")

        pf = ctk.CTkFrame(self.window)
        pf.grid(row=1, column=0, sticky="ew", padx=20, pady=10)
        pf.grid_columnconfigure(1, weight=1)
        ctk.CTkLabel(pf, text="Variant Prefix:",
                     font=self.controller.MEDIUMFONT
                     ).grid(row=0, column=0, sticky="w", padx=10, pady=10)
        self.prefix_var = ctk.StringVar(value="KWE_TA_A")
        ctk.CTkEntry(pf, textvariable=self.prefix_var, width=150
                     ).grid(row=0, column=1, sticky="w", padx=10, pady=10)

        # Scrollable variant table
        table_frame = ctk.CTkFrame(self.window)
        table_frame.grid(row=2, column=0, sticky="nsew", padx=20, pady=10)
        table_frame.grid_columnconfigure(1, weight=1)
        table_frame.grid_rowconfigure(1, weight=1)

        ctk.CTkLabel(table_frame, text="Mutations",
                     font=self.controller.MEDIUMFONT
                     ).grid(row=0, column=1, padx=10, pady=5, sticky="w")

        self.scroll = ScrollableFrameWithWheel(table_frame)
        self.scroll.grid(row=1, column=0, columnspan=2, sticky="nsew",
                         padx=10, pady=5)
        self.scroll.grid_columnconfigure(1, weight=1)

        self._add_variant_row()

        bf = ctk.CTkFrame(self.window)
        bf.grid(row=3, column=0, sticky="ew", padx=20, pady=10)

        def _add_row():
            self._add_variant_row()

        def _save():
            self._save_variants()

        for col, (txt, cmd) in enumerate([
            ("+ Add Row",    _add_row),
            ("Save Variants", _save),
            ("Cancel",       self.window.destroy),
        ]):
            ctk.CTkButton(bf, text=txt, command=cmd
                          ).grid(row=0, column=col, padx=10, pady=10)

    def _add_variant_row(self):
        row = len(self.variant_entries)
        ctk.CTkLabel(self.scroll, text=f"Variant {row + 1}:"
                     ).grid(row=row, column=0, padx=(10, 5), pady=5, sticky="w")
        ev = ctk.StringVar()
        ctk.CTkEntry(self.scroll, textvariable=ev,
                     placeholder_text="e.g. K49R, S50A"
                     ).grid(row=row, column=1, padx=(5, 10), pady=5, sticky="ew")
        self.variant_entries.append(ev)

    def _save_variants(self):
        prefix  = self.prefix_var.get().strip()
        db      = self.controller.get_databank()
        pattern = re.compile(
            r"^[ACDEFGHIKLMNPQRSTVWY]\d+[ACDEFGHIKLMNPQRSTVWY]$")

        existing_ids = [
            int(v[len(prefix):])
            for v in db if v.startswith(prefix) and v[len(prefix):].isdigit()
        ]
        next_id = max(existing_ids, default=0) + 1

        added = 0
        for ev in self.variant_entries:
            raw = ev.get().strip()
            if not raw:
                continue
            muts = [m.strip() for m in raw.split(',') if m.strip()]
            if not all(pattern.match(m) for m in muts):
                messagebox.showerror("Validation Error",
                                     f"Invalid mutation format in: '{raw}'")
                return
            normalized = ','.join(muts)
            new_id = f"{prefix}{next_id:02d}"
            db[new_id] = normalized
            next_id += 1
            added += 1

        if added:
            self.controller.save_databank(db)
            self.controller.frames[DatabankPage].refresh_variants()
            messagebox.showinfo("Saved",
                                f"Added {added} variant(s) to databank.")
            self.window.destroy()
        else:
            messagebox.showwarning("Warning", "No variants entered.")


# ===========================================================================
# VariantEditorWindow
# ===========================================================================

class VariantEditorWindow:
    def __init__(self, parent, controller):
        self.parent     = parent
        self.controller = controller
        self.all_variants:      Dict = {}
        self.filtered_variants: Dict = {}
        self.original_databank: Dict = {}
        self.modified_variants: Dict = {}
        self.variant_entries:   Dict = {}
        self.search_fields: list = []
        self.search_vars:   list = []

        self.window = ctk.CTkToplevel(parent)
        self.window.title("Edit Variants")
        self.window.geometry("1000x700")
        self.window.transient(parent)
        self.window.grid_columnconfigure(0, weight=2)
        self.window.grid_columnconfigure(1, weight=1)
        self.window.grid_rowconfigure(1, weight=1)
        self._setup_ui()
        self._load_variants()
        self.window.after(100, self.window.grab_set)

    def _setup_ui(self):
        # Title
        ctk.CTkLabel(self.window, text="Edit Variants",
                     font=self.controller.LARGEFONT
                     ).grid(row=0, column=0, columnspan=2,
                            pady=20, padx=20, sticky="ew")

        # Scrollable variant table
        self.variant_scroll = ScrollableFrameWithWheel(self.window)
        self.variant_scroll.grid(row=1, column=0, sticky="nsew",
                                  padx=(20, 5), pady=10)
        self.variant_scroll.grid_columnconfigure(1, weight=1)

        # Search panel
        sf = ctk.CTkFrame(self.window)
        sf.grid(row=1, column=1, sticky="nsew", padx=(5, 20), pady=10)
        sf.grid_columnconfigure(0, weight=1)
        sf.grid_rowconfigure(1, weight=1)
        ctk.CTkLabel(sf, text="Search:",
                     font=self.controller.MEDIUMFONT
                     ).grid(row=0, column=0, sticky="w", padx=10, pady=5)
        self.search_scroll = ScrollableFrameWithWheel(sf)
        self.search_scroll.grid(row=1, column=0, sticky="nsew",
                                padx=10, pady=5)
        self.search_scroll.grid_columnconfigure(0, weight=1)
        self.search_logic_var = ctk.StringVar(value="AND")
        ctk.CTkOptionMenu(sf, values=["AND", "OR"],
                          variable=self.search_logic_var,
                          command=self._filter_variants, width=80
                          ).grid(row=2, column=0, padx=10, pady=5, sticky="w")
        self.results_label = ctk.CTkLabel(self.window, text="",
                                          text_color="gray")
        self.results_label.grid(row=2, column=0, columnspan=2,
                                sticky="w", padx=20, pady=5)
        self._add_search_field(is_first=True)

        # Action buttons
        bf = ctk.CTkFrame(self.window)
        bf.grid(row=3, column=0, columnspan=2, sticky="ew",
                padx=20, pady=10)
        self.status_label = ctk.CTkLabel(
            bf, text="Modify mutations, then click Save Changes",
            font=self.controller.MEDIUMFONT)
        self.status_label.grid(row=0, column=0, columnspan=3, pady=10, padx=10)
        for col, (txt, cmd) in enumerate([
            ("Save Changes",     self._save_changes),
            ("Cancel",           self.window.destroy),
            ("Reset All Changes", self._reset_changes),
        ]):
            ctk.CTkButton(bf, text=txt, command=cmd
                          ).grid(row=1, column=col, padx=10, pady=10, sticky="w")

    def _add_search_field(self, is_first=False):
        sv = ctk.StringVar()
        sv.trace("w", self._filter_variants)
        self.search_vars.append(sv)
        ff = ctk.CTkFrame(self.search_scroll)
        ff.grid(row=len(self.search_fields), column=0,
                sticky="ew", padx=5, pady=2)
        ff.grid_columnconfigure(0, weight=1)
        ctk.CTkEntry(ff, textvariable=sv,
                     placeholder_text="Search variants..."
                     ).grid(row=0, column=0, sticky="ew", padx=5, pady=5)
        self.search_fields.append({'frame': ff, 'var': sv})

    def _get_active_search_terms(self) -> List[str]:
        return [v.get().lower().strip()
                for v in self.search_vars if v.get().strip()]

    def _filter_variants(self, *_):
        terms = self._get_active_search_terms()
        logic = self.search_logic_var.get()
        if not terms:
            self.filtered_variants = dict(self.all_variants)
        else:
            self.filtered_variants = {}
            for vid, muts in self.all_variants.items():
                text = f"{vid} {muts}".lower()
                match = (all(t in text for t in terms) if logic == "AND"
                         else any(t in text for t in terms))
                if match:
                    self.filtered_variants[vid] = muts
        self._update_variant_display()

    def _load_variants(self):
        self.all_variants      = self.controller.get_databank()
        self.original_databank = copy.deepcopy(self.all_variants)
        self._filter_variants()

    def _update_variant_display(self):
        for w in self.variant_scroll.winfo_children():
            w.destroy()
        self.variant_entries.clear()

        total    = len(self.all_variants)
        filtered = len(self.filtered_variants)
        self.results_label.configure(
            text=(f"Showing all {total} variants" if not self._get_active_search_terms()
                  else f"Showing {filtered} of {total}"))

        if not self.filtered_variants:
            ctk.CTkLabel(self.variant_scroll,
                         text="No variants found.",
                         text_color="orange"
                         ).grid(row=0, column=0, pady=20)
            return

        for i, (vid, muts) in enumerate(
            sorted(self.filtered_variants.items(),
                   key=lambda x: natural_sort_key(x[0]))
        ):
            ctk.CTkLabel(self.variant_scroll, text=f"{vid}:",
                         font=self.controller.MEDIUMFONT
                         ).grid(row=i, column=0, sticky="w",
                                padx=(10, 5), pady=5)
            ev = ctk.StringVar(value=muts)
            ctk.CTkEntry(self.variant_scroll, textvariable=ev
                         ).grid(row=i, column=1, sticky="ew",
                                padx=(5, 10), pady=5)
            ev.trace("w", lambda *_, vid=vid: self._mark_modified(vid))
            self.variant_entries[vid] = ev

    def _mark_modified(self, vid: str):
        if vid in self.variant_entries:
            cur = self.variant_entries[vid].get()
            orig = self.original_databank.get(vid, "")
            if cur != orig:
                self.modified_variants[vid] = cur
            else:
                self.modified_variants.pop(vid, None)

    def _reset_changes(self):
        for vid, ev in self.variant_entries.items():
            if vid in self.original_databank:
                ev.set(self.original_databank[vid])
        self.modified_variants.clear()
        self.status_label.configure(
            text="All changes reset to original values", text_color="orange")

    def _save_changes(self):
        if not self.modified_variants:
            messagebox.showinfo("Info", "No changes to save")
            return
        pattern = re.compile(
            r"^[ACDEFGHIKLMNPQRSTVWY]\d+[ACDEFGHIKLMNPQRSTVWY]$")
        normalized: Dict = {}
        errors: List[str] = []
        for vid, raw in self.modified_variants.items():
            clean = ','.join(m.strip() for m in raw.split(',') if m.strip()) \
                if raw and raw != "(none)" else "(none)"
            muts = [m for m in clean.split(',') if m and m != "(none)"]
            bad  = [m for m in muts if not pattern.match(m)]
            if bad:
                errors.append(f"{vid}: invalid format {bad}")
            else:
                normalized[vid] = clean
        if errors:
            messagebox.showerror("Validation Error",
                                 "\n".join(errors))
            return
        summary = "\n".join(
            f"{vid}:\n  Old: {self.original_databank.get(vid,'')}\n  New: {v}"
            for vid, v in sorted(normalized.items())
        )
        if not messagebox.askyesno(
            "Confirm Changes",
            f"Save changes to {len(normalized)} variant(s)?\n\n{summary}"
        ):
            return
        db = self.controller.get_databank()
        self.controller.undo_stack.append(copy.deepcopy(db))
        db.update(normalized)
        self.controller.save_databank(db)
        self.controller.frames[DatabankPage].refresh_variants()
        messagebox.showinfo("Saved",
                            f"Saved {len(normalized)} variant(s).")
        self.window.destroy()


# ===========================================================================
# ProtocolResultsPage
# ===========================================================================

class ProtocolResultsPage(ctk.CTkFrame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(2, weight=1)

        hf = ctk.CTkFrame(self)
        hf.grid(row=0, column=0, sticky="ew", padx=10, pady=20)
        hf.grid_columnconfigure(0, weight=1)
        ctk.CTkLabel(hf, text="Protocol Results", corner_radius=10,
                     fg_color="#4a90e2", text_color="white",
                     font=controller.LARGEFONT, height=50
                     ).grid(row=0, column=0, padx=10, pady=10, sticky="ew")
        ctk.CTkButton(hf, text="Back to Input",
                      command=lambda: controller.show_frame(InputPage)
                      ).grid(row=0, column=1, sticky="e", padx=10, pady=10)

        cf = ctk.CTkFrame(self)
        cf.grid(row=1, column=0, sticky="ew", padx=10, pady=5)
        ctk.CTkButton(cf, text="Refresh Results",
                      command=self.refresh_results
                      ).grid(row=0, column=0, sticky="w", padx=10, pady=10)

        self.results_text = ctk.CTkTextbox(self, font=controller.SMALLFONT)
        self.results_text.grid(row=2, column=0, sticky="nsew",
                               padx=10, pady=10)

    def refresh_results(self):
        self.results_text.configure(state="normal")
        self.results_text.delete("1.0", "end")

        protocols_dir = Path(self.controller.output_dir.get()) / "protocols"
        data = io.load_protocol_data(protocols_dir)
        if data is None:
            self.results_text.insert(
                "end", "No protocol_data.json found. Run workflow first.\n")
            self.results_text.configure(state="disabled")
            return

        self.results_text.insert("end",
                                 "=== STEP-BY-STEP MUTAGENESIS PROTOCOL ===\n\n")

        final = data.get("final_variants_overview", {})
        if final:
            self.results_text.insert("end", "FINAL VARIANTS OVERVIEW:\n")
            self.results_text.insert("end", "=" * 80 + "\n")
            for label, info in final.items():
                fv   = info.get("final_variant", "")
                muts = ", ".join(info.get("mutations", []))
                self.results_text.insert(
                    "end", f"{label:<15} -> {fv:<15} | {muts}\n")
            self.results_text.insert("end", "\n")

        steps = data.get("protocol_steps", {})
        if steps:
            self.results_text.insert("end", "STEP-BY-STEP PROTOCOL:\n")
            self.results_text.insert("end", "=" * 80 + "\n")
            for step_num in sorted(map(int, steps.keys())):
                variants = steps.get(str(step_num), [])
                self.results_text.insert("end", f"\nSTEP {step_num}:\n")
                self.results_text.insert("end", "-" * 60 + "\n")
                for vi in variants:
                    nv    = vi.get("new_variant", "Unknown")
                    pv    = vi.get("parent_variant", "Unknown")
                    added = vi.get("mutations_added_str", "")
                    self.results_text.insert(
                        "end", f"Create: {nv:<15} from {pv:<15}\n")
                    self.results_text.insert(
                        "end", f"       Add mutations: {added}\n\n")

        all_vars = data.get("all_variants", {})
        self.results_text.insert("end", "\nALL VARIANTS IN PROTOCOL DATA:\n")
        self.results_text.insert("end", "=" * 80 + "\n")
        for v, muts in sorted(all_vars.items(),
                               key=lambda x: natural_sort_key(x[0])):
            self.results_text.insert("end", f"{v:<15} | {muts}\n")
        self.results_text.insert(
            "end", f"\nTotal variants: {len(all_vars)}\n")
        self.results_text.configure(state="disabled")


# ===========================================================================
# PrimerEditorWindow
# ===========================================================================

class PrimerEditorWindow:
    """Per-primer sequence editor with +/- 5' and 3' base buttons."""

    def __init__(self, parent, controller, primer_data: list):
        self.parent      = parent
        self.controller  = controller
        self.primer_data = copy.deepcopy(primer_data)
        self.current_index = 0
        self.modified = False

        self.window = ctk.CTkToplevel(parent)
        self.window.title("Primer Editor")
        self.window.geometry("900x600")
        self.window.transient(parent)
        self.window.grid_columnconfigure(0, weight=1)
        self.window.grid_rowconfigure(2, weight=1)

        # Get WT sequence from InputPage
        try:
            self.wt_sequence = controller.frames[InputPage].get_dna_sequence()
        except Exception:
            self.wt_sequence = ""

        self._setup_ui()
        self.window.after(100, self.window.grab_set)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _is_reverse_primer(name: str) -> bool:
        return '_rev' in name.lower() or '_reverse' in name.lower()

    def _find_primer_pos(self, primer_seq: str, is_reverse: bool) -> int:
        ref = (reverse_complement(self.wt_sequence)
               if is_reverse else self.wt_sequence)
        if not ref or not primer_seq:
            return -1
        pos = ref.find(primer_seq)
        if pos != -1:
            return pos
        # Approximate match (80% threshold)
        best_pos, best_score = -1, 0.0
        for i in range(len(ref) - len(primer_seq) + 1):
            score = sum(a == b for a, b in zip(
                ref[i:i + len(primer_seq)], primer_seq)
            ) / len(primer_seq)
            if score > best_score and score > 0.8:
                best_score, best_pos = score, i
        return best_pos

    def _next_base(self, seq: str, direction: str,
                   is_reverse: bool) -> Optional[str]:
        ref = (reverse_complement(self.wt_sequence)
               if is_reverse else self.wt_sequence)
        if not ref:
            return None
        pos = self._find_primer_pos(seq, is_reverse)
        if pos == -1:
            return None
        if direction == '5prime' and pos > 0:
            return ref[pos - 1]
        if direction == '3prime' and pos + len(seq) < len(ref):
            return ref[pos + len(seq)]
        return None

    def _get_current_sequence(self) -> str:
        raw = self.sequence_text.get("1.0", "end-1c")
        return raw.replace('\n', '').strip()

    def _set_sequence(self, seq: str):
        self.sequence_text.delete("1.0", "end")
        formatted = '\n'.join(seq[i:i+60] for i in range(0, len(seq), 60))
        self.sequence_text.insert("1.0", formatted)

    def _update_primer_info(self):
        seq  = self._get_current_sequence()
        tm   = calculate_tm(seq)
        gc   = calculate_gc_content(seq)
        name = self.primer_data[self.current_index].get("Primer Name", "")

        self.primer_data[self.current_index]["Primer Sequence"] = seq
        self.primer_data[self.current_index]["Length"]          = str(len(seq))
        self.primer_data[self.current_index]["Tm (C)"]          = str(tm)
        self.primer_data[self.current_index]["GC Content (%)"]  = str(gc)

        self.tm_label.configure(
            text=f"{tm} °C" if tm > 0 else "0.0 °C (primer3 unavailable)")
        self.length_label.configure(text=f"{len(seq)} bp")

        is_rev = self._is_reverse_primer(name)
        left_next  = self._next_base(seq, '5prime', is_rev)
        right_next = self._next_base(seq, '3prime', is_rev)
        self.left_next_base_label.configure(
            text=f"Next: {left_next}" if left_next else "Next: —")
        self.right_next_base_label.configure(
            text=f"Next: {right_next}" if right_next else "Next: —")

    # ------------------------------------------------------------------
    # UI
    # ------------------------------------------------------------------

    def _setup_ui(self):
        # Header / navigation
        hf = ctk.CTkFrame(self.window)
        hf.grid(row=0, column=0, sticky="ew", padx=20, pady=20)
        hf.grid_columnconfigure(1, weight=1)
        ctk.CTkLabel(hf, text="Primer Editor",
                     font=self.controller.LARGEFONT
                     ).grid(row=0, column=0, columnspan=3, pady=10)

        nav = ctk.CTkFrame(hf)
        nav.grid(row=1, column=0, columnspan=3, pady=10)
        self.prev_btn = ctk.CTkButton(nav, text="◄ Previous",
                                      command=self._prev_primer, width=100)
        self.prev_btn.grid(row=0, column=0, padx=5)
        self.primer_label = ctk.CTkLabel(nav, text="Primer 1 of X",
                                         font=self.controller.MEDIUMFONT)
        self.primer_label.grid(row=0, column=1, padx=20)
        self.next_btn = ctk.CTkButton(nav, text="Next ►",
                                      command=self._next_primer, width=100)
        self.next_btn.grid(row=0, column=2, padx=5)

        # Info row
        info = ctk.CTkFrame(self.window)
        info.grid(row=1, column=0, sticky="ew", padx=20, pady=10)
        info.grid_columnconfigure(5, weight=1)
        ctk.CTkLabel(info, text="Primer Name:",
                     font=self.controller.MEDIUMFONT
                     ).grid(row=0, column=0, sticky="w", padx=10, pady=5)
        self.name_label = ctk.CTkLabel(info, text="",
                                       font=self.controller.MEDIUMFONT,
                                       text_color="#4a90e2")
        self.name_label.grid(row=0, column=1, sticky="w", padx=10, pady=5)
        ctk.CTkLabel(info, text="Current Tm:",
                     font=self.controller.MEDIUMFONT
                     ).grid(row=0, column=2, sticky="w", padx=10, pady=5)
        self.tm_label = ctk.CTkLabel(info, text="0.0 °C",
                                     font=self.controller.MEDIUMFONT,
                                     text_color="green")
        self.tm_label.grid(row=0, column=3, sticky="w", padx=10, pady=5)
        ctk.CTkLabel(info, text="Length:",
                     font=self.controller.MEDIUMFONT
                     ).grid(row=0, column=4, sticky="w", padx=10, pady=5)
        self.length_label = ctk.CTkLabel(info, text="0 bp",
                                         font=self.controller.MEDIUMFONT)
        self.length_label.grid(row=0, column=5, sticky="w", padx=10, pady=5)

        # Editor (5' | sequence | 3')
        ef = ctk.CTkFrame(self.window)
        ef.grid(row=2, column=0, sticky="nsew", padx=20, pady=10)
        ef.grid_columnconfigure(1, weight=1)
        ef.grid_rowconfigure(0, weight=1)

        for side, col, dir5, dir3, add_cmd, rem_cmd in [
            ("5' End", 0, None, None,
             self._add_left_base, self._remove_left_base),
            ("3' End", 2, None, None,
             self._add_right_base, self._remove_right_base),
        ]:
            sf = ctk.CTkFrame(ef)
            sf.grid(row=0, column=col, sticky="ns", padx=5, pady=5)
            ctk.CTkLabel(sf, text=side,
                         font=self.controller.MEDIUMFONT).pack(pady=10)
            ctk.CTkButton(sf, text="− Remove", width=100, height=40,
                          command=rem_cmd).pack(pady=10)
            if col == 0:
                self.left_next_base_label = ctk.CTkLabel(
                    sf, text="Next: ?",
                    font=ctk.CTkFont(size=16, weight="bold"),
                    text_color="#4a90e2")
                self.left_next_base_label.pack(pady=10)
            else:
                self.right_next_base_label = ctk.CTkLabel(
                    sf, text="Next: ?",
                    font=ctk.CTkFont(size=16, weight="bold"),
                    text_color="#4a90e2")
                self.right_next_base_label.pack(pady=10)
            ctk.CTkButton(sf, text="+ Add", width=100, height=40,
                          command=add_cmd).pack(pady=10)

        seq_frame = ctk.CTkFrame(ef)
        seq_frame.grid(row=0, column=1, sticky="nsew", padx=5, pady=5)
        seq_frame.grid_rowconfigure(1, weight=1)
        seq_frame.grid_columnconfigure(0, weight=1)
        ctk.CTkLabel(seq_frame, text="Primer Sequence",
                     font=self.controller.MEDIUMFONT
                     ).grid(row=0, column=0, pady=10)
        self.sequence_text = ctk.CTkTextbox(
            seq_frame, font=ctk.CTkFont(family="Courier", size=14), wrap="char")
        self.sequence_text.grid(row=1, column=0, sticky="nsew",
                                padx=10, pady=10)

        # Buttons
        bf = ctk.CTkFrame(self.window)
        bf.grid(row=3, column=0, sticky="ew", padx=20, pady=10)
        self.status_label = ctk.CTkLabel(
            bf, text="Make changes and click Apply",
            font=self.controller.MEDIUMFONT)
        self.status_label.grid(row=0, column=0, columnspan=3, pady=10)
        for col, (txt, cmd, kwargs) in enumerate([
            ("Reset Current",  self._reset_current_primer, {}),
            ("Apply Changes",  self._apply_changes,
             {"fg_color": "green"}),
            ("Close",          self._close_editor,        {}),
        ]):
            ctk.CTkButton(bf, text=txt, command=cmd, **kwargs
                          ).grid(row=1, column=col, padx=10, pady=10)

        self.window.update_idletasks()
        self._load_primer(0)

    # ------------------------------------------------------------------
    # Navigation
    # ------------------------------------------------------------------

    def _load_primer(self, index: int):
        if not (0 <= index < len(self.primer_data)):
            return
        self.current_index = index
        primer = self.primer_data[index]
        name   = primer.get("Primer Name", "")
        seq    = primer.get("Primer Sequence", "")
        is_rev = self._is_reverse_primer(name)

        self.name_label.configure(text=name)
        suffix = " (Reverse)" if is_rev else " (Forward)"
        self.primer_label.configure(
            text=f"Primer {index+1} of {len(self.primer_data)}{suffix}")
        self._set_sequence(seq)
        self._update_primer_info()

        self.prev_btn.configure(state="normal" if index > 0 else "disabled")
        self.next_btn.configure(
            state="normal" if index < len(self.primer_data)-1 else "disabled")

    def _prev_primer(self):
        if self.current_index > 0:
            self._load_primer(self.current_index - 1)

    def _next_primer(self):
        if self.current_index < len(self.primer_data) - 1:
            self._load_primer(self.current_index + 1)

    # ------------------------------------------------------------------
    # Edit actions
    # ------------------------------------------------------------------

    def _add_left_base(self):
        seq     = self._get_current_sequence()
        is_rev  = self._is_reverse_primer(
            self.primer_data[self.current_index].get("Primer Name", ""))
        base    = self._next_base(seq, '5prime', is_rev)
        if not base:
            self.status_label.configure(
                text="Cannot determine next base from WT", text_color="red")
            return
        self._set_sequence(base + seq)
        self._update_primer_info()
        self.modified = True
        self.status_label.configure(
            text=f"Added {base} to 5' end", text_color="green")

    def _remove_left_base(self):
        seq = self._get_current_sequence()
        if len(seq) <= 1:
            self.status_label.configure(
                text="Cannot remove — sequence too short!", text_color="red")
            return
        removed = seq[0]
        self._set_sequence(seq[1:])
        self._update_primer_info()
        self.modified = True
        self.status_label.configure(
            text=f"Removed {removed} from 5' end", text_color="orange")

    def _add_right_base(self):
        seq    = self._get_current_sequence()
        is_rev = self._is_reverse_primer(
            self.primer_data[self.current_index].get("Primer Name", ""))
        base   = self._next_base(seq, '3prime', is_rev)
        if not base:
            self.status_label.configure(
                text="Cannot determine next base from WT", text_color="red")
            return
        self._set_sequence(seq + base)
        self._update_primer_info()
        self.modified = True
        self.status_label.configure(
            text=f"Added {base} to 3' end", text_color="green")

    def _remove_right_base(self):
        seq = self._get_current_sequence()
        if len(seq) <= 1:
            self.status_label.configure(
                text="Cannot remove — sequence too short!", text_color="red")
            return
        removed = seq[-1]
        self._set_sequence(seq[:-1])
        self._update_primer_info()
        self.modified = True
        self.status_label.configure(
            text=f"Removed {removed} from 3' end", text_color="orange")

    def _reset_current_primer(self):
        orig = self.parent.primer_data[self.current_index]
        self.primer_data[self.current_index] = copy.deepcopy(orig)
        self._load_primer(self.current_index)
        self.status_label.configure(
            text="Primer reset to original", text_color="blue")

    def _apply_changes(self):
        if not self.modified:
            self.status_label.configure(
                text="No modifications made", text_color="orange")
            return
        try:
            out_dir = Path(self.controller.output_dir.get())
            io.save_primer_temp_csv(self.primer_data, out_dir)
            self.status_label.configure(
                text="Saved to primer_temp.csv — click 'Save Changes' "
                     "in the Primer page to apply permanently",
                text_color="green")
        except Exception as e:
            self.status_label.configure(
                text=f"Error: {e}", text_color="red")

    def _close_editor(self):
        self.window.destroy()


# ===========================================================================
# PrimerPage
# ===========================================================================

class PrimerPage(ctk.CTkFrame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller  = controller
        self.table_labels: list = []
        self.primer_data:  list = []
        self.wt_sequence    = ""
        self.wt_sequence_rc = ""

        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(2, weight=1)

        # Header
        hf = ctk.CTkFrame(self)
        hf.grid(row=0, column=0, sticky="ew", padx=10, pady=20)
        hf.grid_columnconfigure(0, weight=1)
        ctk.CTkLabel(hf, text="Primer List", corner_radius=10,
                     fg_color="#4a90e2", text_color="white",
                     font=controller.LARGEFONT, height=50
                     ).grid(row=0, column=0, padx=10, pady=10, sticky="ew")
        ctk.CTkButton(hf, text="Back to Input",
                      command=lambda: controller.show_frame(InputPage)
                      ).grid(row=0, column=1, sticky="e", padx=10, pady=10)

        # Controls
        cf = ctk.CTkFrame(self)
        cf.grid(row=1, column=0, sticky="ew", padx=10, pady=5)
        for col, (txt, cmd, kwargs) in enumerate([
            ("Copy Name + Sequence", self._copy_selected,      {}),
            ("Refresh",              self.load_primers,         {}),
            ("Edit Primers",         self._open_primer_editor,  {}),
            ("Save Changes",         self._save_primers_to_disk,
             {"fg_color": "green"}),
        ]):
            ctk.CTkButton(cf, text=txt, command=cmd, **kwargs
                          ).grid(row=0, column=col, sticky="w",
                                 padx=10, pady=10)

        # Table
        self.scroll_frame = ScrollableFrameWithWheel(self)
        self.scroll_frame.grid(row=2, column=0, sticky="nsew",
                               padx=10, pady=10)

    # ------------------------------------------------------------------
    # Sequence helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _is_reverse_primer(name: str) -> bool:
        nl = name.lower()
        return '_rev' in nl or '_reverse' in nl or nl.endswith('rev')

    def _get_mutated_positions(self, primer_name: str,
                               primer_seq: str) -> set:
        ref = (self.wt_sequence_rc
               if self._is_reverse_primer(primer_name)
               else self.wt_sequence)
        if not ref or not primer_seq:
            return set()
        pos = ref.find(primer_seq)
        if pos == -1:
            best_pos, best = -1, 0.0
            for i in range(len(ref) - len(primer_seq) + 1):
                score = sum(a == b for a, b in zip(
                    ref[i:i+len(primer_seq)], primer_seq)) / len(primer_seq)
                if score > best and score > 0.75:
                    best, best_pos = score, i
            pos = best_pos
        if pos == -1:
            return set()
        wt_region = ref[pos:pos + len(primer_seq)]
        return {i for i, (pb, wb) in enumerate(zip(primer_seq, wt_region))
                if pb != wb}

    def _create_highlighted_textbox(self, parent, seq: str,
                                    primer_name: str) -> ctk.CTkTextbox:
        tb = ctk.CTkTextbox(parent, height=30, width=400,
                            font=ctk.CTkFont(family="Courier", size=12))
        tb.grid_propagate(False)
        tb.tag_config("mutation", foreground="red")
        mutated = self._get_mutated_positions(primer_name, seq)
        for i, base in enumerate(seq):
            if i in mutated:
                tb.insert("end", base, "mutation")
            else:
                tb.insert("end", base)
        tb.configure(state="disabled")
        return tb

    # ------------------------------------------------------------------
    # Load
    # ------------------------------------------------------------------

    def load_primers(self):
        for w in self.scroll_frame.winfo_children():
            w.destroy()
        self.table_labels.clear()

        out_dir  = Path(self.controller.output_dir.get())
        csv_path = io.active_primer_csv_path(out_dir)

        if csv_path is None:
            ctk.CTkLabel(self.scroll_frame,
                         text="❌ No primer_list.csv found.",
                         text_color="red").pack(pady=20)
            return

        if csv_path.name == io.PRIMER_TEMP_CSV:
            ctk.CTkLabel(self.scroll_frame,
                         text="⚠️ Showing unsaved changes from primer_temp.csv "
                              "— click 'Save Changes' to apply permanently",
                         text_color="orange"
                         ).grid(row=0, column=0, columnspan=6,
                                pady=5, sticky="ew")

        try:
            try:
                ip_frame = self.controller.frames[InputPage]
                self.wt_sequence = ip_frame.get_dna_sequence()
            except Exception:
                self.wt_sequence = ""

            if self.wt_sequence:
                self.wt_sequence_rc = reverse_complement(self.wt_sequence)
            else:
                ctk.CTkLabel(self.scroll_frame,
                             text="⚠️ No wildtype sequence — mutations won't "
                                  "be highlighted",
                             text_color="orange"
                             ).grid(row=1, column=0, columnspan=6,
                                    pady=10, sticky="ew")

            primers = io.load_primer_csv(csv_path)
            if not primers:
                ctk.CTkLabel(self.scroll_frame,
                             text="No primers found in file.",
                             text_color="orange").pack(pady=20)
                return

            self.primer_data = primers

            headers = ["Primer Name", "Primer Sequence",
                       "Tm (C)", "GC%", "Overlap Tm (C)", "Overlap GC%"]
            for col, minsize in [(0, 180), (1, 420), (2, 80),
                                  (3, 80),  (4, 100), (5, 100)]:
                self.scroll_frame.grid_columnconfigure(
                    col, weight=(1 if col == 1 else 0), minsize=minsize)

            def _gc_color(gc_str: str) -> str:
                try:
                    v = float(gc_str)
                    if 40 <= v <= 60: return "#1A531A"
                    if 60 < v <= 70:  return "#C0A50B"
                    return "#C54343"
                except ValueError:
                    return "#1A531A"

            header_row = 0 if not self.wt_sequence else 1
            for i, h in enumerate(headers):
                ctk.CTkLabel(
                    self.scroll_frame, text=h,
                    font=ctk.CTkFont(family="Verdana", size=14,
                                     weight="bold")
                ).grid(row=header_row, column=i, padx=5, pady=10, sticky="w")

            # Add primer3 warning column if needed
            if not _PRIMER3_OK:
                ctk.CTkLabel(
                    self.scroll_frame,
                    text="⚠ Tm = 0 (primer3 unavailable)",
                    text_color="orange",
                    font=ctk.CTkFont(size=10),
                ).grid(row=header_row, column=2, padx=2, pady=2, sticky="w")

            start_row = header_row + 1
            for r, primer in enumerate(primers, start=start_row):
                row_widgets: list = []
                name  = primer.get("Primer Name", "")
                seq   = primer.get("Primer Sequence", "")
                gc    = primer.get("GC Content (%)", "")
                ovlgc = primer.get("Overlap GC (%)", "")

                # Col 0: Name
                nl = ctk.CTkLabel(self.scroll_frame, text=name,
                                  font=ctk.CTkFont(family="Verdana", size=12))
                nl.grid(row=r, column=0, padx=5, pady=2, sticky="w")
                row_widgets.append(nl)

                # Col 1: Sequence (highlighted if WT available)
                if self.wt_sequence:
                    tb = self._create_highlighted_textbox(
                        self.scroll_frame, seq, name)
                    tb.grid(row=r, column=1, padx=5, pady=2, sticky="ew")
                    row_widgets.append(tb)
                else:
                    sl = ctk.CTkLabel(self.scroll_frame, text=seq,
                                      font=ctk.CTkFont(family="Courier",
                                                       size=12))
                    sl.grid(row=r, column=1, padx=5, pady=2, sticky="w")
                    row_widgets.append(sl)

                # Col 2: Tm
                tm_val = primer.get("Tm (C)", "")
                tm_lbl = ctk.CTkLabel(self.scroll_frame, text=tm_val,
                                      font=ctk.CTkFont(size=12))
                tm_lbl.grid(row=r, column=2, padx=5, pady=2, sticky="w")
                row_widgets.append(tm_lbl)

                # Col 3: GC%
                gc_lbl = ctk.CTkLabel(self.scroll_frame, text=gc,
                                      fg_color=_gc_color(gc),
                                      corner_radius=5, width=60)
                gc_lbl.grid(row=r, column=3, padx=5, pady=2, sticky="w")
                row_widgets.append(gc_lbl)

                # Col 4: Overlap Tm
                ovl_tm = primer.get("Overlap Tm", "")
                ctk.CTkLabel(self.scroll_frame, text=ovl_tm,
                             font=ctk.CTkFont(size=12)
                             ).grid(row=r, column=4, padx=5, pady=2,
                                    sticky="w")

                # Col 5: Overlap GC%
                ovl_gc_lbl = ctk.CTkLabel(self.scroll_frame, text=ovlgc,
                                          fg_color=_gc_color(ovlgc),
                                          corner_radius=5, width=60)
                ovl_gc_lbl.grid(row=r, column=5, padx=5, pady=2, sticky="w")
                row_widgets.append(ovl_gc_lbl)

                self.table_labels.append(row_widgets)

        except Exception as e:
            ctk.CTkLabel(self.scroll_frame,
                         text=f"Error loading primers: {e}",
                         text_color="red").pack(pady=20)

    # ------------------------------------------------------------------
    # Actions
    # ------------------------------------------------------------------

    def _copy_selected(self):
        if not self.table_labels:
            return
        lines = ["Primer Name\tPrimer Sequence"]
        for row in self.table_labels:
            name = row[0].cget("text")
            seq_w = row[1]
            seq = (seq_w.get("1.0", "end-1c").replace('\n', '')
                   if isinstance(seq_w, ctk.CTkTextbox)
                   else seq_w.cget("text"))
            lines.append(f"{name}\t{seq}")
        self.clipboard_clear()
        self.clipboard_append("\n".join(lines))

    def _open_primer_editor(self):
        if not self.primer_data:
            messagebox.showwarning("No Primers",
                                   "Please load primers first (Refresh).")
            return
        PrimerEditorWindow(self, self.controller, self.primer_data)

    def _save_primers_to_disk(self):
        out_dir = Path(self.controller.output_dir.get())
        result  = io.promote_temp_to_main(out_dir)
        if result is None:
            messagebox.showinfo("No Changes",
                                "No pending changes (primer_temp.csv not found).")
            return
        self.load_primers()
        messagebox.showinfo("Saved",
                            f"Changes saved to primer_list.csv\n"
                            f"Backup: primer_list.csv.backup")


# ===========================================================================
# MutationExtractorPage
# ===========================================================================

class MutationExtractorPage(ctk.CTkFrame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(3, weight=1)

        hf = ctk.CTkFrame(self)
        hf.grid(row=0, column=0, sticky="ew", padx=10, pady=20)
        hf.grid_columnconfigure(0, weight=1)
        ctk.CTkLabel(hf, text="Mutation Extractor", corner_radius=10,
                     fg_color="#4a90e2", text_color="white",
                     font=controller.LARGEFONT, height=50
                     ).grid(row=0, column=0, padx=10, pady=10, sticky="ew")
        ctk.CTkButton(hf, text="Back to Input",
                      command=lambda: controller.show_frame(InputPage)
                      ).grid(row=0, column=1, padx=10, pady=10, sticky="e")

        ref_f = ctk.CTkFrame(self)
        ref_f.grid(row=1, column=0, sticky="ew", padx=10, pady=10)
        ref_f.grid_columnconfigure(1, weight=1)
        ctk.CTkLabel(ref_f, text="Reference Sequence:",
                     font=controller.MEDIUMFONT
                     ).grid(row=0, column=0, sticky="w", padx=10, pady=10)
        self.ref_sequence_var = ctk.StringVar(value="")
        self.ref_dropdown = ctk.CTkOptionMenu(
            ref_f, values=["No sequences available"],
            variable=self.ref_sequence_var, width=200)
        self.ref_dropdown.grid(row=0, column=1, sticky="w", padx=5, pady=10)
        ctk.CTkButton(ref_f, text="Refresh",
                      command=self.refresh_reference_sequences, width=80
                      ).grid(row=0, column=2, padx=10, pady=10, sticky="w")

        fasta_f = ctk.CTkFrame(self)
        fasta_f.grid(row=2, column=0, sticky="ew", padx=10, pady=10)
        fasta_f.grid_columnconfigure(1, weight=1)
        ctk.CTkLabel(fasta_f, text="FASTA File:",
                     font=controller.MEDIUMFONT
                     ).grid(row=0, column=0, sticky="w", padx=10, pady=10)
        self.fasta_path_var = ctk.StringVar(value="")
        ctk.CTkEntry(fasta_f, textvariable=self.fasta_path_var,
                     state="readonly"
                     ).grid(row=0, column=1, sticky="ew", padx=5, pady=10)
        ctk.CTkButton(fasta_f, text="Browse...",
                      command=self._browse_fasta, width=80
                      ).grid(row=0, column=2, padx=10, pady=10, sticky="w")

        self.fasta_preview = ctk.CTkTextbox(
            self, font=ctk.CTkFont(family="Courier", size=11))
        self.fasta_preview.grid(row=3, column=0, sticky="nsew",
                                padx=10, pady=(0, 10))
        self.fasta_preview.configure(state="disabled")

        af = ctk.CTkFrame(self)
        af.grid(row=4, column=0, sticky="ew", padx=10, pady=10)
        ctk.CTkButton(af, text="Extract Mutations",
                      command=self._extract_mutations,
                      font=controller.MEDIUMFONT
                      ).grid(row=0, column=0, padx=10, pady=10, sticky="w")

        self.status_label = ctk.CTkLabel(
            self,
            text="Select reference sequence and FASTA file to begin",
            font=controller.SMALLFONT, text_color="gray")
        self.status_label.grid(row=5, column=0, sticky="w", padx=10, pady=5)

        self.refresh_reference_sequences()

    def refresh_reference_sequences(self):
        seqs = io.load_wt_sequences(
            Path(self.controller.output_dir.get()))
        if seqs:
            names = list(seqs.keys())
            self.ref_dropdown.configure(values=names)
            self.ref_sequence_var.set(names[0])
            self.status_label.configure(
                text=f"Loaded {len(names)} reference sequence(s)",
                text_color="green")
        else:
            self.ref_dropdown.configure(values=["No sequences available"])
            self.ref_sequence_var.set("")
            self.status_label.configure(
                text="No saved sequences. Save sequences in Input tab first.",
                text_color="orange")

    def _browse_fasta(self):
        path = filedialog.askopenfilename(
            title="Select FASTA File",
            filetypes=[("FASTA files", "*.fasta *.fas *.fa"),
                       ("All files", "*.*")])
        if path:
            self.fasta_path_var.set(path)
            self._preview_fasta(path)

    def _preview_fasta(self, path: str):
        try:
            with open(path, 'r') as f:
                content = f.read()
            preview = content[:1000] + (
                "\n\n... (truncated)" if len(content) > 1000 else "")
            self.fasta_preview.configure(state="normal")
            self.fasta_preview.delete("1.0", "end")
            self.fasta_preview.insert("1.0", preview)
            self.fasta_preview.configure(state="disabled")
            self.status_label.configure(
                text=f"FASTA loaded ({content.count('>')} sequences)",
                text_color="green")
        except Exception as e:
            self.status_label.configure(
                text=f"Error reading FASTA: {e}", text_color="red")

    def _extract_mutations(self):
        ref_name = self.ref_sequence_var.get()
        fasta_path = self.fasta_path_var.get()

        if not ref_name or ref_name == "No sequences available":
            messagebox.showerror("Error",
                                 "Please select a reference sequence.")
            return
        if not fasta_path:
            messagebox.showerror("Error", "Please select a FASTA file.")
            return

        try:
            seqs = io.load_wt_sequences(
                Path(self.controller.output_dir.get()))
            dna  = seqs.get(ref_name, "")
            if not dna or len(dna) <= 60:
                raise ValueError("Reference DNA too short.")
            ref_protein = translate_dna_to_protein(dna[30:-30])

            fasta_seqs = io.parse_fasta(Path(fasta_path))
            if not fasta_seqs:
                raise ValueError("No sequences found in FASTA file.")

            self.status_label.configure(
                text=f"Aligning {len(fasta_seqs)} sequences...",
                text_color="blue")
            self.update_idletasks()

            mutation_lines, errors = [], []
            for header, seq in fasta_seqs:
                try:
                    muts, _, _ = find_mutations_from_alignment(ref_protein, seq)
                    mutation_lines.append(", ".join(muts) if muts else "")
                except Exception as e:
                    errors.append(f"{header}: {e}")

            if errors:
                messagebox.showerror("Alignment Errors", "\n".join(errors))
                if not mutation_lines:
                    return

            non_empty = [l for l in mutation_lines if l.strip()]
            if non_empty:
                ip = self.controller.frames[InputPage]
                ip.var_text.delete("1.0", "end")
                ip.var_text.insert("1.0", "\n".join(non_empty))
                self.status_label.configure(
                    text=f"Extracted mutations from {len(non_empty)} "
                         "variant(s) — transferred to Input tab",
                    text_color="green")
                self.controller.show_frame(InputPage)
            else:
                self.status_label.configure(
                    text="No mutations found — all sequences identical to ref",
                    text_color="orange")

        except Exception as e:
            self.status_label.configure(
                text=f"Error: {e}", text_color="red")


# ===========================================================================
# MutationFailureHandler (Mutation Eraser page)
# ===========================================================================

class MutationFailureHandler(ctk.CTkFrame):
    """Mark failed mutations per variant and produce repair protocols."""

    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller          = controller
        self.variant_rows:       Dict = {}
        self.all_variants:       Dict = {}
        self.filtered_variants:  Dict = {}
        self.search_fields:      list = []
        self.search_vars:        list = []
        self.failure_undo_stack: list = []
        self.min_mutations_var   = ctk.StringVar(value="")
        self.max_mutations_var   = ctk.StringVar(value="")
        self.search_logic_var    = ctk.StringVar(value="AND")

        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)

        self._setup_ui()
        self.add_search_field(is_first=True)
        self.load_variants()

    # ------------------------------------------------------------------
    # Layout
    # ------------------------------------------------------------------

    def _setup_ui(self):
        # Header
        hf = ctk.CTkFrame(self)
        hf.grid(row=0, column=0, sticky="ew", padx=10, pady=20)
        hf.grid_columnconfigure(0, weight=1)
        ctk.CTkLabel(hf, text="Mutation Failure Handler",
                     corner_radius=10, fg_color="#4a90e2",
                     text_color="white", font=self.controller.LARGEFONT,
                     height=50
                     ).grid(row=0, column=0, padx=10, pady=10, sticky="ew")
        ctk.CTkButton(hf, text="Back to Databank",
                      command=lambda: self.controller.show_frame(DatabankPage)
                      ).grid(row=0, column=1, padx=10, pady=10, sticky="e")

        # Main content: variant table left, search right
        content = ctk.CTkFrame(self)
        content.grid(row=1, column=0, sticky="nsew", padx=10, pady=5)
        content.grid_columnconfigure(0, weight=2)
        content.grid_columnconfigure(1, weight=1)
        content.grid_rowconfigure(0, weight=1)

        self.create_variant_table(content)
        self.create_search_section(content)

        # Results label
        self.results_label = ctk.CTkLabel(self, text="", text_color="gray")
        self.results_label.grid(row=2, column=0, sticky="w", padx=10, pady=(5, 0))

        # Action buttons
        bf = ctk.CTkFrame(self)
        bf.grid(row=3, column=0, sticky="ew", padx=10, pady=10)

        ctk.CTkButton(bf, text="Validate Selection",
                      command=self.validate_selection,
                      font=self.controller.MEDIUMFONT,
                      fg_color="orange", width=150
                      ).grid(row=0, column=0, padx=10, pady=5, sticky="w")
        ctk.CTkButton(bf, text="Process Failed Mutations",
                      command=self.process_failures,
                      font=self.controller.MEDIUMFONT,
                      fg_color="green", width=200
                      ).grid(row=0, column=1, padx=10, pady=5, sticky="w")
        ctk.CTkButton(bf, text="Clear All Checks",
                      command=self.clear_all_checks,
                      font=self.controller.MEDIUMFONT, width=120
                      ).grid(row=0, column=2, padx=10, pady=5, sticky="w")
        ctk.CTkButton(bf, text="Undo Last Processing",
                      command=self.undo_last_processing,
                      font=self.controller.MEDIUMFONT,
                      fg_color="gray", width=150
                      ).grid(row=1, column=0, padx=10, pady=5, sticky="w")
        ctk.CTkButton(bf, text="Refresh Variants",
                      command=self.load_variants,
                      font=self.controller.MEDIUMFONT, width=120
                      ).grid(row=1, column=1, padx=10, pady=5, sticky="w")

        self.status_label = ctk.CTkLabel(
            bf, text="", font=self.controller.MEDIUMFONT)
        self.status_label.grid(row=2, column=0, columnspan=3,
                               sticky="w", padx=10, pady=5)

    def create_variant_table(self, parent):
        """Left panel — canvas with horizontal + vertical scrollbars."""
        variant_frame = ctk.CTkFrame(parent)
        variant_frame.grid(row=0, column=0, sticky="nsew", padx=(5, 2), pady=5)
        variant_frame.grid_columnconfigure(0, weight=1)
        variant_frame.grid_rowconfigure(1, weight=1)

        ctk.CTkLabel(variant_frame,
                     text="✓ Check mutations that FAILED (missing)",
                     font=self.controller.SMALLFONT, text_color="orange"
                     ).grid(row=0, column=0, padx=10, pady=5, sticky="w")

        # Container for canvas + two scrollbars
        scroll_container = ctk.CTkFrame(variant_frame)
        scroll_container.grid(row=1, column=0, sticky="nsew", padx=10, pady=5)
        scroll_container.grid_columnconfigure(0, weight=1)
        scroll_container.grid_rowconfigure(0, weight=1)

        import tkinter as tk
        bg = self._apply_appearance_mode(
            ctk.ThemeManager.theme["CTkFrame"]["fg_color"])
        canvas = tk.Canvas(scroll_container, bg=bg, highlightthickness=0)
        canvas.grid(row=0, column=0, sticky="nsew")

        v_scrollbar = ctk.CTkScrollbar(scroll_container,
                                       orientation="vertical",
                                       command=canvas.yview)
        v_scrollbar.grid(row=0, column=1, sticky="ns")

        h_scrollbar = ctk.CTkScrollbar(scroll_container,
                                       orientation="horizontal",
                                       command=canvas.xview)
        h_scrollbar.grid(row=1, column=0, sticky="ew")

        canvas.configure(yscrollcommand=v_scrollbar.set,
                         xscrollcommand=h_scrollbar.set)

        # Inner scrollable frame
        self.scroll_frame = ctk.CTkFrame(canvas)
        canvas_window = canvas.create_window(
            (0, 0), window=self.scroll_frame, anchor="nw")

        def _update_scrollregion(event=None):
            canvas.configure(scrollregion=canvas.bbox("all"))

        self.scroll_frame.bind("<Configure>", _update_scrollregion)

        def _resize_canvas_window(event):
            min_w = event.width
            canvas.itemconfig(canvas_window,
                              width=max(min_w,
                                        self.scroll_frame.winfo_reqwidth()))

        canvas.bind("<Configure>", _resize_canvas_window)

        # Mouse-wheel bindings
        def _on_wheel(event):
            canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        def _on_up(event):   canvas.yview_scroll(-1, "units")
        def _on_down(event): canvas.yview_scroll(1,  "units")

        canvas.bind("<Enter>", lambda e: [
            canvas.bind_all("<MouseWheel>", _on_wheel),
            canvas.bind_all("<Button-4>",   _on_up),
            canvas.bind_all("<Button-5>",   _on_down),
        ])
        canvas.bind("<Leave>", lambda e: [
            canvas.unbind_all("<MouseWheel>"),
            canvas.unbind_all("<Button-4>"),
            canvas.unbind_all("<Button-5>"),
        ])

        self._canvas     = canvas
        self._h_scrollbar = h_scrollbar
        self._v_scrollbar = v_scrollbar

    def create_search_section(self, parent):
        """Right panel — search & filter."""
        sf = ctk.CTkFrame(parent)
        sf.grid(row=0, column=1, sticky="nsew", padx=(2, 5), pady=5)
        sf.grid_columnconfigure(0, weight=1)
        sf.grid_rowconfigure(1, weight=1)

        ctk.CTkLabel(sf, text="Search & Filter:",
                     font=self.controller.MEDIUMFONT
                     ).grid(row=0, column=0, sticky="w", padx=10, pady=(10, 5))

        self.search_scroll = ScrollableFrameWithWheel(sf)
        self.search_scroll.grid(row=1, column=0, sticky="nsew",
                                padx=10, pady=5)
        self.search_scroll.grid_columnconfigure(0, weight=1)

        self.create_search_options()

    def create_search_options(self):
        of = ctk.CTkFrame(self.search_scroll)
        of.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
        of.grid_columnconfigure(5, weight=1)

        ctk.CTkLabel(of, text="Mutations:",
                     font=self.controller.SMALLFONT
                     ).grid(row=0, column=0, padx=5, pady=5, sticky="w")

        min_e = ctk.CTkEntry(of, textvariable=self.min_mutations_var,
                             placeholder_text="Min", width=50)
        min_e.grid(row=0, column=1, padx=2, pady=5, sticky="w")
        self.min_mutations_var.trace("w", self.on_search_change)

        ctk.CTkLabel(of, text="to",
                     font=self.controller.SMALLFONT
                     ).grid(row=0, column=2, padx=2, pady=5, sticky="w")

        max_e = ctk.CTkEntry(of, textvariable=self.max_mutations_var,
                             placeholder_text="Max", width=50)
        max_e.grid(row=0, column=3, padx=2, pady=5, sticky="w")
        self.max_mutations_var.trace("w", self.on_search_change)

        ctk.CTkButton(of, text="Clear",
                      command=self.clear_mutation_filter, width=60
                      ).grid(row=0, column=4, padx=10, pady=5, sticky="w")

        ctk.CTkLabel(of, text="Logic:",
                     font=self.controller.SMALLFONT
                     ).grid(row=1, column=0, padx=5, pady=5, sticky="w")

        ctk.CTkOptionMenu(of, values=["AND", "OR"],
                          variable=self.search_logic_var,
                          command=self.on_search_change, width=70
                          ).grid(row=1, column=1, padx=5, pady=5, sticky="w")

        self.add_btn = ctk.CTkButton(of, text="+", width=30, height=30,
                                     command=self.add_search_field)
        self.add_btn.grid(row=1, column=2, padx=2, pady=5, sticky="w")

        self.remove_btn = ctk.CTkButton(of, text="−", width=30, height=30,
                                        command=self.remove_search_field,
                                        state="disabled")
        self.remove_btn.grid(row=1, column=3, padx=2, pady=5, sticky="w")

        ctk.CTkButton(of, text="Clear All",
                      command=self.clear_all_search, width=80
                      ).grid(row=1, column=4, padx=10, pady=5, sticky="w")

        ctk.CTkLabel(of,
                     text="Use AND to match all terms, OR for any.\n"
                          "Filter by mutation count.",
                     font=self.controller.SMALLFONT, text_color="gray"
                     ).grid(row=2, column=0, columnspan=6,
                            padx=5, pady=(0, 5), sticky="w")

    # ------------------------------------------------------------------
    # Search helpers
    # ------------------------------------------------------------------

    def add_search_field(self, is_first=False):
        sv = ctk.StringVar()
        sv.trace("w", self.on_search_change)
        self.search_vars.append(sv)

        ff = ctk.CTkFrame(self.search_scroll)
        ff.grid(row=len(self.search_fields) + 2, column=0,
                sticky="ew", padx=5, pady=2)
        ff.grid_columnconfigure(1, weight=1)

        ctk.CTkLabel(ff, text="Search:",
                     font=self.controller.SMALLFONT
                     ).grid(row=0, column=0, sticky="w", padx=5, pady=5)
        ctk.CTkEntry(ff, textvariable=sv,
                     placeholder_text="Search by variant ID or mutations..."
                     ).grid(row=0, column=1, sticky="ew", padx=5, pady=5)

        self.search_fields.append({'frame': ff, 'var': sv})
        self.remove_btn.configure(
            state="normal" if len(self.search_fields) > 1 else "disabled")
        self.filter_variants()

    def remove_search_field(self):
        if len(self.search_fields) <= 1:
            return
        last = self.search_fields.pop()
        self.search_vars.pop()
        last['frame'].destroy()
        self.remove_btn.configure(
            state="normal" if len(self.search_fields) > 1 else "disabled")
        self.filter_variants()

    def clear_all_search(self):
        for v in self.search_vars:
            v.set("")

    def clear_mutation_filter(self):
        self.min_mutations_var.set("")
        self.max_mutations_var.set("")

    def on_search_change(self, *_):
        if hasattr(self, '_filter_timer'):
            self.after_cancel(self._filter_timer)
        self._filter_timer = self.after(1000, self.filter_variants)

    def get_active_search_terms(self) -> List[str]:
        return [v.get().lower().strip()
                for v in self.search_vars if v.get().strip()]

    def get_mutation_count_filter(self):
        min_str = self.min_mutations_var.get().strip()
        max_str = self.max_mutations_var.get().strip()
        min_val = int(min_str) if min_str.isdigit() else None
        max_val = int(max_str) if max_str.isdigit() else None
        return min_val, max_val

    def count_mutations(self, mutation_string: str) -> int:
        if not mutation_string or mutation_string == "(none)":
            return 0
        return len([m for m in mutation_string.split(',') if m.strip()])

    # ------------------------------------------------------------------
    # Data loading & filtering
    # ------------------------------------------------------------------

    def load_variants(self):
        db = self.controller.get_databank()
        self.all_variants = {
            vid: muts
            for vid, muts in sorted(db.items(), key=lambda x: natural_sort_key(x[0]))
            if muts and muts != "(none)"
        }
        if not self.all_variants:
            for w in self.scroll_frame.winfo_children():
                w.destroy()
            self.variant_rows.clear()
            ctk.CTkLabel(self.scroll_frame,
                         text="No variants in databank. Run workflow first.",
                         font=self.controller.MEDIUMFONT, text_color="orange"
                         ).grid(row=0, column=0, pady=20)
            self.status_label.configure(
                text="No variants found", text_color="orange")
            return
        self.filter_variants()

    def filter_variants(self):
        terms    = self.get_active_search_terms()
        min_muts, max_muts = self.get_mutation_count_filter()

        if not terms and min_muts is None and max_muts is None:
            self.filtered_variants = dict(self.all_variants)
            self.update_variant_display()
            return

        if len(self.all_variants) > 100:
            self.status_label.configure(text="Filtering...", text_color="blue")
            self.update_idletasks()

        logic = self.search_logic_var.get()
        self.filtered_variants = {}
        for vid, muts in self.all_variants.items():
            text = f"{vid} {muts}".lower()
            if terms:
                match = (all(t in text for t in terms) if logic == "AND"
                         else any(t in text for t in terms))
                if not match:
                    continue
            if min_muts is not None or max_muts is not None:
                count = muts.count(',') + 1 if muts else 0
                if min_muts is not None and count < min_muts:
                    continue
                if max_muts is not None and count > max_muts:
                    continue
            self.filtered_variants[vid] = muts

        self.update_variant_display()

    # ------------------------------------------------------------------
    # Display
    # ------------------------------------------------------------------

    def update_variant_display(self):
        if hasattr(self, '_display_timer'):
            self.after_cancel(self._display_timer)
        if len(self.filtered_variants) > 200:
            self._display_timer = self.after(
                100, self._update_variant_display_now)
            return
        self._update_variant_display_now()

    def _update_variant_display_now(self):
        for w in self.scroll_frame.winfo_children():
            w.destroy()
        self.variant_rows.clear()

        total    = len(self.all_variants)
        filtered = len(self.filtered_variants)
        terms    = self.get_active_search_terms()
        min_muts, max_muts = self.get_mutation_count_filter()

        filter_desc: List[str] = []
        if terms:
            logic = self.search_logic_var.get()
            filter_desc.append(
                f"text: {f' {logic} '.join(repr(t) for t in terms)}")
        if min_muts is not None and max_muts is not None:
            filter_desc.append(f"mutations: {min_muts}-{max_muts}")
        elif min_muts is not None:
            filter_desc.append(f"mutations: ≥{min_muts}")
        elif max_muts is not None:
            filter_desc.append(f"mutations: ≤{max_muts}")

        if not filter_desc:
            self.results_label.configure(
                text=f"Showing all {total} variants with mutations")
        else:
            self.results_label.configure(
                text=f"Showing {filtered} of {total} "
                     f"matching [{'; '.join(filter_desc)}]")

        if not self.filtered_variants:
            ctk.CTkLabel(self.scroll_frame,
                         text="No variants match your search criteria.",
                         font=self.controller.MEDIUMFONT, text_color="orange"
                         ).grid(row=0, column=0, pady=20)
            self._update_scroll_region()
            return

        # Column widths: ID col fixed, one col per mutation
        self.scroll_frame.grid_columnconfigure(0, weight=0, minsize=150)
        max_mutations = max(
            len(m.split(',')) for m in self.filtered_variants.values())
        for col in range(1, max_mutations + 1):
            self.scroll_frame.grid_columnconfigure(col, weight=0, minsize=120)

        # Header row
        hdr_font = ctk.CTkFont(family="Verdana", size=14, weight="bold")
        ctk.CTkLabel(self.scroll_frame, text="Variant ID",
                     font=hdr_font
                     ).grid(row=0, column=0, padx=10, pady=10, sticky="w")
        ctk.CTkLabel(self.scroll_frame,
                     text="Mutations (Check if FAILED/MISSING)",
                     font=hdr_font
                     ).grid(row=0, column=1, padx=10, pady=10, sticky="w",
                            columnspan=max_mutations)

        medium_font = self.controller.MEDIUMFONT
        small_font  = self.controller.SMALLFONT

        for row_num, vid in enumerate(
            sorted(self.filtered_variants, key=natural_sort_key), start=1
        ):
            mut_str   = self.filtered_variants[vid]
            mutations = [m.strip() for m in mut_str.split(',')]

            id_label = ctk.CTkLabel(self.scroll_frame, text=vid,
                                    font=medium_font)
            id_label.grid(row=row_num, column=0, padx=10, pady=5, sticky="w")

            checkboxes: List[tuple] = []
            for col, mutation in enumerate(mutations, start=1):
                var = ctk.BooleanVar(value=False)
                cb  = ctk.CTkCheckBox(self.scroll_frame, text=mutation,
                                      variable=var, font=small_font)
                cb.grid(row=row_num, column=col, padx=5, pady=5, sticky="w")
                checkboxes.append((mutation, var))

            self.variant_rows[vid] = {
                'mutations':  mutations,
                'checkboxes': checkboxes,
                'label':      id_label,
                'row':        row_num,
            }

            if row_num % 20 == 0 and len(self.filtered_variants) > 50:
                self.scroll_frame.update_idletasks()
                self._update_scroll_region()

        self._update_scroll_region()
        self.status_label.configure(
            text=f"Displaying {len(self.variant_rows)} variants",
            text_color="green")

    def _update_scroll_region(self):
        if hasattr(self, '_canvas'):
            self.scroll_frame.update_idletasks()
            self._canvas.configure(
                scrollregion=self._canvas.bbox("all"))

    # ------------------------------------------------------------------
    # Actions
    # ------------------------------------------------------------------

    def clear_all_checks(self):
        for data in self.variant_rows.values():
            for _, var in data['checkboxes']:
                var.set(False)
        self.status_label.configure(
            text="All checks cleared", text_color="blue")

    def get_failed_mutations(self) -> Dict:
        result: Dict = {}
        for vid, data in self.variant_rows.items():
            failed  = [m for m, v in data['checkboxes'] if v.get()]
            present = [m for m, v in data['checkboxes'] if not v.get()]
            if failed:
                result[vid] = {
                    'failed':   failed,
                    'present':  present,
                    'original': data['mutations'],
                }
        return result

    def validate_selection(self):
        failed = self.get_failed_mutations()
        if not failed:
            messagebox.showinfo("No Failures",
                                "No failed mutations marked.")
            return

        db        = self.controller.get_databank()
        hierarchy = build_variant_hierarchy(db)
        marked    = {vid: set(info['failed']) for vid, info in failed.items()}

        def _task():
            variant_sets = {
                vid: frozenset(m.strip() for m in muts.split(',') if m.strip())
                for vid, muts in db.items()
                if muts and muts != "(none)"
            }
            warnings: List[str] = []
            for mv, failed_muts in marked.items():
                current  = mv
                unchecked: list = []
                while current in hierarchy:
                    parent      = hierarchy[current]
                    parent_muts = variant_sets.get(parent, frozenset())
                    overlap     = failed_muts & parent_muts
                    if overlap:
                        if parent in marked:
                            missing = overlap - marked[parent]
                            if missing:
                                unchecked.append((parent, missing))
                        else:
                            unchecked.append((parent, overlap))
                    current = parent
                if unchecked:
                    txt = f"\n⚠️ Variant {mv}:\n"
                    txt += f"   Failed: {', '.join(sorted(failed_muts))}\n"
                    txt += "   Potentially missed parents:\n"
                    for pid, pmuts in unchecked:
                        txt += (f"      • {pid}: check "
                                f"{', '.join(sorted(pmuts))}\n")
                    warnings.append(txt)
            self.after(0, lambda: _show(warnings))

        def _show(warnings):
            if warnings:
                win = ctk.CTkToplevel(self)
                win.title("Validation Warnings")
                win.geometry("700x500")
                tb = ctk.CTkTextbox(win, wrap="word")
                tb.pack(fill="both", expand=True, padx=20, pady=20)
                tb.insert("1.0",
                          "⚠️ VALIDATION WARNINGS ⚠️\n" + "=" * 60
                          + "\n" + "".join(warnings))
                tb.configure(state="disabled")
                ctk.CTkButton(win, text="Close",
                              command=win.destroy).pack(pady=10)
                win.transient(self)
                win.focus_set()
                self.status_label.configure(
                    text=f"⚠ {len(warnings)} potential issue(s)",
                    text_color="orange")
            else:
                messagebox.showinfo("Validation Passed",
                                    "✓ No inconsistencies detected!")
                self.status_label.configure(
                    text="✓ Validation passed", text_color="green")

        threading.Thread(target=_task, daemon=True).start()

    def process_failures(self):
        failed = self.get_failed_mutations()
        if not failed:
            messagebox.showinfo("No Failures",
                                "No failed mutations marked.")
            return

        summary = "\n".join(
            f"{vid}: {len(info['failed'])} failed mutation(s)"
            for vid, info in failed.items()
        )
        confirm = (
            f"Found failures in the following variants:\n\n{summary}\n\n"
            "This will:\n"
            "1. Update marked variants (remove failed mutations)\n"
            "2. Create new repair variants (with all mutations)\n"
            "3. Generate repair protocol\n\n"
            "⚠️ Only manually marked variants will be updated.\n"
            "Run 'Validate Selection' first to check for missed variants.\n\n"
            "Continue?"
        )
        if not messagebox.askyesno("Confirm Processing", confirm):
            return

        try:
            self.status_label.configure(
                text="Processing failures...", text_color="blue")
            db          = self.controller.get_databank()
            undo_backup = copy.deepcopy(db)

            ip_frame       = self.controller.frames[InputPage]
            variant_prefix = ip_frame.var_prefix.get().strip()
            existing_ids   = [
                int(v[len(variant_prefix):])
                for v in db
                if v.startswith(variant_prefix)
                and v[len(variant_prefix):].isdigit()
            ]
            next_id = max(existing_ids, default=0) + 1

            repair_variants:  list = []
            updated_variants: list = []
            updated_ids:      list = []
            repair_ids:       list = []

            for vid, info in failed.items():
                present_str = (','.join(info['present'])
                               if info['present'] else '(none)')
                db[vid] = present_str
                updated_variants.append((vid, present_str))
                updated_ids.append(vid)

                for fm in info['failed']:
                    new_id   = f"{variant_prefix}{next_id:02d}"
                    next_id += 1
                    all_muts = ','.join(info['original'])
                    db[new_id] = all_muts
                    repair_variants.append({
                        'new_id':           new_id,
                        'parent_id':        vid,
                        'missing_mutation': fm,
                        'all_mutations':    all_muts,
                    })
                    repair_ids.append(new_id)

            self.failure_undo_stack.append({
                'databank':    undo_backup,
                'updated_ids': updated_ids,
                'repair_ids':  repair_ids,
            })
            if len(self.failure_undo_stack) > 10:
                self.failure_undo_stack.pop(0)

            self.controller.save_databank(db)
            self.controller.mark_variants_modified(updated_ids, repair_ids)
            self._generate_repair_protocol(
                repair_variants, updated_variants, variant_prefix)

            self.load_variants()
            if DatabankPage in self.controller.frames:
                self.controller.frames[DatabankPage].refresh_variants()

            messagebox.showinfo(
                "Success",
                f"Updated {len(updated_variants)} variant(s)\n"
                f"Created {len(repair_variants)} repair variant(s)\n"
                "Repair protocol: protocols/repair_protocol.pdf")
            self.status_label.configure(
                text=f"Processed {len(failed)} variants",
                text_color="green")

        except Exception as e:
            messagebox.showerror("Error", f"Failed: {e}")
            self.status_label.configure(
                text=f"Error: {e}", text_color="red")

    def _generate_repair_protocol(self, repair_variants: list,
                                  updated_variants: list,
                                  variant_prefix: str):
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.styles import getSampleStyleSheet
        from reportlab.lib import colors
        from reportlab.platypus import (SimpleDocTemplate, Paragraph,
                                        Spacer, Table, TableStyle)

        out_dir  = Path(self.controller.output_dir.get())
        pdf_path = out_dir / "protocols" / "repair_protocol.pdf"
        pdf_path.parent.mkdir(parents=True, exist_ok=True)

        doc      = SimpleDocTemplate(str(pdf_path), pagesize=A4)
        styles   = getSampleStyleSheet()
        elements = [
            Paragraph("<b>Repair Mutagenesis Protocol</b>", styles['Title']),
            Spacer(1, 20),
            Paragraph("<b>1. Updated Variants</b>", styles['Heading2']),
            Spacer(1, 12),
        ]

        upd_data = [['Variant ID', 'Remaining Mutations']]
        for vid, muts in updated_variants:
            upd_data.append([vid, muts])
        upd_tbl = Table(upd_data, colWidths=[150, 350])
        upd_tbl.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.lightgrey),
            ('GRID',       (0, 0), (-1, -1), 0.5, colors.grey),
        ]))
        elements += [upd_tbl, Spacer(1, 20),
                     Paragraph("<b>2. Repair Mutagenesis Steps</b>",
                               styles['Heading2']),
                     Spacer(1, 12)]

        rep_data = [['Step', 'New Variant', 'Parent',
                     'Add Missing', 'Final Mutations']]
        for i, r in enumerate(repair_variants, 1):
            rep_data.append([str(i), r['new_id'], r['parent_id'],
                             r['missing_mutation'], r['all_mutations']])
        rep_tbl = Table(rep_data, colWidths=[40, 100, 100, 120, 150])
        rep_tbl.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.lightgrey),
            ('GRID',       (0, 0), (-1, -1), 0.5, colors.grey),
            ('FONTSIZE',   (0, 0), (-1, -1), 9),
        ]))
        elements.append(rep_tbl)
        doc.build(elements)

    def undo_last_processing(self):
        if not self.failure_undo_stack:
            messagebox.showinfo("Nothing to Undo",
                                "No previous failure processing to undo.")
            self.status_label.configure(
                text="No undo available", text_color="orange")
            return
        try:
            undo_data = self.failure_undo_stack.pop()
            self.controller.save_databank(undo_data['databank'])
            repair_pdf = (Path(self.controller.output_dir.get())
                          / "protocols" / "repair_protocol.pdf")
            if repair_pdf.exists():
                repair_pdf.unlink()
            self.load_variants()
            if DatabankPage in self.controller.frames:
                self.controller.frames[DatabankPage].refresh_variants()
            messagebox.showinfo(
                "Undo Successful",
                f"Reverted:\n"
                f"• {len(undo_data['updated_ids'])} updated variant(s)\n"
                f"• {len(undo_data['repair_ids'])} repair variant(s)")
            self.status_label.configure(
                text="Undo successful - previous state restored",
                text_color="green")
        except Exception as e:
            messagebox.showerror("Undo Error", f"Failed to undo:\n{e}")


# ===========================================================================
# LabelSelectionWindow
# ===========================================================================

class LabelSelectionWindow:
    def __init__(self, parent, selected_variants: list,
                 output_dir: str, controller):
        self.controller       = controller
        self.selected_variants = selected_variants
        self.output_dir       = Path(output_dir)
        self.label_system     = LabelPrintingSystem()

        self.window = ctk.CTkToplevel(parent)
        self.window.title("Print Labels")
        self.window.geometry("400x300")
        self.window.transient(parent)
        self._setup_ui()
        self.window.after(100, self.window.grab_set)

    def _setup_ui(self):
        ctk.CTkLabel(self.window,
                     text=f"Print labels for {len(self.selected_variants)} "
                          "variant(s)",
                     font=self.controller.MEDIUMFONT
                     ).pack(pady=20)

        pf = ctk.CTkFrame(self.window)
        pf.pack(fill="x", padx=20, pady=10)
        pf.grid_columnconfigure(1, weight=1)

        ctk.CTkLabel(pf, text="Start column (0-6):").grid(
            row=0, column=0, padx=10, pady=5, sticky="w")
        self.start_col_var = ctk.IntVar(value=0)
        ctk.CTkEntry(pf, textvariable=self.start_col_var, width=60
                     ).grid(row=0, column=1, padx=10, pady=5, sticky="w")

        ctk.CTkLabel(pf, text="Start row (0-26):").grid(
            row=1, column=0, padx=10, pady=5, sticky="w")
        self.start_row_var = ctk.IntVar(value=0)
        ctk.CTkEntry(pf, textvariable=self.start_row_var, width=60
                     ).grid(row=1, column=1, padx=10, pady=5, sticky="w")

        bf = ctk.CTkFrame(self.window)
        bf.pack(fill="x", padx=20, pady=10)
        ctk.CTkButton(bf, text="Generate PDF",
                      command=self._generate, fg_color="green"
                      ).grid(row=0, column=0, padx=10, pady=10)
        ctk.CTkButton(bf, text="Cancel",
                      command=self.window.destroy
                      ).grid(row=0, column=1, padx=10, pady=10)

    def _generate(self):
        out = (self.output_dir / "protocols"
               / "selected_variant_labels.pdf")
        out.parent.mkdir(parents=True, exist_ok=True)
        try:
            self.label_system.generate_labels_pdf(
                self.selected_variants,
                self.start_row_var.get(),
                self.start_col_var.get(),
                out,
            )
            messagebox.showinfo("Success", f"Labels saved:\n{out}")
            self.window.destroy()
        except Exception as e:
            messagebox.showerror("Error", f"Failed:\n{e}")


# ===========================================================================
# SLiCEDesignerPage
# ===========================================================================

class SLiCEDesignerPage(ctk.CTkFrame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        self.sequence   = ""
        self.fragments: List[SLiCEFragment] = []
        self.primers:   List[SLiCEPrimer]   = []
        self.designer:  Optional[SLiCEPrimerDesigner] = None

        self.grid_columnconfigure(0, weight=2)
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(1, weight=1)

        self._setup_ui()
        self._update_designer()

    def _setup_ui(self):
        hf = ctk.CTkFrame(self)
        hf.grid(row=0, column=0, columnspan=2, sticky="ew",
                padx=10, pady=20)
        hf.grid_columnconfigure(0, weight=1)
        ctk.CTkLabel(hf, text="SLiCE Assembly Designer",
                     corner_radius=10, fg_color="#4a90e2",
                     text_color="white", font=self.controller.LARGEFONT,
                     height=50
                     ).grid(row=0, column=0, padx=10, pady=10, sticky="ew")
        ctk.CTkButton(hf, text="Back to Input",
                      command=lambda: self.controller.show_frame(InputPage)
                      ).grid(row=0, column=1, padx=10, pady=10, sticky="e")

        self._create_sequence_viewer()
        self._create_control_panel()

    def _create_sequence_viewer(self):
        vf = ctk.CTkFrame(self)
        vf.grid(row=1, column=0, sticky="nsew", padx=(10, 5), pady=10)
        vf.grid_columnconfigure(0, weight=1)
        vf.grid_rowconfigure(2, weight=1)

        tf = ctk.CTkFrame(vf)
        tf.grid(row=0, column=0, sticky="ew", padx=10, pady=10)
        tf.grid_columnconfigure(0, weight=1)
        ctk.CTkLabel(tf, text="Plasmid Sequence",
                     font=self.controller.MEDIUMFONT
                     ).grid(row=0, column=0, sticky="w")
        ctk.CTkButton(tf, text="Load Plasmid File",
                      command=self._load_plasmid
                      ).grid(row=0, column=1, padx=5)
        ctk.CTkButton(tf, text="Use WT Sequence",
                      command=self._use_wildtype
                      ).grid(row=0, column=2, padx=5)

        info = ctk.CTkFrame(vf)
        info.grid(row=1, column=0, sticky="ew", padx=10, pady=5)
        self.seq_info_label = ctk.CTkLabel(info, text="No sequence loaded",
                                           font=self.controller.SMALLFONT)
        self.seq_info_label.pack(side="left", padx=10)

        self.seq_text = ctk.CTkTextbox(
            vf, wrap="char",
            font=ctk.CTkFont(family="Courier", size=10))
        self.seq_text.grid(row=2, column=0, sticky="nsew", padx=10, pady=5)

        sel = ctk.CTkFrame(vf)
        sel.grid(row=3, column=0, sticky="ew", padx=10, pady=10)
        ctk.CTkLabel(sel, text="Select Fragment:",
                     font=self.controller.SMALLFONT
                     ).grid(row=0, column=0, padx=5)
        ctk.CTkLabel(sel, text="Start:").grid(row=0, column=1, padx=5)
        self.start_entry = ctk.CTkEntry(sel, width=80, placeholder_text="1")
        self.start_entry.grid(row=0, column=2, padx=5)
        ctk.CTkLabel(sel, text="End:").grid(row=0, column=3, padx=5)
        self.end_entry = ctk.CTkEntry(sel, width=80, placeholder_text="100")
        self.end_entry.grid(row=0, column=4, padx=5)
        ctk.CTkButton(sel, text="Add Fragment", fg_color="green",
                      command=self._add_fragment_manual
                      ).grid(row=0, column=5, padx=10)

    def _create_control_panel(self):
        cf = ctk.CTkFrame(self)
        cf.grid(row=1, column=1, sticky="nsew", padx=(5, 10), pady=10)
        cf.grid_columnconfigure(0, weight=1)
        cf.grid_rowconfigure(2, weight=1)
        cf.grid_rowconfigure(4, weight=2)

        # Parameters
        pf = ctk.CTkFrame(cf)
        pf.grid(row=0, column=0, sticky="ew", padx=10, pady=10)
        pf.grid_columnconfigure(1, weight=1)

        self.overlap_min_var = ctk.IntVar(value=15)
        self.overlap_max_var = ctk.IntVar(value=25)
        self.overlap_tm_min  = ctk.IntVar(value=55)
        self.overlap_tm_max  = ctk.IntVar(value=65)
        self.circular_var    = ctk.BooleanVar(value=True)

        for row, (lbl, var) in enumerate([
            ("Min overlap (bp):",   self.overlap_min_var),
            ("Max overlap (bp):",   self.overlap_max_var),
            ("Min overlap Tm (°C):", self.overlap_tm_min),
            ("Max overlap Tm (°C):", self.overlap_tm_max),
        ]):
            ctk.CTkLabel(pf, text=lbl,
                         font=self.controller.SMALLFONT
                         ).grid(row=row, column=0, sticky="w", padx=5, pady=3)
            ctk.CTkEntry(pf, textvariable=var, width=60
                         ).grid(row=row, column=1, sticky="w",
                                padx=5, pady=3)

        ctk.CTkCheckBox(pf, text="Circular assembly",
                        variable=self.circular_var
                        ).grid(row=4, column=0, columnspan=2,
                               padx=5, pady=5, sticky="w")

        # Fragment list
        ctk.CTkLabel(cf, text="Fragments",
                     font=self.controller.MEDIUMFONT
                     ).grid(row=1, column=0, sticky="nw", padx=10, pady=(10, 5))
        self.fragment_list = ctk.CTkTextbox(cf, height=100)
        self.fragment_list.grid(row=2, column=0, sticky="nsew",
                                padx=10, pady=(0, 5))
        self.fragment_list.configure(state="disabled")

        # Actions
        af = ctk.CTkFrame(cf)
        af.grid(row=3, column=0, sticky="ew", padx=10, pady=10)
        ctk.CTkButton(af, text="Design Primers",
                      command=self._design_primers,
                      font=self.controller.MEDIUMFONT,
                      fg_color="blue", height=40
                      ).pack(fill="x", pady=5)
        ctk.CTkButton(af, text="Validate Assembly",
                      command=self._validate_assembly
                      ).pack(fill="x", pady=5)
        ctk.CTkButton(af, text="Export Primers CSV",
                      command=self._export_primers
                      ).pack(fill="x", pady=5)
        ctk.CTkButton(af, text="Clear All Fragments",
                      command=self._clear_fragments
                      ).pack(fill="x", pady=5)

        # Results
        ctk.CTkLabel(cf, text="Primer Results",
                     font=self.controller.MEDIUMFONT
                     ).grid(row=4, column=0, sticky="nw", padx=10, pady=(10, 5))
        self.results_text = ctk.CTkTextbox(
            cf, font=ctk.CTkFont(family="Courier", size=9))
        self.results_text.grid(row=4, column=0, sticky="nsew",
                               padx=10, pady=(0, 10))

        self.status_label = ctk.CTkLabel(cf, text="",
                                         font=self.controller.SMALLFONT)
        self.status_label.grid(row=5, column=0, sticky="w", padx=10, pady=5)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _update_designer(self):
        self.designer = SLiCEPrimerDesigner(
            overlap_length_range=(self.overlap_min_var.get(),
                                  self.overlap_max_var.get()),
            overlap_tm_range=(self.overlap_tm_min.get(),
                              self.overlap_tm_max.get()),
        )

    def _load_plasmid(self):
        path = filedialog.askopenfilename(
            title="Load Plasmid Sequence",
            filetypes=[("FASTA/GenBank", "*.fasta *.fa *.gb *.gbk"),
                       ("All files", "*.*")])
        if not path:
            return
        try:
            with open(path, 'r') as f:
                content = f.read()
            if content.startswith('>'):
                seq = ''.join(line.strip()
                              for line in content.splitlines()
                              if not line.startswith('>'))
            else:
                seq = ''.join(content.split())
            seq = re.sub(r'[^ATGCatgc]', '', seq).upper()
            self._load_sequence(seq)
        except Exception as e:
            messagebox.showerror("Error", f"Could not load plasmid:\n{e}")

    def _use_wildtype(self):
        try:
            dna = self.controller.frames[InputPage].get_dna_sequence()
            if dna:
                self._load_sequence(dna)
            else:
                messagebox.showwarning("No Sequence",
                                       "Please enter a sequence in the Input page.")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def _load_sequence(self, seq: str):
        self.sequence = seq.upper()
        annotations   = PlasmidAnnotator.annotate_sequence(self.sequence)
        self.seq_info_label.configure(
            text=f"Sequence: {len(self.sequence)} bp | "
                 f"Annotations: {len(annotations)}")
        self.seq_text.delete("1.0", "end")
        for i in range(0, len(self.sequence), 60):
            chunk = self.sequence[i:i+60]
            fmt   = ' '.join(chunk[j:j+10] for j in range(0, len(chunk), 10))
            self.seq_text.insert("end", f"{i+1:>6}: {fmt}\n")

    def _add_fragment_manual(self):
        if not self.sequence:
            messagebox.showwarning("No Sequence",
                                   "Please load a sequence first.")
            return
        try:
            start = int(self.start_entry.get()) - 1
            end   = int(self.end_entry.get())
            if start < 0 or end > len(self.sequence) or start >= end:
                raise ValueError("Invalid range")
            name = ctk.CTkInputDialog(
                text=f"Enter fragment name ({start+1}-{end}, "
                     f"{end-start} bp):",
                title="Fragment Name",
            ).get_input()
            if not name:
                return
            self.fragments.append(SLiCEFragment(
                name=name, start=start, end=end,
                sequence=self.sequence[start:end],
            ))
            self._update_fragment_list()
            self.status_label.configure(
                text=f"Added fragment: {name}", text_color="green")
        except ValueError:
            messagebox.showerror("Invalid Input",
                                 "Please enter valid start and end positions.")

    def _update_fragment_list(self):
        self.fragment_list.configure(state="normal")
        self.fragment_list.delete("1.0", "end")
        if not self.fragments:
            self.fragment_list.insert("1.0", "No fragments defined")
        else:
            for i, f in enumerate(self.fragments, 1):
                self.fragment_list.insert(
                    "end",
                    f"{i}. {f.name}\n"
                    f"   {f.start+1}-{f.end} ({len(f.sequence)} bp)\n\n")
        self.fragment_list.configure(state="disabled")

    def _clear_fragments(self):
        if self.fragments and messagebox.askyesno(
            "Confirm", "Clear all fragments?"
        ):
            self.fragments.clear()
            self._update_fragment_list()
            self.status_label.configure(
                text="Fragments cleared", text_color="blue")

    def _design_primers(self):
        if not self.fragments:
            messagebox.showwarning("No Fragments",
                                   "Please define fragments first.")
            return
        self._update_designer()
        try:
            self.primers = self.designer.design_primers_for_fragments(
                self.fragments, circular=self.circular_var.get())
            self._display_results()
            self.status_label.configure(
                text=f"Designed {len(self.primers)} primers",
                text_color="green")
        except Exception as e:
            messagebox.showerror("Error", f"Primer design failed:\n{e}")

    def _display_results(self):
        self.results_text.delete("1.0", "end")
        if not self.primers:
            self.results_text.insert("1.0", "No primers designed")
            return
        by_frag: Dict = {}
        for p in self.primers:
            by_frag.setdefault(p.fragment_name, []).append(p)
        for frag_name, prims in by_frag.items():
            self.results_text.insert("end", f"=== {frag_name} ===\n")
            for p in prims:
                dir_txt = "Fwd" if p.is_forward else "Rev"
                self.results_text.insert(
                    "end",
                    f"{p.name} ({dir_txt}):\n"
                    f"  {p.sequence}\n"
                    f"  {p.length}bp, Tm:{p.tm:.1f}°C"
                )
                if p.overlap_length > 0:
                    self.results_text.insert(
                        "end", f", Ovl:{p.overlap_length}bp")
                self.results_text.insert("end", "\n\n")

    def _validate_assembly(self):
        if not self.fragments:
            messagebox.showwarning("No Fragments",
                                   "Please define fragments first.")
            return
        self._update_designer()
        v = self.designer.validate_assembly(
            self.fragments, circular=self.circular_var.get())
        if v['valid']:
            messagebox.showinfo("Valid",
                                f"✓ Assembly valid!\n\n"
                                f"Fragments: {v['total_fragments']}\n"
                                f"Total: {v['total_length']} bp")
        else:
            issues = '\n'.join(f"• {i}" for i in v['issues'])
            messagebox.showwarning("Issues",
                                   f"⚠ Problems found:\n\n{issues}")

    def _export_primers(self):
        if not self.primers:
            messagebox.showwarning("No Primers",
                                   "Design primers first.")
            return
        path = filedialog.asksaveasfilename(
            title="Save Primers CSV",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        if path:
            try:
                io.save_slice_primer_csv(self.primers, Path(path))
                messagebox.showinfo("Success", f"Exported to:\n{path}")
            except Exception as e:
                messagebox.showerror("Error", f"Export failed:\n{e}")


# ===========================================================================
# Entry point
# ===========================================================================

if __name__ == "__main__":
    app = MutagenesisApp()
    app.mainloop()
