# Installation
Download the executable directly from the Releases section on the right-hand side of the repository. Simply open the latest release and select the appropriate file for your system.

Or to create the .exe file yourself. Download the .spec file and run:
```
pyinstaller primer_windows.spec
```


# PriMut: Primer Design, Protocol Generation, and Protein Variant Management
This repository provides a complete workflow for site-directed mutagenesis experiment planning with a modern Windows GUI — from primer design and mutation step planning to generating printable protocols and tube labels. It is tailored for molecular biology workflows involving multi-site mutagenesis with automated primer design and detailed PDF output.

✨ Key Features
- Automated Primer Design

  - Validates wildtype sequence start/stop codons and flanking regions.
  - Groups nearby mutations into single primer sets (adjustable distance threshold).
  - Uses primer3‑py to optimize melting temperatures and primer lengths.
  - Outputs TXT, CSV, and JSON primer files.

- Smart Stepwise Mutagenesis Planning
  - Calculates minimal mutation pathways, reusing intermediate variants where possible.
  - Supports configurable maximum mutations per step.
  - Maintains a persistent variant_databank.json to avoid duplicate work.

- Protocol PDF Generation
  - Rich formatting with overview tables, per-step details, and all variants generated.
  - Includes notice of pre-existing variants already in databank.

- Label Sheet PDF Generation
  - Automatically prints labels for all variants from the current run (final + intermediates).
  - Customizable start column/row to reuse partially used sheets.
  - Uses standard Hema 10900 label layout.

- Windows GUI with CustomTkinter
  - Interactive desktop application for inputting sequences and mutations.
  - Real-time validation and display of wildtype and variant protein sequences.
  - Easy-to-use workflow controls for primer design and protocol generation.
  - Output management directly integrated into GUI.

📂 Outputs
- primer_list.txt / .csv / .json – primer sequences and statistics.
- mutagenesis_protocol.pdf – full mutagenesis plan, with steps and overview tables.
- herma_10900_labels.pdf – printable tube labels.
- variant_databank.json – persistent variant registry.

🛠 Requirements
- Python 3.x
  - primer3-py (primer design)
  - reportlab (PDF generation)
  - customtkinter (Windows GUI toolkit)

© 2026 Konstantin Weigmann — licensed under GPL-3.0
