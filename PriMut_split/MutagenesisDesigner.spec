# -*- mode: python ; coding: utf-8 -*-
#
# MutagenesisDesigner.spec
# ------------------------
# Updated from the original spec to work with the refactored three-file layout:
#   primer_core.py  — pure logic
#   primer_io.py    — file I/O
#   primer_gui.py   — GUI entry point  ← PyInstaller starts here
#
# Usage (from the folder containing all .py files):
#   pyinstaller MutagenesisDesigner.spec
#
# Output: dist/MutagenesisDesigner.exe

import sys
import os
from PyInstaller.utils.hooks import (
    collect_data_files,
    collect_dynamic_libs,
    collect_submodules,
)

block_cipher = None

datas        = []
binaries     = []
hiddenimports = []

# ===========================================================================
# PRIMER3
# ===========================================================================
print("Collecting primer3 files...")
try:
    primer3_binaries = collect_dynamic_libs('primer3')
    binaries.extend(primer3_binaries)
    print(f"  Found {len(primer3_binaries)} primer3 binaries")

    primer3_datas = collect_data_files('primer3')
    datas.extend(primer3_datas)
    print(f"  Found {len(primer3_datas)} primer3 data files")

    primer3_imports = collect_submodules('primer3')
    hiddenimports.extend(primer3_imports)
    print(f"  Found {len(primer3_imports)} primer3 submodules")
except Exception as e:
    print(f"  WARNING: Could not collect primer3 files: {e}")
    print("  Tm values will show as 0.0 — reinstall primer3-py to fix")

# ===========================================================================
# CUSTOMTKINTER
# ===========================================================================
print("Collecting customtkinter files...")
try:
    import customtkinter
    customtkinter_path = os.path.dirname(customtkinter.__file__)

    # Assets folder (themes, icons, images)
    assets_path = os.path.join(customtkinter_path, 'assets')
    if os.path.exists(assets_path):
        datas.append((assets_path, 'customtkinter/assets'))
        print(f"  Added customtkinter assets from: {assets_path}")

    ctk_imports = collect_submodules('customtkinter')
    hiddenimports.extend(ctk_imports)
    print(f"  Found {len(ctk_imports)} customtkinter submodules")
except Exception as e:
    print(f"  WARNING: Could not collect customtkinter files: {e}")

# ===========================================================================
# REPORTLAB
# ===========================================================================
print("Collecting reportlab files...")
try:
    reportlab_datas = collect_data_files('reportlab')
    datas.extend(reportlab_datas)

    reportlab_imports = collect_submodules('reportlab')
    hiddenimports.extend(reportlab_imports)
    print(f"  Found {len(reportlab_datas)} reportlab data files")
    print(f"  Found {len(reportlab_imports)} reportlab submodules")
except Exception as e:
    print(f"  WARNING: Could not collect reportlab files: {e}")

# ===========================================================================
# MANUAL HIDDEN IMPORTS
# ===========================================================================
hiddenimports.extend([
    # --- our own modules (PyInstaller finds them via import, but be explicit) ---
    'primer_core',
    'primer_io',
    # --- primer3 ---
    'primer3',
    # --- GUI ---
    'customtkinter',
    'tkinter',
    'tkinter.ttk',
    'tkinter.filedialog',
    'tkinter.messagebox',
    'tkinter.simpledialog',
    # --- PIL (needed by customtkinter for image handling) ---
    'PIL',
    'PIL._tkinter_finder',
    'PIL.Image',
    'PIL.ImageTk',
    # --- reportlab (runtime-loaded sub-modules) ---
    'reportlab',
    'reportlab.pdfgen',
    'reportlab.pdfgen.canvas',
    'reportlab.lib',
    'reportlab.lib.pagesizes',
    'reportlab.lib.styles',
    'reportlab.lib.colors',
    'reportlab.lib.units',
    'reportlab.platypus',
    'reportlab.pdfbase.pdfdoc',
    'reportlab.pdfbase._fontdata',
    'reportlab.pdfbase.pdfmetrics',
    'reportlab.graphics.barcode.code39',
    'reportlab.graphics.barcode.code128',
    # --- stdlib (always safe to be explicit) ---
    'json',
    'csv',
    'pathlib',
    're',
    'collections',
    'webbrowser',
    'threading',
    'copy',
    'shutil',
])

# Remove duplicates
hiddenimports = list(set(hiddenimports))

print(f"\nTotal collected:")
print(f"  - {len(binaries)} binaries")
print(f"  - {len(datas)} data files")
print(f"  - {len(hiddenimports)} hidden imports")
print()

# ===========================================================================
# ANALYSIS
# Entry point is now primer_gui.py — PyInstaller follows the imports to
# primer_core and primer_io automatically. No need to list them separately.
# ===========================================================================
a = Analysis(
    ['primer_gui.py'],          # ← changed from primer_windows.py
    pathex=['.'],               # look in current directory for the other modules
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='MutagenesisDesigner',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,                  # False = faster startup (your original setting)
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,              # no CMD window
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon='primut.ico',          # your original icon
)
