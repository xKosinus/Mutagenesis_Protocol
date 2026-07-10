# -*- mode: python ; coding: utf-8 -*-
import sys
import os
from PyInstaller.utils.hooks import collect_data_files, collect_dynamic_libs, collect_submodules

block_cipher = None

# Sammle alle notwendigen Dateien und Bibliotheken
datas = []
binaries = []
hiddenimports = []

# ===== PRIMER3 =====
print("Collecting primer3 files...")
try:
    # Sammle alle primer3 Bibliotheken (DLLs, .pyd Dateien)
    primer3_binaries = collect_dynamic_libs('primer3')
    binaries.extend(primer3_binaries)
    print(f"  Found {len(primer3_binaries)} primer3 binaries")
    
    # Sammle alle primer3 Daten
    primer3_datas = collect_data_files('primer3')
    datas.extend(primer3_datas)
    print(f"  Found {len(primer3_datas)} primer3 data files")
    
    # Sammle alle primer3 Submodule
    primer3_imports = collect_submodules('primer3')
    hiddenimports.extend(primer3_imports)
    print(f"  Found {len(primer3_imports)} primer3 submodules")
    
except Exception as e:
    print(f"  WARNING: Could not collect primer3 files: {e}")
    print("  Program will use fallback Tm calculation")

# ===== CUSTOMTKINTER =====
print("Collecting customtkinter files...")
try:
    import customtkinter
    customtkinter_path = os.path.dirname(customtkinter.__file__)
    
    # Assets Ordner
    assets_path = os.path.join(customtkinter_path, 'assets')
    if os.path.exists(assets_path):
        datas.append((assets_path, 'customtkinter/assets'))
        print(f"  Added customtkinter assets from: {assets_path}")
    
    # Submodule
    ctk_imports = collect_submodules('customtkinter')
    hiddenimports.extend(ctk_imports)
    print(f"  Found {len(ctk_imports)} customtkinter submodules")
    
except Exception as e:
    print(f"  WARNING: Could not collect customtkinter files: {e}")

# ===== REPORTLAB =====
print("Collecting reportlab files...")
try:
    reportlab_datas = collect_data_files('reportlab')
    datas.extend(reportlab_datas)
    reportlab_imports = collect_submodules('reportlab')
    hiddenimports.extend(reportlab_imports)
    print(f"  Found {len(reportlab_datas)} reportlab data files")
except Exception as e:
    print(f"  WARNING: Could not collect reportlab files: {e}")

# ===== MANUELLE HIDDEN IMPORTS =====
hiddenimports.extend([
    'primer3',
    'customtkinter',
    'tkinter',
    'tkinter.ttk',
    'tkinter.filedialog',
    'tkinter.messagebox',
    'tkinter.simpledialog',
    'PIL',
    'PIL._tkinter_finder',
    'PIL.Image',
    'PIL.ImageTk',
    'reportlab',
    'reportlab.pdfgen',
    'reportlab.pdfgen.canvas',
    'reportlab.lib',
    'reportlab.lib.pagesizes',
    'reportlab.lib.styles',
    'reportlab.lib.colors',
    'reportlab.lib.units',
    'reportlab.platypus',
    'json',
    'csv',
    'pathlib',
    're',
    'collections',
    'webbrowser',
])

# Entferne Duplikate
hiddenimports = list(set(hiddenimports))

print(f"\nTotal collected:")
print(f"  - {len(binaries)} binaries")
print(f"  - {len(datas)} data files")
print(f"  - {len(hiddenimports)} hidden imports")
print()

a = Analysis(
    ['primer_windows.py'],
    pathex=[],
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
    upx=False,  # GEÄNDERT: False für schnelleren Start
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,  # GEÄNDERT: False um CMD-Fenster zu verstecken
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon='primut.ico',
)
