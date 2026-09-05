# scpecg.spec - build per macOS
#
#   pip install pyinstaller
#   pyinstaller --clean --noconfirm scpecg.spec
#
# Produce dist/scpecg.app (doppio clic) e dist/scpecg/ (cartella).
# Da riga di comando l'eseguibile vero sta dentro il bundle:
#   ./dist/scpecg.app/Contents/MacOS/scpecg --info 13.SCP
#
# NOTA: PyInstaller non compila per un sistema diverso da quello su cui gira.
# Da macOS esce solo il .app; per Windows serve una macchina Windows.

import sys
from PyInstaller.utils.hooks import collect_data_files

# matplotlib porta con se' font e file di stile che non vengono raccolti da
# soli: senza questi il programma parte e crasha al primo disegno.
datas = collect_data_files("matplotlib")

hiddenimports = [
    "matplotlib.backends.backend_tkagg",
    "matplotlib.backends.backend_agg",
]

# Tutto quello che non serve: senza escluderlo il bundle raddoppia. Le GUI
# alternative in particolare vengono trascinate dentro da matplotlib, che
# supporta piu' backend di quanti ne usiamo.
excludes = [
    "PyQt5", "PyQt6", "PySide2", "PySide6", "wx",
    "IPython", "jupyter", "notebook", "pytest",
    "scipy", "pandas", "sympy", "PIL.ImageQt",
    "matplotlib.backends.backend_qt5agg",
    "matplotlib.backends.backend_qtagg",
    "matplotlib.backends.backend_wxagg",
    "matplotlib.backends.backend_webagg",
    "matplotlib.backends.backend_gtk3agg",
    "matplotlib.backends.backend_gtk4agg",
]

a = Analysis(
    ["scpecg.py"],
    pathex=[],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=excludes,
    noarchive=False,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="scpecg",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,          # UPX rompe la firma su macOS e insospettisce gli antivirus
    console=False,      # niente finestra di terminale al doppio clic
    disable_windowed_traceback=False,
    icon="icona.icns",
    argv_emulation=False,
    target_arch=None,   # nativo: arm64 su Apple Silicon
    codesign_identity=None,
    entitlements_file=None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    name="scpecg",
)

if sys.platform == "darwin":
    app = BUNDLE(
        coll,
        name="scpecg.app",
        icon="icona.icns",
        bundle_identifier="it.silfox.scpecg",
        info_plist={
            "CFBundleName": "scpecg",
            "CFBundleDisplayName": "scpecg",
            "CFBundleShortVersionString": "1.0.0",
            "CFBundleVersion": "1.0.0",
            "NSHighResolutionCapable": True,
            "LSMinimumSystemVersion": "11.0",
            "NSHumanReadableCopyright": "Copyright 2026 Silvestro Scuderi - GPL-3.0-or-later",
        },
    )
