#!/usr/bin/env python3
"""Genera l'icona di scpecg: icona.png, icona.icns (macOS), icona.ico (Windows).

Il disegno e' un complesso QRS bianco su fondo scuro, con un accenno di carta
millimetrata. Deve restare leggibile a 32 pixel, quindi niente dettagli fini:
tratto spesso e una sola forma riconoscibile.
"""

import struct
from io import BytesIO

from PIL import Image, ImageDraw

S = 1024
FONDO = (32, 42, 56)
GRIGLIA = (58, 72, 90)
GRIGLIA_FORTE = (78, 96, 118)
TRATTO = (240, 244, 248)


def disegna(size=S):
    img = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    r = int(size * 0.22)
    d.rounded_rectangle([0, 0, size - 1, size - 1], radius=r, fill=FONDO)

    # carta millimetrata appena accennata
    passo = size // 20
    for i in range(1, 20):
        x = i * passo
        col = GRIGLIA_FORTE if i % 5 == 0 else GRIGLIA
        w = max(1, size // 340) if i % 5 == 0 else max(1, size // 680)
        d.line([(x, 0), (x, size)], fill=col, width=w)
        d.line([(0, x), (size, x)], fill=col, width=w)

    # un battito: linea di base, P, QRS, T
    u = size / 100.0
    pts = [
        (6, 50), (18, 50), (22, 44), (26, 50), (36, 50),
        (40, 56), (46, 18), (52, 74), (56, 50),
        (66, 50), (72, 40), (80, 50), (94, 50),
    ]
    linea = [(x * u, y * u) for x, y in pts]
    d.line(linea, fill=TRATTO, width=int(size * 0.055),
           joint="curve")

    # gli angoli arrotondati devono restare trasparenti
    maschera = Image.new("L", (size, size), 0)
    ImageDraw.Draw(maschera).rounded_rectangle(
        [0, 0, size - 1, size - 1], radius=r, fill=255)
    img.putalpha(maschera)
    return img


def scrivi_icns(img, path):
    """Compone un .icns con voci PNG. Il formato e' un contenitore semplice:
    magic 'icns', lunghezza totale, poi blocchi tipo+lunghezza+dati."""
    tipi = {
        b"ic11": 32, b"ic12": 64, b"ic07": 128,
        b"ic08": 256, b"ic09": 512, b"ic10": 1024,
    }
    blocchi = b""
    for tipo, lato in tipi.items():
        buf = BytesIO()
        img.resize((lato, lato), Image.LANCZOS).save(buf, format="PNG")
        dati = buf.getvalue()
        blocchi += tipo + struct.pack(">I", len(dati) + 8) + dati
    with open(path, "wb") as fh:
        fh.write(b"icns" + struct.pack(">I", len(blocchi) + 8) + blocchi)


if __name__ == "__main__":
    base = disegna()
    base.save("icona.png")
    scrivi_icns(base, "icona.icns")
    base.save("icona.ico", sizes=[(16, 16), (32, 32), (48, 48),
                                  (64, 64), (128, 128), (256, 256)])
    print("scritti: icona.png, icona.icns, icona.ico")
