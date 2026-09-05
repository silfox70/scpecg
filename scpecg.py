#!/usr/bin/env python3
"""
scpecg.py - lettore, convertitore, visualizzatore e analizzatore SCP-ECG
            (EN 1064 / ISO 11073-91064)

    ./scpecg.py                        finestra, chiede il file
    ./scpecg.py 11.SCP                 apre il file nella finestra
    ./scpecg.py 13.SCP 14.SCP 15.SCP   unisce i segmenti e li apre come uno solo
    ./scpecg.py --info 13.SCP          struttura e analisi, senza produrre file
    ./scpecg.py --batch *.SCP          converte ogni file separatamente

In misura continua il PC-80B spezza la registrazione in record da 30 secondi:
passandoli tutti insieme vengono uniti, dopo aver verificato sui timestamp che
la sequenza sia completa. Se manca un segmento il programma si ferma.

Nella finestra: cursore a mirino con lettura di tempo e millivolt, caliper a
due punti, battiti colorati per tipo di intervallo, guadagno e velocita' carta
regolabili, esportazione in CSV e PNG.

Scritto per un Creative/Gima PC-80B, che scrive SCP-ECG v1.3 con due liberta'
rispetto allo standard: codice derivazione proprietario (101) e bit 15 dei
campioni usato come marcatore dei picchi R. Entrambe riconosciute da sole.

L'analisi del ritmo usa il dizionario SOGLIE, in parte dal manuale e in parte
scelto qui: NON riproduce il verdetto dell'apparecchio e non e' una diagnosi.

Dipendenze: numpy, matplotlib; per l'interfaccia anche tkinter
(su macOS con Python di Homebrew: brew install python-tk@<versione>).

Autore:  Silvestro Scuderi  <https://github.com/silfox70>
Licenza: MIT (vedi il file LICENSE)
Sorgenti e documentazione del formato: https://github.com/silfox70/scpecg
"""

__author__ = "Silvestro Scuderi"
__license__ = "MIT"
__version__ = "1.0.0"
__url__ = "https://github.com/silfox70/scpecg"

import argparse
import datetime
import os
import struct
import sys

import numpy as np


# --------------------------------------------------------------------------
# Tabelle dello standard
# --------------------------------------------------------------------------

# Codici di identificazione delle derivazioni (SCP-ECG, sezione 3)
LEAD_NAMES = {
    0: "sconosciuta", 1: "I", 2: "II", 3: "V1", 4: "V2", 5: "V3",
    6: "V4", 7: "V5", 8: "V6", 9: "V7", 10: "V2R", 11: "V3R",
    12: "V4R", 13: "V5R", 14: "V6R", 15: "V7R", 16: "X", 17: "Y",
    18: "Z", 19: "CC5", 20: "CM5", 21: "left arm", 22: "right arm",
    23: "left leg", 24: "I1", 25: "E", 26: "C", 27: "A", 28: "M",
    29: "F", 30: "H", 31: "I-cal", 32: "II-cal", 33: "V1-cal",
    61: "III", 62: "aVR", 63: "aVL", 64: "aVF",
    # codice fuori standard usato da alcuni palmari a due elettrodi
    101: "palmare (~DI)",
}

# Tag della sezione 1 che ci interessa interpretare
SECTION1_TAGS = {
    0: "cognome paziente",
    1: "primo nome paziente",
    2: "ID paziente",
    3: "secondo nome paziente",
    4: "eta",
    5: "data di nascita",
    6: "altezza",
    7: "peso",
    8: "sesso",
    14: "ID apparecchio acquisizione",
    15: "ID apparecchio analisi",
    16: "istituzione",
    20: "medico richiedente",
    25: "data di acquisizione",
    26: "ora di acquisizione",
    27: "filtro passa-alto",
    28: "filtro passa-basso",
    29: "bitmap filtri",
    34: "fuso orario",
    255: "terminatore",
}

SEX = {0: "non noto", 1: "maschile", 2: "femminile", 9: "non specificato"}

# Marcatore convenzionale: "uso la tabella di Huffman di default"
HUFFMAN_DEFAULT_TABLE = 19999

SECTION_HEADER_LEN = 16

# valori convenzionali offerti dall'interfaccia
GAINS = [2.5, 5, 10, 20, 40]          # mm/mV
SPEEDS = [12.5, 25, 50]               # mm/s
WINDOWS = [2, 5, 10, 15, 30]          # secondi visibili


# --------------------------------------------------------------------------
# Utilita' di basso livello
# --------------------------------------------------------------------------

def u16(buf, off):
    return struct.unpack_from("<H", buf, off)[0]


def u32(buf, off):
    return struct.unpack_from("<I", buf, off)[0]


def crc_ccitt(data):
    """CRC-CCITT come definito dallo standard SCP-ECG (poly 0x1021, init 0xFFFF)."""
    crc = 0xFFFF
    for byte in data:
        crc ^= byte << 8
        for _ in range(8):
            if crc & 0x8000:
                crc = ((crc << 1) ^ 0x1021) & 0xFFFF
            else:
                crc = (crc << 1) & 0xFFFF
    return crc


# --------------------------------------------------------------------------
# Parsing delle sezioni
# --------------------------------------------------------------------------

class ScpFile:
    def __init__(self, path):
        self.path = path
        self.invert = False
        with open(path, "rb") as fh:
            self.raw = fh.read()

        if len(self.raw) < 8:
            raise ValueError("file troppo corto per essere un SCP-ECG")

        self.file_crc = u16(self.raw, 0)
        self.record_len = u32(self.raw, 2)

        self.sections = {}      # id -> dict(offset, length, version, protocol)
        self.demographics = {}
        self.leads = []
        self.rhythm = None      # dict con avm, interval_us, samples ...
        self.rpeak_flag_used = False
        self.segmenti = None       # popolato solo da merge_segments
        self.giunzioni = []

        self._parse_section0()
        self._parse_section1()
        self._parse_section2()
        self._parse_section3()
        self._parse_section6()

    # -- verifiche ---------------------------------------------------------

    def integrity(self):
        """Ritorna (crc_ok, len_ok, crc_calcolato, dimensione_reale)."""
        actual = len(self.raw)
        calc = crc_ccitt(self.raw[2:self.record_len]) if self.record_len <= actual else None
        return (calc == self.file_crc if calc is not None else False,
                self.record_len == actual, calc, actual)

    # -- sezione 0: mappa del file ----------------------------------------

    def _parse_section0(self):
        base = 6                       # gli offset SCP sono 1-based
        sec_id = u16(self.raw, base + 2)
        if sec_id != 0:
            raise ValueError("sezione 0 non trovata: non sembra un SCP-ECG")

        marker = self.raw[base + 10:base + 16].rstrip(b"\x00")
        if marker != b"SCPECG":
            print("attenzione: firma 'SCPECG' assente nell'header della sezione 0",
                  file=sys.stderr)

        length = u32(self.raw, base + 4)
        self.scp_version = self.raw[base + 8]
        self.scp_protocol = self.raw[base + 9]

        # i dati della sezione 0 sono una lista di puntatori da 10 byte
        ptr_area = self.raw[base + SECTION_HEADER_LEN:base + length]
        for i in range(0, len(ptr_area) - 9, 10):
            sid = u16(ptr_area, i)
            slen = u32(ptr_area, i + 2)
            soff = u32(ptr_area, i + 6)
            if slen == 0 or soff == 0:
                continue
            self.sections[sid] = {
                "offset": soff - 1,     # da 1-based a offset di file
                "length": slen,
            }

    def _section_data(self, sid):
        """Ritorna (bytes dei dati, dict header) per la sezione richiesta."""
        if sid not in self.sections:
            return None, None
        off = self.sections[sid]["offset"]
        length = self.sections[sid]["length"]
        head = {
            "crc": u16(self.raw, off),
            "id": u16(self.raw, off + 2),
            "length": u32(self.raw, off + 4),
            "version": self.raw[off + 8],
            "protocol": self.raw[off + 9],
        }
        data = self.raw[off + SECTION_HEADER_LEN:off + length]
        return data, head

    # -- sezione 1: anagrafica e acquisizione -----------------------------

    def _parse_section1(self):
        data, _ = self._section_data(1)
        if data is None:
            return
        i = 0
        while i + 3 <= len(data):
            tag = data[i]
            vlen = u16(data, i + 1)
            if tag == 255:
                break
            value = data[i + 3:i + 3 + vlen]
            self._store_tag(tag, value)
            i += 3 + vlen

    def _store_tag(self, tag, value):
        name = SECTION1_TAGS.get(tag, f"tag {tag}")
        if tag == 25 and len(value) >= 4:            # data
            y, m, d = u16(value, 0), value[2], value[3]
            self.demographics[name] = f"{y:04d}-{m:02d}-{d:02d}"
            self.demographics["_date_raw"] = (y, m, d)
        elif tag == 26 and len(value) >= 3:          # ora
            self.demographics[name] = f"{value[0]:02d}:{value[1]:02d}:{value[2]:02d}"
            self.demographics["_time_raw"] = (value[0], value[1], value[2])
        elif tag == 8 and len(value) >= 1:
            self.demographics[name] = SEX.get(value[0], str(value[0]))
        elif tag in (27, 28) and len(value) >= 2:
            self.demographics[name] = f"{u16(value, 0)} (unita' dello standard)"
        elif tag in (0, 1, 2, 3, 14, 15, 16, 20):
            txt = value.rstrip(b"\x00").decode("latin-1", "replace").strip()
            if txt:
                self.demographics[name] = txt
        elif value:
            self.demographics[name] = value.hex(" ")

    # -- sezione 2: tabelle di Huffman ------------------------------------

    def _parse_section2(self):
        self.huffman_tables = None
        self.huffman_default = False
        self.compressed = False
        data, _ = self._section_data(2)
        if data is None or len(data) < 2:
            return
        n = u16(data, 0)
        self.huffman_tables = n
        if n == HUFFMAN_DEFAULT_TABLE:
            self.huffman_default = True
            self.compressed = True
        elif n == 1 and len(data) >= 4 and u16(data, 2) == 1:
            # una tabella con una sola struttura di codice: convenzione per
            # "nessuna compressione effettiva"
            self.compressed = False
        else:
            self.compressed = True

    # -- sezione 3: definizione delle derivazioni -------------------------

    def _parse_section3(self):
        data, _ = self._section_data(3)
        if data is None or len(data) < 2:
            return
        n_leads = data[0]
        flags = data[1]
        self.leads_simultaneous = bool(flags & 0x04)
        self.n_simultaneous = (flags >> 3) & 0x1F
        for k in range(n_leads):
            base = 2 + k * 9
            if base + 9 > len(data):
                break
            start = u32(data, base)
            end = u32(data, base + 4)
            lid = data[base + 8]
            self.leads.append({
                "id": lid,
                "name": LEAD_NAMES.get(lid, f"lead {lid}"),
                "start": start,
                "end": end,
                "n": end - start + 1,
            })

    # -- sezione 6: dati di ritmo -----------------------------------------

    def _parse_section6(self):
        data, _ = self._section_data(6)
        if data is None or len(data) < 8:
            return

        avm = u16(data, 0)                 # nanovolt per unita'
        interval = u16(data, 2)            # microsecondi fra campioni
        diff = data[4]                     # 0=assoluti 1=diff prima 2=diff seconda
        bimodal = data[5]

        n_leads = max(1, len(self.leads))
        block_lengths = [u16(data, 6 + 2 * i) for i in range(n_leads)]

        payload_start = 6 + 2 * n_leads
        payload = data[payload_start:]

        if self.compressed:
            raise NotImplementedError(
                "la sezione 6 e' compressa con Huffman: questo script gestisce "
                "solo dati non compressi (decoder non implementato)"
            )

        signals = []
        rpeaks = []
        cursor = 0
        for blen in block_lengths:
            chunk = payload[cursor:cursor + blen]
            cursor += blen
            chunk = chunk[: (len(chunk) // 2) * 2]

            # Lo standard prescrive interi con segno. Alcuni apparecchi
            # (fuori standard) usano pero' il bit 15 come marcatore del picco
            # R rilevato a bordo: letti come int16 quei campioni diventano
            # enormi valori negativi. Li riconosco e li tratto a parte.
            raw_u = np.frombuffer(chunk, dtype="<u2").astype(np.int64)
            flagged = np.where(raw_u & 0x8000)[0]
            if self._looks_like_rpeak_flags(raw_u, flagged):
                self.rpeak_flag_used = True
                arr = (raw_u & 0x7FFF).astype(np.float64)
                rpeaks.append(flagged)
            else:
                arr = np.frombuffer(chunk, dtype="<i2").astype(np.float64)
                rpeaks.append(np.array([], dtype=int))

            if diff == 1:
                arr = np.cumsum(arr)
            elif diff == 2:
                arr = np.cumsum(np.cumsum(arr))
            signals.append(arr)

        self.rhythm = {
            "avm_nv": avm,
            "interval_us": interval,
            "fs": 1e6 / interval if interval else None,
            "diff": diff,
            "bimodal": bimodal,
            "block_lengths": block_lengths,
            "signals": signals,
            "rpeaks": rpeaks,
        }

    @staticmethod
    def _looks_like_rpeak_flags(raw_u, flagged):
        """Decide se il bit 15 e' un marcatore o il segno di un int16.

        Criterio: i campioni marcati devono essere pochi (un ECG con segnale
        davvero negativo ne avrebbe a migliaia) e, tolto il bit, i loro valori
        devono cadere nello stesso intervallo degli altri campioni. Se invece
        sono valori con segno legittimi, il mascheramento produce numeri
        scollegati dal resto del tracciato.
        """
        n = len(raw_u)
        if not len(flagged) or len(flagged) > 0.02 * n:
            return False
        plain = raw_u[(raw_u & 0x8000) == 0]
        if not len(plain):
            return False
        lo, hi = np.percentile(plain, [0.5, 99.5])
        span = max(hi - lo, 1.0)
        masked = raw_u[flagged] & 0x7FFF
        # devono stare entro il range dei campioni normali, allargato del 100%
        return bool(np.all(masked > lo - span) and np.all(masked < hi + span))

    def acquisition_datetime(self):
        """Istante di inizio della registrazione, dalla sezione 1.

        E' l'ora di INIZIO del segmento, non di fine: verificato su tre
        segmenti consecutivi distanti 30 s l'uno dall'altro.
        Ritorna None se il file non porta data e ora.
        """
        d = self.demographics.get("_date_raw")
        t = self.demographics.get("_time_raw")
        if not d or not t:
            return None
        try:
            return datetime.datetime(d[0], d[1], d[2], t[0], t[1], t[2])
        except ValueError:
            return None

    def duration_s(self):
        if not self.rhythm or not self.rhythm["signals"]:
            return 0.0
        return len(self.rhythm["signals"][0]) * self.rhythm["interval_us"] / 1e6

    def heart_rate(self, lead=0):
        """Frequenza dai marcatori di picco R, se presenti.

        Ritorna None se l'apparecchio non li ha scritti. La risoluzione e'
        limitata dall'intervallo di campionamento: a 150 Hz sono 6.7 ms per
        intervallo RR, sufficienti per la frequenza media ma non per
        un'analisi fine della variabilita'.
        """
        if not self.rhythm or not self.rhythm.get("rpeaks"):
            return None
        pk = self.rhythm["rpeaks"][lead]
        if len(pk) < 3:
            return None
        rr = np.diff(pk) * self.rhythm["interval_us"] / 1e6     # secondi
        bpm = 60.0 / rr
        return {
            "n_peaks": len(pk),
            "rr_mean_s": float(rr.mean()),
            "rr_min_s": float(rr.min()),
            "rr_max_s": float(rr.max()),
            "bpm_mean": float(bpm.mean()),
            "bpm_min": float(bpm.min()),
            "bpm_max": float(bpm.max()),
            "rr_sd_ms": float(rr.std() * 1000),
        }

    # -- conversione in millivolt -----------------------------------------

    def millivolts(self, baseline="auto"):
        """Ritorna lista di array in mV, uno per derivazione.

        baseline='auto'  toglie l'offset se i campioni sembrano codici ADC
                         non centrati sullo zero (es. offset 2048 a 12 bit)
        baseline='none'  lascia i valori come stanno
        baseline=<int>   toglie il valore indicato
        """
        if not self.rhythm:
            return []
        gain = self.rhythm["avm_nv"] / 1e6          # nV -> mV
        out = []
        for arr in self.rhythm["signals"]:
            a = arr
            if baseline == "auto":
                med = float(np.median(a))
                # se la mediana e' lontana da zero, e' un offset di mezza scala
                if abs(med) > 100:
                    a = a - med
            elif baseline == "none":
                pass
            else:
                a = a - float(baseline)
            out.append(-a * gain if self.invert else a * gain)
        return out

    def detected_offset(self):
        if not self.rhythm or not self.rhythm["signals"]:
            return None
        return float(np.median(self.rhythm["signals"][0]))


# --------------------------------------------------------------------------
# Stampa della struttura
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# Unione di segmenti consecutivi
#
# In misura continua il PC-80B non salva un file unico: spezza la
# registrazione in record da 30 secondi. Analizzarli separatamente perde
# l'intervallo RR a cavallo di ogni taglio e fa ripartire da uno il conteggio
# dei battiti. Qui li si rimette insieme.
#
# L'ordinamento viene dai timestamp della sezione 1, non dal nome del file
# (l'ordine alfabetico metterebbe 10 prima di 2). La contiguita' e' verificata:
# se fra due segmenti manca un pezzo, ci si ferma invece di incollare
# comunque, perche' un buco silenzioso produrrebbe un intervallo RR falso e
# quindi un'anomalia inventata.
# --------------------------------------------------------------------------

# il timestamp ha risoluzione di un secondo, la durata reale di un segmento e'
# 29.997 s: serve una tolleranza, altrimenti ogni giunzione sembrerebbe rotta
TOLLERANZA_GIUNZIONE_S = 1.5


class SegmentiNonContigui(Exception):
    pass


def merge_segments(paths):
    """Unisce piu' file .scp consecutivi in un unico oggetto analizzabile.

    Solleva SegmentiNonContigui se manca un segmento o se i file non sono
    omogenei. L'oggetto restituito ha in piu' due attributi:
      .segmenti   lista di (nome file, istante di inizio, n campioni)
      .giunzioni  indici dei campioni in cui finisce un segmento
    """
    if len(paths) == 1:
        return ScpFile(paths[0])

    scps = []
    for p in paths:
        s = ScpFile(p)
        if not s.rhythm:
            raise SegmentiNonContigui(f"{p}: nessun dato di ritmo")
        if s.acquisition_datetime() is None:
            raise SegmentiNonContigui(
                f"{p}: manca data o ora nella sezione 1, impossibile ordinare "
                "i segmenti in modo affidabile")
        scps.append(s)

    scps.sort(key=lambda s: s.acquisition_datetime())

    rif = scps[0].rhythm
    for s in scps[1:]:
        if s.rhythm["interval_us"] != rif["interval_us"]:
            raise SegmentiNonContigui(
                f"{os.path.basename(s.path)}: frequenza di campionamento "
                "diversa dagli altri segmenti")
        if s.rhythm["avm_nv"] != rif["avm_nv"]:
            raise SegmentiNonContigui(
                f"{os.path.basename(s.path)}: AVM diverso dagli altri segmenti")
        if len(s.rhythm["signals"]) != len(rif["signals"]):
            raise SegmentiNonContigui(
                f"{os.path.basename(s.path)}: numero di derivazioni diverso")

    # contiguita': l'inizio di ogni segmento deve cadere alla fine del precedente
    for a, b in zip(scps, scps[1:]):
        atteso = a.duration_s()
        reale = (b.acquisition_datetime() - a.acquisition_datetime()).total_seconds()
        scarto = reale - atteso
        if abs(scarto) > TOLLERANZA_GIUNZIONE_S:
            manca = scarto
            raise SegmentiNonContigui(
                f"sequenza incompleta fra {os.path.basename(a.path)} "
                f"({a.acquisition_datetime():%H:%M:%S}) e "
                f"{os.path.basename(b.path)} "
                f"({b.acquisition_datetime():%H:%M:%S}): "
                f"mancano {manca:.0f} s di registrazione.\n"
                "Recupera i segmenti mancanti, oppure uniscine un "
                "sottoinsieme davvero consecutivo.")

    base = scps[0]
    n_leads = len(rif["signals"])
    segnali, picchi, giunzioni, segmenti = [], [], [], []
    offset = 0
    for li in range(n_leads):
        segnali.append([])
        picchi.append([])
    for s in scps:
        n = len(s.rhythm["signals"][0])
        for li in range(n_leads):
            segnali[li].append(s.rhythm["signals"][li])
            pk = s.rhythm["rpeaks"][li] if li < len(s.rhythm["rpeaks"]) else []
            picchi[li].append(np.asarray(pk, dtype=int) + offset)
        segmenti.append((os.path.basename(s.path), s.acquisition_datetime(), n))
        offset += n
        giunzioni.append(offset)

    base.rhythm["signals"] = [np.concatenate(x) for x in segnali]
    base.rhythm["rpeaks"] = [np.concatenate(x) if any(len(y) for y in x)
                             else np.array([], dtype=int) for x in picchi]
    base.segmenti = segmenti
    base.giunzioni = giunzioni[:-1]      # l'ultimo confine e' la fine del file
    base.rpeak_flag_used = any(s.rpeak_flag_used for s in scps)
    return base


# --------------------------------------------------------------------------
# Analisi del ritmo
#
# PROVENIENZA DELLE SOGLIE. Il manuale del PC-80B elenca i 17 verdetti ma non
# pubblica i criteri numerici. Qui sotto ogni soglia dice da dove viene:
#
#   [manuale]   regola definita nella documentazione Creative/Gima
#   [nostra]    valore scelto da noi, plausibile ma NON quello dell'apparecchio
#
# Cambiando questi numeri si ritara tutto: sono l'unico punto da toccare se un
# giorno le soglie vere saltano fuori.
# --------------------------------------------------------------------------

SOGLIE = {
    # [manuale] 60-100 bpm e' il ritmo normale; sotto e' lento, sopra e' rapido
    "bpm_lento": 60.0,
    "bpm_rapido": 100.0,

    # [nostra] fasce intermedie "a little slow" / "a little fast": il PC-80B le
    # riporta ma non le definisce
    "bpm_poco_lento": 50.0,
    "bpm_poco_rapido": 120.0,

    # [nostra] quanto un intervallo deve essere piu' corto della mediana per
    # contare come battito anticipato
    "anticipo_frazione": 0.20,

    # [manuale] battito mancante: intervallo pari al doppio della media dei
    # precedenti. Uso 1.8 per lasciare margine alla variabilita' respiratoria
    "mancante_fattore": 1.8,

    # [nostra] irregolarita': coefficiente di variazione degli intervalli RR
    "cv_irregolare": 0.10,

    # [manuale] short run: extrasistoli consecutive piu' di tre volte
    "short_run_min": 3,

    # [nostra] deriva della linea di base, in mV picco-picco della componente
    # lenta
    "deriva_mv": 0.15,

    # [manuale] rumore interno dello strumento <= 30 uVpp; oltre questa soglia
    # [nostra] il segnale e' considerato inutilizzabile
    "rumore_uv_poor": 120.0,
}


def signal_quality(scp, lead=0):
    """Indici di qualita': rumore ad alta frequenza e deriva della linea di base.

    Il rumore lo stimo come scarto fra il segnale e la sua versione lisciata su
    tre campioni: cattura la componente veloce senza toccare il QRS, che e'
    coerente e sopravvive alla lisciatura. La deriva la stimo con una media
    mobile lunga un secondo, che tiene solo la componente respiratoria e di
    contatto.
    """
    sigs = scp.millivolts()
    if not sigs or lead >= len(sigs):
        return None
    x = sigs[lead]
    fs = scp.rhythm["fs"]

    liscio = np.convolve(x, np.ones(3) / 3.0, mode="same")
    rumore_uv = float(np.std(x - liscio) * 1000.0)

    w = max(3, int(round(fs)))
    deriva = np.convolve(x, np.ones(w) / w, mode="valid")
    deriva_mv = float(np.percentile(deriva, 99) - np.percentile(deriva, 1))

    return {
        "rumore_uv": rumore_uv,
        "deriva_mv": deriva_mv,
        "wander": deriva_mv > SOGLIE["deriva_mv"],
        "poor": rumore_uv > SOGLIE["rumore_uv_poor"],
    }


def analyze_rhythm(scp, lead=0):
    """Classifica gli intervalli RR e compone un verdetto in stile PC-80B.

    Lavora sui marcatori di picco R scritti dall'apparecchio. Non li ricalcola:
    se l'apparecchio non li ha scritti (misura a palmo) l'analisi non parte.
    """
    if not scp.rhythm or not scp.rhythm.get("rpeaks"):
        return None
    pk = scp.rhythm["rpeaks"][lead] if lead < len(scp.rhythm["rpeaks"]) else []
    if len(pk) < 4:
        return None

    pk = np.asarray(pk)
    t = pk / scp.rhythm["fs"]
    rr = np.diff(t)

    # Un intervallo che scavalca il confine fra due segmenti non e' misurato:
    # fra la fine di un file e l'inizio del successivo l'apparecchio puo' aver
    # perso qualche campione, e classificarlo produrrebbe un'anomalia inventata
    # dal montaggio. Lo si marca e lo si esclude dalle statistiche.
    su_giunzione = np.zeros(len(rr), dtype=bool)
    for g in getattr(scp, "giunzioni", []) or []:
        su_giunzione |= (pk[:-1] < g) & (pk[1:] >= g)

    buoni = rr[~su_giunzione]
    if len(buoni) < 2:
        buoni = rr
    mediana = float(np.median(buoni))
    bpm_medio = 60.0 / float(np.mean(buoni))
    cv = float(np.std(buoni) / np.mean(buoni))

    # classificazione intervallo per intervallo
    etichette = []
    for i, v in enumerate(rr):
        if su_giunzione[i]:
            etichette.append("giunzione")
        elif v > SOGLIE["mancante_fattore"] * mediana:
            etichette.append("mancante")
        elif v < (1.0 - SOGLIE["anticipo_frazione"]) * mediana:
            etichette.append("anticipato")
        else:
            etichette.append("normale")

    anticipati = [i for i, e in enumerate(etichette) if e == "anticipato"]
    mancanti = [i for i, e in enumerate(etichette) if e == "mancante"]
    giunzioni_rr = [i for i, e in enumerate(etichette) if e == "giunzione"]

    # sequenze di anticipi consecutivi: short run
    run_max, run = 0, 0
    for e in etichette:
        run = run + 1 if e == "anticipato" else 0
        run_max = max(run_max, run)

    # bigeminismo e trigeminismo: anticipi a passo regolare di 2 o 3
    def passo_regolare(idx, passo):
        if len(idx) < 3:
            return False
        d = np.diff(idx)
        return bool(np.all(d == passo))

    bigeminismo = passo_regolare(anticipati, 2)
    trigeminismo = passo_regolare(anticipati, 3)

    # asse frequenza
    if bpm_medio < SOGLIE["bpm_poco_lento"]:
        freq = "lento"
    elif bpm_medio < SOGLIE["bpm_lento"]:
        freq = "poco lento"
    elif bpm_medio <= SOGLIE["bpm_rapido"]:
        freq = "normale"
    elif bpm_medio <= SOGLIE["bpm_poco_rapido"]:
        freq = "poco rapido"
    else:
        freq = "rapido"

    # asse regolarita'
    if run_max >= SOGLIE["short_run_min"]:
        reg = "serie di battiti ravvicinati"
    elif bigeminismo:
        reg = "sospetto bigeminismo"
    elif trigeminismo:
        reg = "sospetto trigeminismo"
    elif mancanti:
        reg = "battiti mancanti"
    elif cv > SOGLIE["cv_irregolare"]:
        reg = "intervalli irregolari"
    elif anticipati:
        reg = "battiti anticipati occasionali"
    else:
        reg = "regolare"

    q = signal_quality(scp, lead) or {}

    # composizione del messaggio, nello stile della tabella del manuale
    pezzi = []
    if freq != "normale":
        pezzi.append({"lento": "battito lento",
                      "poco lento": "battito lievemente lento",
                      "poco rapido": "battito lievemente rapido",
                      "rapido": "battito rapido"}[freq])
    if reg != "regolare":
        pezzi.append(reg)
    if q.get("wander"):
        pezzi.append("deriva della linea di base")

    if q.get("poor"):
        verdetto = "Segnale scarso, ripetere la misura"
    elif not pezzi:
        verdetto = "Nessuna irregolarita' evidente"
    else:
        verdetto = "Sospetto: " + ", ".join(pezzi)

    return {
        "bpm_medio": bpm_medio,
        "bpm_min": float(60.0 / buoni.max()),
        "bpm_max": float(60.0 / buoni.min()),
        "rr_mediana_s": mediana,
        "cv": cv,
        "n_battiti": len(t),
        "etichette": etichette,
        "tempi": t,
        "anticipati": anticipati,
        "mancanti": mancanti,
        "giunzioni_rr": giunzioni_rr,
        "run_max": run_max,
        "qualita": q,
        "verdetto": verdetto,
    }


def print_analysis(scp):
    a = analyze_rhythm(scp)
    q = signal_quality(scp)

    if q:
        print("\nQualita' del segnale")
        print(f"  Rumore alta frequenza : {q['rumore_uv']:.0f} uV "
              f"(rumore interno dichiarato dello strumento: 30 uVpp)")
        print(f"  Deriva linea di base  : {q['deriva_mv']:.3f} mV"
              + ("   [oltre soglia]" if q["wander"] else ""))

    if a is None:
        print("\nAnalisi del ritmo non eseguita: il file non contiene i "
              "marcatori\ndei picchi R (tipico della misura a palmo, dove "
              "l'apparecchio\nfornisce il proprio verdetto sul display).")
        return

    print("\nAnalisi del ritmo")
    print(f"  Battiti               : {a['n_battiti']}")
    print(f"  Frequenza media       : {a['bpm_medio']:.1f} bpm "
          f"(da {a['bpm_min']:.0f} a {a['bpm_max']:.0f})")
    print(f"  Coeff. variazione RR  : {a['cv']*100:.1f} %")
    if a["anticipati"]:
        istanti = ", ".join(f"{a['tempi'][i+1]:.1f}s" for i in a["anticipati"][:8])
        print(f"  Intervalli anticipati : {len(a['anticipati'])}  a {istanti}"
              + (" ..." if len(a["anticipati"]) > 8 else ""))
    if a.get("giunzioni_rr"):
        print(f"  Intervalli su giunzione: {len(a['giunzioni_rr'])}  "
              "(non classificati, cadono fra due segmenti)")
    if a["mancanti"]:
        istanti = ", ".join(f"{a['tempi'][i+1]:.1f}s" for i in a["mancanti"][:8])
        print(f"  Intervalli lunghi     : {len(a['mancanti'])}  a {istanti}")
    print(f"\n  >>> {a['verdetto']}")
    print("\n  Verdetto calcolato da questo programma con soglie in parte "
          "scelte da noi\n  (vedi il dizionario SOGLIE): NON e' il verdetto "
          "dell'apparecchio e non\n  e' una diagnosi. L'interpretazione "
          "spetta al medico.")



# --------------------------------------------------------------------------
# Simulatore: crea copie alterate di un file reale per collaudare l'analisi
#
# Serve a verificare che il classificatore reagisca come previsto: si parte da
# un tracciato vero, si inietta un'anomalia nota e si guarda cosa dice.
# L'originale non viene MAI modificato: si scrive sempre un file nuovo, con il
# nome che dichiara l'alterazione.
# --------------------------------------------------------------------------

SIMULAZIONI = ("anticipo", "mancante", "deriva", "rumore", "rapido", "lento")


def _span_sezione6(raw):
    """Ritorna (offset sezione, lunghezza, offset primo campione, n campioni)."""
    off6 = ln6 = None
    for i in range(22, 22 + 60, 10):
        sid, ln, off = u16(raw, i), u32(raw, i + 2), u32(raw, i + 6)
        if sid == 6 and ln:
            off6, ln6 = off - 1, ln
    if off6 is None:
        raise ValueError("sezione 6 assente")
    n_leads = 1
    dati = off6 + SECTION_HEADER_LEN
    blk = u16(raw, dati + 6)
    return off6, ln6, dati + 6 + 2 * n_leads, blk // 2


def simulate(path, kind, out=None):
    """Scrive una copia di `path` con l'anomalia `kind` iniettata."""
    if kind not in SIMULAZIONI:
        raise ValueError(f"simulazione sconosciuta: {kind}")

    raw = bytearray(open(path, "rb").read())
    off6, ln6, start, n = _span_sezione6(raw)

    u = np.frombuffer(bytes(raw[start:start + 2 * n]), dtype="<u2").astype(np.int64)
    flag = (u & 0x8000) != 0
    val = (u & 0x7FFF).astype(np.float64)
    base = float(np.median(val))
    pk = np.where(flag)[0]

    if kind in ("anticipo", "mancante") and len(pk) < 5:
        raise ValueError("servono i marcatori dei picchi R per questa "
                         "simulazione (file registrato a palmo?)")

    if kind == "anticipo":
        # deformazione locale dell'asse dei tempi: un battito viene spostato
        # in avanti e il successivo si allontana, come in una extrasistole
        # con pausa compensatoria. La forma del QRS si comprime un poco:
        # e' una simulazione, non un vero battito ectopico.
        k = len(pk) // 2
        a, c = int(pk[k - 1]), int(pk[k + 1])
        p = int(pk[k])
        delta = int(0.28 * (p - a))
        q = p - delta
        dest = np.arange(a, c + 1)
        src = np.where(
            dest <= q,
            a + (dest - a) * (p - a) / max(q - a, 1),
            p + (dest - q) * (c - p) / max(c - q, 1))
        val[a:c + 1] = np.interp(src, np.arange(len(val)), val)
        flag[p] = False
        flag[q] = True

    elif kind == "mancante":
        # un battito viene proprio tolto: campioni riportati alla linea di
        # base e marcatore rimosso. L'intervallo risultante e' doppio.
        k = len(pk) // 2
        p = int(pk[k])
        a = max(0, p - int(0.20 * (p - pk[k - 1])))
        c = min(len(val) - 1, p + int(0.55 * (pk[k + 1] - p)))
        val[a:c + 1] = np.linspace(val[a], val[c], c - a + 1)
        flag[a:c + 1] = False

    elif kind == "deriva":
        t = np.arange(n) / n
        val += 250.0 * np.sin(2 * np.pi * 6.0 * t)      # ~0.2 mV di oscillazione lenta

    elif kind == "rumore":
        rng = np.random.default_rng(0)
        val += rng.normal(0, 190.0, n)                  # ~0.15 mV RMS

    elif kind in ("rapido", "lento"):
        # ricampionamento dell'intero tracciato: comprimere l'asse dei tempi
        # aumenta la frequenza, allungarlo la riduce. Il numero di campioni
        # resta 4500, quindi il file mantiene la stessa struttura.
        f = 1.5 if kind == "rapido" else 0.46
        v2 = np.tile(val, 3)
        src = np.arange(n) * f
        val = np.interp(src, np.arange(len(v2)), v2)
        nuovo = np.zeros(n, dtype=bool)
        for p in np.where(np.tile(flag, 3))[0]:
            j = int(round(p / f))
            if 0 <= j < n:
                nuovo[j] = True
        flag = nuovo

    val = np.clip(np.round(val - np.median(val) + base), 0, 4095).astype(np.int64)
    u2 = (val | np.where(flag, 0x8000, 0)).astype("<u2")
    raw[start:start + 2 * n] = u2.tobytes()

    # CRC della sezione 6, poi CRC del record: l'ordine conta
    raw[off6:off6 + 2] = struct.pack("<H", crc_ccitt(bytes(raw[off6 + 2:off6 + ln6])))
    rec_len = u32(raw, 2)
    raw[0:2] = struct.pack("<H", crc_ccitt(bytes(raw[2:rec_len])))

    if out is None:
        stem, ext = os.path.splitext(path)
        out = f"{stem}_sim_{kind}{ext}"
    if os.path.abspath(out) == os.path.abspath(path):
        raise ValueError("il file di uscita coincide con l'originale")
    with open(out, "wb") as fh:
        fh.write(raw)
    print(f"simulazione '{kind}' scritta in: {out}")
    return out


def print_segments(scp):
    if not getattr(scp, "segmenti", None):
        return
    tot = sum(n for _, _, n in scp.segmenti) * scp.rhythm["interval_us"] / 1e6
    print(f"\nSegmenti uniti: {len(scp.segmenti)}  ->  {tot:.1f} s complessivi")
    for nome, quando, n in scp.segmenti:
        print(f"  {nome:<16} {quando:%Y-%m-%d %H:%M:%S}  {n} campioni")
    print("  contiguita' verificata sui timestamp della sezione 1")


def print_info(scp):
    crc_ok, len_ok, calc, actual = scp.integrity()
    print(f"File           : {scp.path}")
    print(f"Dimensione     : {actual} byte")
    print(f"Record dichiar.: {scp.record_len} byte  "
          f"[{'ok' if len_ok else 'NON CORRISPONDE'}]")
    print(f"CRC header     : 0x{scp.file_crc:04X}  calcolato 0x{calc:04X}  "
          f"[{'ok' if crc_ok else 'NON CORRISPONDE'}]"
          if calc is not None else "CRC        : non verificabile")
    print(f"Versione SCP   : {scp.scp_version / 10:.1f} "
          f"(protocollo {scp.scp_protocol / 10:.1f})")

    print("\nSezioni presenti")
    print(f"  {'id':>3}  {'lunghezza':>10}  {'offset':>8}")
    for sid in sorted(scp.sections):
        s = scp.sections[sid]
        print(f"  {sid:>3}  {s['length']:>10}  {s['offset'] + 1:>8}")

    if scp.demographics:
        print("\nSezione 1 - anagrafica e acquisizione")
        for k, v in scp.demographics.items():
            if k.startswith("_"):
                continue
            print(f"  {k:<28}: {v}")

    print("\nSezione 2 - compressione")
    if scp.huffman_tables is None:
        print("  assente")
    elif scp.huffman_default:
        print(f"  tabella di Huffman di default ({HUFFMAN_DEFAULT_TABLE})")
    else:
        print(f"  {scp.huffman_tables} tabella/e dichiarate -> "
              f"{'compressione attiva' if scp.compressed else 'nessuna compressione'}")

    if scp.leads:
        print("\nSezione 3 - derivazioni")
        for ld in scp.leads:
            print(f"  {ld['name']:<14} (codice {ld['id']:>3})  "
                  f"campioni {ld['start']}-{ld['end']}  ({ld['n']})")
            if ld["id"] > 64:
                print("      codice fuori dalla tabella standard SCP-ECG "
                      "(che arriva a 64): valore proprietario del costruttore")

    r = scp.rhythm
    if r:
        n = len(r["signals"][0]) if r["signals"] else 0
        dur = n / r["fs"] if r["fs"] else 0
        print("\nSezione 6 - dati di ritmo")
        print(f"  AVM                 : {r['avm_nv']} nV per unita'")
        print(f"  Intervallo campioni : {r['interval_us']} us "
              f"-> {r['fs']:.2f} Hz")
        print(f"  Codifica differenze : {r['diff']} "
              f"({'valori assoluti' if r['diff'] == 0 else 'differenze'})")
        print(f"  Bimodale            : {r['bimodal']}")
        print(f"  Blocchi (byte)      : {r['block_lengths']}")
        print(f"  Campioni            : {n}  -> {dur:.2f} s")
        if scp.rpeak_flag_used:
            print("  Marcatori picco R   : presenti (bit 15 dei campioni, "
                  "convenzione fuori standard)")
        off = scp.detected_offset()
        if off is not None and abs(off) > 100:
            # l'offset e' tipicamente meta' scala: cerco la potenza di due piu'
            # vicina, senza pretendere che la mediana ci caschi esatta
            bits = min(range(8, 25),
                       key=lambda b: abs(abs(off) - (1 << (b - 1))))
            print(f"  Offset ADC rilevato : {off:.0f} "
                  f"(mezza scala di un ADC a {bits} bit = {1 << (bits - 1)}; "
                  f"fondo scala {(1 << bits) * r['avm_nv'] / 1e6:.2f} mV)")


# --------------------------------------------------------------------------
# Uscite
# --------------------------------------------------------------------------

def print_rate(scp):
    hr = scp.heart_rate()
    if not hr:
        return
    print("\nFrequenza dai marcatori dell'apparecchio")
    print(f"  Battiti marcati     : {hr['n_peaks']}")
    print(f"  Intervallo RR medio : {hr['rr_mean_s']:.3f} s "
          f"(min {hr['rr_min_s']:.3f}, max {hr['rr_max_s']:.3f})")
    print(f"  Frequenza media     : {hr['bpm_mean']:.1f} bpm "
          f"(min {hr['bpm_min']:.1f}, max {hr['bpm_max']:.1f})")
    print(f"  Deviazione std RR   : {hr['rr_sd_ms']:.1f} ms")
    res = scp.rhythm["interval_us"] / 1000.0
    print(f"  Nota: risoluzione temporale {res:.1f} ms per campione; adeguata")
    print("        alla frequenza media, non a un'analisi fine della "
          "variabilita'.")
    print("        Il rilevamento dei picchi e' quello dell'apparecchio, "
          "non ricalcolato.")


def write_csv(scp, path, baseline):
    sigs = scp.millivolts(baseline)
    if not sigs:
        print("nessun dato da esportare", file=sys.stderr)
        return
    fs = scp.rhythm["fs"]
    n = max(len(s) for s in sigs)
    names = [ld["name"] for ld in scp.leads] or [f"lead{i+1}" for i in range(len(sigs))]
    with open(path, "w") as fh:
        fh.write("tempo_s," + ",".join(f"{nm}_mV" for nm in names) + "\n")
        for i in range(n):
            row = [f"{i / fs:.6f}"]
            for s in sigs:
                row.append(f"{s[i]:.4f}" if i < len(s) else "")
            fh.write(",".join(row) + "\n")
    print(f"CSV scritto: {path}  ({n} righe)")


def choose_gain(sigs):
    """Sceglie un guadagno standard (10/20/40 mm/mV) in base all'ampiezza.

    Gli ECG clinici usano 10 mm/mV. Con derivazioni palmari il segnale e' cosi'
    piccolo da risultare illeggibile e si sale a 20 o 40; con elettrodi veri
    puo' invece uscire dalla carta e si scende a 5 o 2.5. Sono gli stessi
    valori del tasto del guadagno in reparto, quindi il tracciato resta
    misurabile col righello sapendo la scala (che finisce nel titolo).
    """
    peak = max(float(np.percentile(np.abs(s), 99.9)) for s in sigs)
    # dal guadagno piu' alto al piu' basso: prendo il primo che ci sta
    for gain in (40.0, 20.0, 10.0, 5.0, 2.5):
        if peak * gain <= 12.0:      # sta in +/- 12 mm di finestra
            return gain
    return 2.5


def write_plot(scp, path, baseline, seconds_per_row, mm_per_s=25.0, mm_per_mv=0.0):
    # Figura costruita senza pyplot: cosi' la stampa su file convive con
    # l'interfaccia Tk nello stesso processo, senza cambiare backend.
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.ticker import FuncFormatter, MultipleLocator

    # la griglia principale sta ogni 0.2 s, ma un'etichetta ogni 0.2 s
    # renderebbe l'asse illeggibile: ne stampo una al secondo
    def one_label_per_second(x, _pos):
        return f"{x:.0f}" if abs(x - round(x)) < 1e-6 else ""

    sigs = scp.millivolts(baseline)
    if not sigs:
        print("nessun dato da disegnare", file=sys.stderr)
        return

    fs = scp.rhythm["fs"]
    names = [ld["name"] for ld in scp.leads] or [f"lead{i+1}" for i in range(len(sigs))]

    if not mm_per_mv:
        mm_per_mv = choose_gain(sigs)
        if mm_per_mv != 10.0:
            motivo = ("segnale piccolo, scala amplificata" if mm_per_mv > 10
                      else "segnale ampio, scala ridotta per non tagliarlo")
            print(f"guadagno automatico: {mm_per_mv:g} mm/mV ({motivo}; "
                  f"lo standard clinico e' 10)")

    # mezza finestra verticale in mV: 15 mm di carta al guadagno scelto
    half_mv = 15.0 / mm_per_mv

    # rete di sicurezza: se il segnale non ci sta comunque (guadagno imposto a
    # mano, o ampiezza fuori scala), allargo la finestra invece di tagliare.
    # Un tracciato troncato in silenzio e' peggio di un grafico brutto.
    true_peak = max(float(np.max(np.abs(s))) for s in sigs)
    if true_peak > half_mv:
        half_mv = true_peak * 1.1
        print(f"finestra verticale allargata a +/-{half_mv:.2f} mV "
              f"per non troncare il tracciato "
              f"(picco reale {true_peak:.2f} mV)")

    # una riga per ogni spezzone di seconds_per_row secondi, per ogni derivazione
    peaks_by_lead = scp.rhythm.get("rpeaks") or [np.array([], dtype=int)] * len(sigs)

    strips = []
    for li, (sig, nm) in enumerate(zip(sigs, names)):
        total = len(sig) / fs
        n_rows = int(np.ceil(total / seconds_per_row))
        for r in range(n_rows):
            a = int(r * seconds_per_row * fs)
            b = min(len(sig), int((r + 1) * seconds_per_row * fs))
            pk = peaks_by_lead[li] if li < len(peaks_by_lead) else np.array([], dtype=int)
            pk_here = pk[(pk >= a) & (pk < b)] - a
            strips.append((nm, r * seconds_per_row, sig[a:b], pk_here))

    # dimensioni fisiche reali: carta ECG vera
    width_mm = seconds_per_row * mm_per_s
    height_mm = 30.0                       # 30 mm di carta, sempre
    inch = 1 / 25.4
    fig = Figure(figsize=(width_mm * inch + 1.2,
                          len(strips) * height_mm * inch + 1.0))
    FigureCanvasAgg(fig)
    axes = fig.subplots(len(strips), 1, squeeze=False)[:, 0]

    for ax, (nm, t0, seg, pk) in zip(axes, strips):
        t = t0 + np.arange(len(seg)) / fs

        # griglia: 1 mm = 0.04 s e 0.1 mV ; 5 mm = 0.2 s e 0.5 mV
        ax.xaxis.set_minor_locator(MultipleLocator(0.04))
        ax.xaxis.set_major_locator(MultipleLocator(0.2))
        ax.xaxis.set_major_formatter(FuncFormatter(one_label_per_second))
        # 1 mm e 5 mm di carta, tradotti in mV secondo il guadagno scelto
        ax.yaxis.set_minor_locator(MultipleLocator(1.0 / mm_per_mv))
        ax.yaxis.set_major_locator(MultipleLocator(5.0 / mm_per_mv))
        ax.grid(which="minor", color="#f4b8b8", linewidth=0.4)
        ax.grid(which="major", color="#e06666", linewidth=0.8)

        ax.plot(t, seg, color="#101010", linewidth=0.8)
        if len(pk):
            # tacca sotto ogni picco R segnalato dall'apparecchio
            ax.plot(t0 + pk / fs, np.full(len(pk), -half_mv * 0.85),
                    marker="|", linestyle="none", color="#1f6fb2",
                    markersize=5, markeredgewidth=1.0)
        ax.set_xlim(t0, t0 + seconds_per_row)
        ax.set_ylim(-half_mv, half_mv)
        ax.set_aspect((1.0 / mm_per_s) / (1.0 / mm_per_mv))   # 25 mm/s x 10 mm/mV
        short = nm if len(nm) <= 14 else nm[:13] + "."
        ax.set_ylabel(f"{short}\n{t0:.0f}s", fontsize=7)
        ax.tick_params(labelsize=6)
        for side in ax.spines.values():
            side.set_visible(False)

    dur = len(sigs[0]) / fs
    fig.suptitle(
        f"{os.path.basename(scp.path)} - {fs:.1f} Hz, {dur:.1f} s, "
        f"{mm_per_s:.0f} mm/s, {mm_per_mv:.0f} mm/mV",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=200, facecolor="white")
    print(f"Grafico scritto: {path}")


# --------------------------------------------------------------------------

# ==========================================================================
# Interfaccia grafica
#
# Tkinter e il backend TkAgg sono importati qui dentro e non in testa al file:
# cosi' la parte batch continua a funzionare dove Tk non c'e'.
# ==========================================================================

def run_gui(path=None):
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk

    import matplotlib
    matplotlib.use("TkAgg")
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib.figure import Figure
    from matplotlib.ticker import FuncFormatter, MultipleLocator
    from matplotlib.widgets import Cursor

    class ScpViewer(tk.Tk):

        def __init__(self, path=None):
            super().__init__()
            self.title("scpview - visualizzatore SCP-ECG")
            self.geometry("1180x680")

            self.scp = None
            self.sig = None            # segnale corrente in mV
            self.fs = None
            self.peaks = np.array([], dtype=int)
            self.lead_idx = 0

            self.analisi = None
            self.caliper = []          # 0, 1 o 2 punti (t, mV)
            self._caliper_art = []

            self._build_toolbar()
            self._build_verdict()
            self._build_canvas()
            self._build_statusbar()

            if path:
                self.load(path)
            else:
                self.after(200, self.on_open)

        # ------------------------------------------------------------------ UI

        def _build_toolbar(self):
            bar = ttk.Frame(self, padding=(8, 6))
            bar.pack(side=tk.TOP, fill=tk.X)

            ttk.Button(bar, text="Apri...", command=self.on_open).pack(side=tk.LEFT)
            ttk.Button(bar, text="Salva vista", command=self.on_save).pack(
                side=tk.LEFT, padx=(4, 0))
            ttk.Button(bar, text="Esporta CSV", command=self.on_export_csv).pack(
                side=tk.LEFT, padx=(4, 0))
            ttk.Button(bar, text="Esporta pagina", command=self.on_export_page).pack(
                side=tk.LEFT, padx=(4, 12))

            ttk.Label(bar, text="Derivazione").pack(side=tk.LEFT)
            self.var_lead = tk.StringVar()
            self.cmb_lead = ttk.Combobox(bar, textvariable=self.var_lead, width=16,
                                         state="readonly", values=[])
            self.cmb_lead.pack(side=tk.LEFT, padx=(4, 12))
            self.cmb_lead.bind("<<ComboboxSelected>>", lambda e: self.on_lead())

            ttk.Label(bar, text="Guadagno mm/mV").pack(side=tk.LEFT)
            self.var_gain = tk.StringVar(value="auto")
            cmb_g = ttk.Combobox(bar, textvariable=self.var_gain, width=6,
                                 state="readonly",
                                 values=["auto"] + [f"{g:g}" for g in GAINS])
            cmb_g.pack(side=tk.LEFT, padx=(4, 12))
            cmb_g.bind("<<ComboboxSelected>>", lambda e: self.redraw())

            ttk.Label(bar, text="Velocita' mm/s").pack(side=tk.LEFT)
            self.var_speed = tk.StringVar(value="25")
            cmb_s = ttk.Combobox(bar, textvariable=self.var_speed, width=6,
                                 state="readonly",
                                 values=[f"{s:g}" for s in SPEEDS])
            cmb_s.pack(side=tk.LEFT, padx=(4, 12))
            cmb_s.bind("<<ComboboxSelected>>", lambda e: self.redraw())

            ttk.Label(bar, text="Finestra s").pack(side=tk.LEFT)
            self.var_win = tk.StringVar(value="10")
            cmb_w = ttk.Combobox(bar, textvariable=self.var_win, width=5,
                                 state="readonly",
                                 values=[str(w) for w in WINDOWS])
            cmb_w.pack(side=tk.LEFT, padx=(4, 12))
            cmb_w.bind("<<ComboboxSelected>>", lambda e: self.on_window())

            self.var_marks = tk.BooleanVar(value=True)
            ttk.Checkbutton(bar, text="Battiti", variable=self.var_marks,
                            command=self.redraw).pack(side=tk.LEFT)

            # Con le proporzioni reali il rapporto fra assi rispetta davvero
            # mm/s e mm/mV: giusto per la stampa, ma sullo schermo lascia il
            # tracciato schiacciato in una striscia. Di default riempio la
            # finestra; le misure restano esatte perche' si leggono in numeri.
            self.var_aspect = tk.BooleanVar(value=False)
            ttk.Checkbutton(bar, text="Proporzioni reali", variable=self.var_aspect,
                            command=self.redraw).pack(side=tk.LEFT, padx=(8, 0))

            ttk.Button(bar, text="Azzera caliper",
                       command=self.clear_caliper).pack(side=tk.RIGHT)

        def _build_verdict(self):
            row = ttk.Frame(self, padding=(10, 2))
            row.pack(side=tk.TOP, fill=tk.X)
            self.lbl_verdetto = ttk.Label(row, text="", anchor="w",
                                          font=("Menlo", 12))
            self.lbl_verdetto.pack(side=tk.LEFT)

        def _build_canvas(self):
            self.fig = Figure(figsize=(11, 4.2), dpi=100)
            self.ax = self.fig.add_subplot(111)
            self.fig.subplots_adjust(left=0.06, right=0.99, top=0.94, bottom=0.14)

            self.canvas = FigureCanvasTkAgg(self.fig, master=self)
            self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

            self.scroll = ttk.Scale(self, from_=0, to=1, orient=tk.HORIZONTAL,
                                    command=self.on_scroll)
            self.scroll.pack(side=tk.TOP, fill=tk.X, padx=8, pady=(0, 4))

            self.canvas.mpl_connect("motion_notify_event", self.on_motion)
            self.canvas.mpl_connect("button_press_event", self.on_click)
            self.cursor = None

        def _build_statusbar(self):
            st = ttk.Frame(self, padding=(8, 4))
            st.pack(side=tk.BOTTOM, fill=tk.X)

            self.lbl_pos = ttk.Label(st, text="-", width=34, anchor="w",
                                     font=("Menlo", 11))
            self.lbl_pos.pack(side=tk.LEFT)

            self.lbl_beat = ttk.Label(st, text="-", width=46, anchor="w",
                                      font=("Menlo", 11))
            self.lbl_beat.pack(side=tk.LEFT)

            self.lbl_cal = ttk.Label(st, text="caliper: clic per il primo punto",
                                     anchor="w", font=("Menlo", 11))
            self.lbl_cal.pack(side=tk.LEFT)

            self.lbl_file = ttk.Label(st, text="", anchor="e")
            self.lbl_file.pack(side=tk.RIGHT)

        # --------------------------------------------------------------- dati

        def load(self, path):
            paths = [path] if isinstance(path, str) else list(path)
            try:
                scp = merge_segments(paths)
            except SegmentiNonContigui as exc:
                messagebox.showerror("Segmenti non contigui", str(exc))
                return
            except Exception as exc:
                messagebox.showerror("Errore di lettura", str(exc))
                return
            if not scp.rhythm:
                messagebox.showerror("Errore",
                                     "il file non contiene dati di ritmo (sezione 6)")
                return

            self.scp = scp
            self.fs = scp.rhythm["fs"]
            names = [ld["name"] for ld in scp.leads] or \
                    [f"lead{i+1}" for i in range(len(scp.rhythm["signals"]))]
            self.cmb_lead["values"] = names
            self.lead_idx = 0
            self.var_lead.set(names[0])

            self.analisi = analyze_rhythm(scp)
            hr = scp.heart_rate()
            if getattr(scp, "segmenti", None):
                info = (f"{len(scp.segmenti)} segmenti uniti   "
                        f"{scp.duration_s():.0f} s   {self.fs:.1f} Hz")
            else:
                info = f"{os.path.basename(paths[0])}   {self.fs:.1f} Hz"
            if hr:
                info += f"   {hr['bpm_mean']:.0f} bpm medi su {hr['n_peaks']} battiti"
            if scp.rpeak_flag_used:
                info += "   [marcatori a bordo]"
            self.lbl_file.config(text=info)
            if self.analisi:
                q = self.analisi["qualita"]
                self.lbl_verdetto.config(
                    text=f"{self.analisi['verdetto']}    "
                         f"[rumore {q['rumore_uv']:.0f} uV, "
                         f"deriva {q['deriva_mv']:.2f} mV]   "
                         f"soglie in parte nostre - non e' una diagnosi")
            else:
                self.lbl_verdetto.config(
                    text="analisi del ritmo non disponibile: il file non "
                         "contiene marcatori di battito")

            self.clear_caliper()
            self.on_lead()

        def on_lead(self):
            names = list(self.cmb_lead["values"])
            self.lead_idx = names.index(self.var_lead.get()) if self.var_lead.get() in names else 0
            self.sig = self.scp.millivolts()[self.lead_idx]
            pk = self.scp.rhythm.get("rpeaks") or []
            self.peaks = pk[self.lead_idx] if self.lead_idx < len(pk) else np.array([], dtype=int)
            self.on_window()

        def duration(self):
            return len(self.sig) / self.fs if self.sig is not None else 0.0

        def window(self):
            return float(self.var_win.get())

        def on_window(self):
            span = max(0.0, self.duration() - self.window())
            self.scroll.config(to=span if span > 0 else 1)
            if float(self.scroll.get()) > span:
                self.scroll.set(0)
            self.redraw()

        def on_scroll(self, _value):
            self.redraw()

        # -------------------------------------------------------------- disegno

        def current_gain(self):
            g = self.var_gain.get()
            if g == "auto":
                return choose_gain([self.sig])
            return float(g)

        def redraw(self):
            if self.sig is None:
                return
            ax = self.ax
            ax.clear()
            self._caliper_art = []

            gain = self.current_gain()
            speed = float(self.var_speed.get())
            win = self.window()
            t0 = float(self.scroll.get())
            t1 = t0 + win

            half = 15.0 / gain
            peak = float(np.max(np.abs(self.sig)))
            if peak > half:
                half = peak * 1.1

            # griglia: 1 mm e 5 mm di carta tradotti nelle unita' degli assi
            ax.xaxis.set_minor_locator(MultipleLocator(1.0 / speed))
            ax.xaxis.set_major_locator(MultipleLocator(5.0 / speed))
            ax.xaxis.set_major_formatter(
                FuncFormatter(lambda x, p: f"{x:.1f}" if abs(x * 2 - round(x * 2)) < 1e-6 else ""))
            ax.yaxis.set_minor_locator(MultipleLocator(1.0 / gain))
            ax.yaxis.set_major_locator(MultipleLocator(5.0 / gain))
            ax.grid(which="minor", color="#f4b8b8", linewidth=0.4)
            ax.grid(which="major", color="#e06666", linewidth=0.8)

            i0 = max(0, int(t0 * self.fs) - 2)
            i1 = min(len(self.sig), int(t1 * self.fs) + 2)
            t = np.arange(i0, i1) / self.fs
            ax.plot(t, self.sig[i0:i1], color="#101010", linewidth=0.9)

            for g in getattr(self.scp, "giunzioni", []) or []:
                tg = g / self.fs
                if t0 <= tg <= t1:
                    ax.axvline(tg, color="#8e24aa", linewidth=1.0,
                               linestyle=(0, (4, 3)), zorder=2)

            if self.var_marks.get() and len(self.peaks):
                # ogni battito prende il colore dell'intervallo che lo precede:
                # blu normale, rosso anticipato, arancione intervallo lungo
                colori = {"normale": "#1f6fb2", "anticipato": "#d81b60",
                          "mancante": "#f57c00", "giunzione": "#8e8c84"}
                et = self.analisi["etichette"] if self.analisi else []
                for j, k in enumerate(self.peaks):
                    if not (i0 <= k < i1):
                        continue
                    lab = et[j - 1] if 0 < j <= len(et) else "normale"
                    ax.plot([k / self.fs], [-half * 0.88], marker="^",
                            color=colori.get(lab, "#1f6fb2"),
                            markersize=8 if lab != "normale" else 5, zorder=3)

            ax.set_xlim(t0, t1)
            ax.set_ylim(-half, half)
            if self.var_aspect.get():
                ax.set_aspect((1.0 / speed) / (1.0 / gain))
            else:
                ax.set_aspect("auto")
            ax.set_xlabel("secondi")
            ax.set_ylabel("mV")
            scala = "proporzioni reali" if self.var_aspect.get() else "adattato allo schermo"
            ax.set_title(f"{speed:g} mm/s   {gain:g} mm/mV   "
                         f"quadrato grande = {5/speed:g} s x {5/gain:g} mV   "
                         f"[{scala}]", fontsize=9)
            for side in ax.spines.values():
                side.set_visible(False)

            self._draw_caliper()

            # il mirino va ricreato a ogni ridisegno perche' usa il blitting
            self.cursor = Cursor(ax, useblit=True, color="#2f7d32",
                                 linewidth=0.8, alpha=0.8)
            self.canvas.draw()

        def _draw_caliper(self):
            for t, _v in self.caliper:
                self._caliper_art.append(
                    self.ax.axvline(t, color="#8e24aa", linewidth=1.2, linestyle="--"))
            if len(self.caliper) == 2:
                (ta, va), (tb, vb) = self.caliper
                y = self.ax.get_ylim()[1] * 0.82
                self.ax.annotate("", xy=(tb, y), xytext=(ta, y),
                                 arrowprops=dict(arrowstyle="<->", color="#8e24aa"))
                dt = abs(tb - ta)
                txt = f"{dt*1000:.0f} ms"
                if dt > 0:
                    txt += f"  =  {60.0/dt:.1f} bpm"
                self.ax.text((ta + tb) / 2, y * 1.06, txt, ha="center",
                             fontsize=9, color="#8e24aa")

        # -------------------------------------------------------------- eventi

        def on_motion(self, event):
            if event.inaxes is not self.ax or self.sig is None:
                self.lbl_pos.config(text="-")
                self.lbl_beat.config(text="-")
                return

            t = event.xdata
            idx = int(round(t * self.fs))
            if 0 <= idx < len(self.sig):
                self.lbl_pos.config(
                    text=f"t {t:7.3f} s   segnale {self.sig[idx]:+7.3f} mV")
            else:
                self.lbl_pos.config(text=f"t {t:7.3f} s")

            self.lbl_beat.config(text=self._beat_text(t))

        def _beat_text(self, t):
            if not len(self.peaks):
                return "nessun marcatore di battito nel file"
            pt = self.peaks / self.fs
            k = int(np.argmin(np.abs(pt - t)))
            parts = [f"battito {k+1}/{len(pt)} a {pt[k]:.3f} s"]
            if k > 0:
                rr = pt[k] - pt[k - 1]
                parts.append(f"RR {rr*1000:.0f} ms = {60/rr:.1f} bpm")
            return "   ".join(parts)

        def on_click(self, event):
            if event.inaxes is not self.ax or event.button != 1:
                return
            idx = int(round(event.xdata * self.fs))
            v = float(self.sig[idx]) if 0 <= idx < len(self.sig) else float("nan")

            if len(self.caliper) >= 2:
                self.caliper = []
            self.caliper.append((event.xdata, v))

            if len(self.caliper) == 1:
                self.lbl_cal.config(text="caliper: clic per il secondo punto")
            else:
                (ta, va), (tb, vb) = self.caliper
                dt = abs(tb - ta)
                dv = vb - va
                msg = f"caliper: {dt*1000:.0f} ms   {dv:+.3f} mV"
                if dt > 0:
                    msg += f"   equivale a {60/dt:.1f} bpm"
                self.lbl_cal.config(text=msg)
            self.redraw()

        def clear_caliper(self):
            self.caliper = []
            self.lbl_cal.config(text="caliper: clic per il primo punto")
            if self.sig is not None:
                self.redraw()

        # ------------------------------------------------------------- comandi

        def on_open(self):
            paths = filedialog.askopenfilenames(
                title="Apri uno o piu' segmenti SCP-ECG "
                      "(selezione multipla per unirli)",
                filetypes=[("SCP-ECG", "*.scp *.SCP"), ("Tutti i file", "*.*")])
            if paths:
                self.load(list(paths))

        def on_export_csv(self):
            """Scrive il CSV completo accanto al file SCP."""
            if self.scp is None:
                return
            stem = os.path.splitext(self.scp.path)[0]
            out = (stem + "_unito.csv") if getattr(self.scp, "segmenti", None) \
                else (stem + ".csv")
            try:
                write_csv(self.scp, out, "auto")
            except Exception as exc:
                messagebox.showerror("Errore", str(exc))
                return
            self.lbl_cal.config(text=f"CSV scritto: {os.path.basename(out)}")

        def on_export_page(self):
            """Scrive il PNG a pagina intera (tutti i 30 s su piu' righe)."""
            if self.scp is None:
                return
            stem = os.path.splitext(self.scp.path)[0]
            out = (stem + "_unito.png") if getattr(self.scp, "segmenti", None) \
                else (stem + ".png")
            gain = 0.0 if self.var_gain.get() == "auto" else float(self.var_gain.get())
            try:
                write_plot(self.scp, out, "auto", 10.0,
                           float(self.var_speed.get()), gain)
            except Exception as exc:
                messagebox.showerror("Errore", str(exc))
                return
            self.lbl_cal.config(text=f"PNG scritto: {os.path.basename(out)}")

        def on_save(self):
            if self.sig is None:
                return
            path = filedialog.asksaveasfilename(
                defaultextension=".png", filetypes=[("PNG", "*.png")],
                title="Salva la vista corrente")
            if path:
                self.fig.savefig(path, dpi=200, facecolor="white")
                self.lbl_cal.config(text=f"salvato: {os.path.basename(path)}")

    ScpViewer(path).mainloop()


# ==========================================================================

def main():
    ap = argparse.ArgumentParser(
        description="lettore, convertitore e visualizzatore SCP-ECG; piu' file "
                    "insieme vengono uniti come segmenti consecutivi")
    ap.add_argument("file", nargs="*", help="uno o piu' file .scp")
    ap.add_argument("--batch", action="store_true",
                    help="elabora ogni file separatamente, senza unirli e "
                         "senza aprire la finestra")
    ap.add_argument("--info", action="store_true",
                    help="stampa solo struttura e analisi, non produce file")
    ap.add_argument("--csv", metavar="PATH", help="percorso del CSV")
    ap.add_argument("--plot", metavar="PATH", help="percorso del PNG")
    ap.add_argument("--invert", action="store_true",
                    help="ribalta il segnale (verifica di elettrodi invertiti)")
    ap.add_argument("--baseline", default="auto",
                    help="'auto' (default), 'none' oppure un valore numerico")
    ap.add_argument("--seconds-per-row", type=float, default=10.0,
                    help="secondi per riga nel grafico (default 10)")
    ap.add_argument("--mm-per-s", type=float, default=25.0,
                    help="velocita' della carta in mm/s (default 25)")
    ap.add_argument("--mm-per-mv", type=float, default=0.0,
                    help="guadagno in mm/mV (default 0 = automatico)")
    ap.add_argument("--version", action="version",
                    version=f"scpecg {__version__} - {__author__} - "
                            f"licenza {__license__}\n{__url__}")
    args = ap.parse_args()

    if not (args.batch or args.info or args.csv or args.plot):
        try:
            run_gui(args.file or None)
        except ImportError as exc:
            print(f"interfaccia non disponibile ({exc}).", file=sys.stderr)
            print("Su macOS con Python di Homebrew serve: "
                  "brew install python-tk@<versione>", file=sys.stderr)
            print("In alternativa usa --info o --batch.", file=sys.stderr)
            return 1
        return 0

    if not args.file:
        ap.error("serve almeno un file")

    baseline = args.baseline
    if baseline not in ("auto", "none"):
        baseline = float(baseline)

    def elabora(scp, stem):
        scp.invert = args.invert
        print_info(scp)
        print_segments(scp)
        print_rate(scp)
        print_analysis(scp)
        if args.invert:
            print("\nATTENZIONE: segnale ribaltato via software (--invert). "
                  "Utile per capire se gli elettrodi erano invertiti,\n"
                  "            ma il tracciato da conservare e' quello "
                  "registrato con gli elettrodi nella posizione corretta.")
        if args.info:
            return
        write_csv(scp, args.csv or stem + ".csv", baseline)
        write_plot(scp, args.plot or stem + ".png", baseline,
                   args.seconds_per_row, args.mm_per_s, args.mm_per_mv)

    if args.batch:
        if len(args.file) > 1 and (args.csv or args.plot):
            ap.error("--csv e --plot valgono per un file solo")
        for path in args.file:
            if len(args.file) > 1:
                print(f"\n===== {path} =====")
            try:
                elabora(ScpFile(path), os.path.splitext(path)[0])
            except Exception as exc:
                print(f"{path}: {exc}", file=sys.stderr)
        return 0

    try:
        scp = merge_segments(args.file)
    except SegmentiNonContigui as exc:
        print(f"\nSegmenti non contigui:\n{exc}", file=sys.stderr)
        return 2
    stem = os.path.splitext(args.file[0])[0]
    if len(args.file) > 1:
        stem += "_unito"
    elabora(scp, stem)
    return 0


if __name__ == "__main__":
    sys.exit(main())
