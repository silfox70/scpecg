# Il formato dei file .SCP del Creative / Gima PC-80B

Note di reverse engineering, ricavate leggendo i byte di tracciati reali e
verificate ricostruendo file sintetici byte per byte.

Lo standard di riferimento è SCP-ECG, definito da ANSI/AAMI EC71:2001 e
CEN EN 1064:2005, poi ripreso come ISO 11073-91064.

## Struttura del record

Preambolo di 6 byte, poi le sezioni.

| offset | campo |
|---|---|
| 0–1 | CRC-CCITT dell'intero record (polinomio 0x1021, init 0xFFFF), calcolato dal byte 2 in poi |
| 2–5 | lunghezza totale del record, little-endian |
| 6… | sezione 0 |

Ogni sezione ha un header di 16 byte: CRC (2), ID (2), lunghezza (4),
versione (1), protocollo (1), riservato (6). I dati cominciano al byte 16.

**Gli offset nei puntatori della sezione 0 sono 1-based**: vanno decrementati
di uno per ottenere l'offset di file. Verifica rapida: la sezione 0 dichiara
offset 7 e comincia al byte 6.

## Sezioni presenti in un file tipico (9796 byte)

| id | lunghezza | offset dichiarato | contenuto |
|---|---|---|---|
| 0 | 76 | 7 | puntatori (6 voci da 10 byte) |
| 1 | 32 | 83 | solo data e ora, poi terminatore |
| 2 | 30 | 115 | tabelle di Huffman |
| 3 | 28 | 145 | definizione delle derivazioni |
| 6 | 9024 | 173 | campioni del ritmo |
| 9 | 600 | 9197 | dati proprietari del costruttore |

Assenti le sezioni 4, 5, 7 e 8: nessun battito di riferimento, nessun
complesso mediato, nessuna misura o diagnosi automatica scritta nel file.
L'apparecchio, quando analizza, tiene il verdetto sul display.

## Sezione 1 — anagrafica

Contiene **solo** i tag 25 (data) e 26 (ora), poi il terminatore 255. Nessun
nome, nessun identificativo del paziente. Buono per la privacy.

Attenzione: se l'orologio interno non è stato impostato, la data di
acquisizione esce come `2000-01-04` o simile, uguale su tutti i file. Conviene
impostare l'RTC o rinominare i file allo scarico.

## Sezione 2 — compressione

I primi due byte sono il numero di tabelle di Huffman, i due successivi il
numero di strutture di codice.

- `19999` (0x4E1F) significa "uso la tabella di default dello standard"
- **`1` tabella con `1` struttura di codice significa nessuna compressione**

Il PC-80B scrive il secondo caso: `01 00 01 00`. I campioni sono quindi interi
grezzi, e non serve implementare un decoder Huffman.

Verifica dei conti: 2 + 2 + 9 (una struttura) = 13, più 1 di padding perché
SCP-ECG vuole lunghezze pari = 14 byte di dati, più 16 di header = 30. Che è
esattamente la lunghezza dichiarata.

## Sezione 3 — derivazioni

Un byte con il numero di derivazioni, un byte di flag, poi 9 byte per ogni
derivazione: campione iniziale (4), campione finale (4), codice (1).

Il PC-80B dichiara **una** derivazione con **codice 101**, che non esiste nella
tabella dello standard (che arriva a 64, aVF). È un valore proprietario e non
distingue la modalità di misura.

Secondo il manuale del costruttore:

| modalità | derivazione equivalente |
|---|---|
| palmo contro palmo | I |
| cavetti con elettrodi | II (predefinita) |
| appoggio al torace | simile a una derivazione V |

Il file non registra quale hai usato: va annotato a parte, altrimenti la
morfologia non è interpretabile.

## Sezione 6 — i campioni

Struttura dei dati, dopo i 16 byte di header:

| offset nei dati | campo | valore tipico |
|---|---|---|
| 0–1 | AVM, nanovolt per unità | 806 |
| 2–3 | intervallo di campionamento, µs | 6666 → 150,02 Hz |
| 4 | codifica differenze (0 assoluti, 1 prima, 2 seconda) | 0 |
| 5 | flag compressione bimodale | 0 |
| 6–7 | lunghezza in byte del blocco della derivazione 1 | 9000 |
| 8… | campioni, 16 bit little-endian | |

9000 byte = 4500 campioni = 29,997 s.

Coerenza: 16 + 8 + 9000 = 9024, la lunghezza dichiarata dal puntatore.

### Conversione in millivolt

I campioni **non** sono centrati sullo zero come vorrebbe lo standard: sono i
codici grezzi di un ADC a 12 bit con offset di mezza scala.

```
mV = (campione - 2048) * 806e-6
```

Controprova sull'elettronica: 806 nV × 4096 livelli = 3,30 mV di fondo scala,
cioè un ADC a 12 bit con riferimento 3,3 V e amplificatore a guadagno 1000.

In pratica conviene sottrarre la **mediana** invece del 2048 fisso: la mediana
di un ECG cade sulla linea isoelettrica (i QRS sono picchi brevi e non la
spostano) e così si elimina gratis anche la deriva della linea di base.

### Il bit 15 è un marcatore, non un segno

**Questa è la trappola principale.** In alcune registrazioni, circa un campione
per battito ha il bit 15 acceso. Letti come `int16` con segno, quei campioni
diventano valori attorno a −30000, il grafico sembra ribaltato e il picco
calcolato arriva a decine di millivolt — impossibile su un ADC che ha 3,3 mV
di fondo scala.

Sono i **picchi R rilevati a bordo dall'apparecchio**. Vanno trattati così:

```python
u = np.frombuffer(chunk, dtype="<u2").astype(np.int64)
flag = (u & 0x8000) != 0        # posizioni dei picchi R
val  = (u & 0x7FFF)             # campioni veri
```

Come riconoscere il caso senza sapere a priori che apparecchio ha scritto il
file, e senza rovinare un file SCP scritto correttamente con veri interi con
segno? Due condizioni insieme:

1. i campioni col bit 15 acceso sono meno del 2% del totale — un ECG con
   valori davvero negativi ne avrebbe migliaia
2. una volta mascherati, i loro valori cadono nello stesso intervallo degli
   altri campioni

Verifica di conferma: le distanze fra marcatori consecutivi devono
corrispondere agli intervalli RR. Su un file a 150 Hz e 100 bpm sono uscite
90–93 campioni, media 90,4 → 0,602 s → 99,7 bpm. Coerente.

I marcatori compaiono solo in alcune registrazioni, e sembrano legati alla
modalità di acquisizione: presenti nelle misure con i cavetti, assenti in
quelle a palmo, dove l'apparecchio esegue invece la propria analisi e mostra
il verdetto sul display.

## Sezione 9 — proprietaria

600 byte non documentati. Non ancora decodificati.

## Caratteristiche dichiarate dal costruttore

| parametro | valore |
|---|---|
| banda passante | 0,5 – 40 Hz |
| rumore interno | ≤ 30 µV picco-picco |
| range frequenza cardiaca | 30 – 240 bpm |
| accuratezza frequenza | ±2 bpm oppure 2% |

Nota: con la banda tagliata a 40 Hz, campionare a 150 Hz è abbondante
(Nyquist chiederebbe 80 Hz). Il QRS appare arrotondato per via del **filtro**,
non del campionamento.

## Verdetti dell'analisi automatica

Il manuale elenca 17 messaggi, combinazioni di tre assi indipendenti:

- **frequenza**: normale, lievemente rapida, rapida, lievemente lenta, lenta
- **regolarità**: regolare, intervallo occasionalmente corto, intervalli
  irregolari, serie di battiti ravvicinati
- **qualità**: pulito, deriva della linea di base, segnale scarso

Nessuna delle categorie richiede analisi morfologica: bastano gli intervalli
RR e un indice di rumore. Le soglie numeriche **non sono pubblicate** in
nessuna versione del manuale.

Definizioni che il manuale fornisce:

- ritmo normale 60–100 bpm
- battito mancante: intervallo pari al doppio della media dei precedenti
- bigeminismo: un battito normale accoppiato a uno prematuro
- trigeminismo: due normali accoppiati a uno prematuro
- serie: extrasistoli più di tre volte consecutive
