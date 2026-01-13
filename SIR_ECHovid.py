import random
import math
from numbers import Number

infettivi = None
gamma = None
sigma = None
tau = None
distanziamento = None
guariti = None
deceduti = None
infettivi_storico = None
deceduti_storico = None
guariti_storico = None

def upRange(start, stop, step):
  while start <= stop:
    yield start
    start += abs(step)

def downRange(start, stop, step):
  while start >= stop:
    yield start
    start -= abs(step)

"""Crea una matrice di c colonne e r righe.
La matrice è implementata come una lista di liste in cui
ogni elemento è una riga, ogni riga è una lista di elementi.
Parametri:
c, r: il numero di colonne e righe della matrice (due numeri naturali)
Valore restituito:
    la matrice (una lista)
"""
def CreaMatrice(c, r):
  global param_sigma, param_lambda, param_gamma, param_tau, param_distanziamento, param_interazioni, m, infettivi, durata, output_parziale, mu, sigma, valore, riga, colonna, x1, y1, x2, y2, interazioni, giorno, prima_riga, c_v_a, n_c, tizio, d, i, stringa_matrice, lambda2, segno, argomento_radice, n_r, riga_interessata, guariti, deceduti, gamma, infettivi_storico, stato_individuo, nuova_riga, tau, risultati_giorno, coordinate, j, distanziamento, deceduti_storico, caio, guariti_storico
  m = []
  if r > 1 and c > 1:
    for i in (1 <= float(r)) and upRange(1, float(r), 1) or downRange(1, float(r), 1):
      nuova_riga = []
      for j in (1 <= float(c)) and upRange(1, float(c), 1) or downRange(1, float(c), 1):
        nuova_riga.append(0)
      m.append(nuova_riga)
  return m

"""Inizializza i parametri e le variabili per la simulazione. Inizializza
le liste con lo storico dei valori di infettivi, deceduti e guariti.
Parametri:
param_sigma: varianza, usato per determinare la distanza
dell'individuo con cui interagire (numero reale)
param_lambda: tasso di contagio, probabilità che un
infettivo contagi un suscettibile (numero reale)
param_gamma: tasso di guarigione, probabilità
che un infettivo guarisca (numero reale)
param_tau: tasso di mortalità, probabilità che un infettivo muoia (numero reale)
param_distanziamento: valore che incide sul distanziamento sociale (numero reale)
param_interazioni: numero di interazioni per ciascun individuo al giorno
Valore resituito:
    nessuno
"""
def Inizializzazione(param_sigma, param_lambda, param_gamma, param_tau, param_distanziamento, param_interazioni):
  global c, r, m, infettivi, durata, output_parziale, mu, sigma, valore, riga, colonna, x1, y1, x2, y2, interazioni, giorno, prima_riga, c_v_a, n_c, tizio, d, i, stringa_matrice, lambda2, segno, argomento_radice, n_r, riga_interessata, guariti, deceduti, gamma, infettivi_storico, stato_individuo, nuova_riga, tau, risultati_giorno, coordinate, j, distanziamento, deceduti_storico, caio, guariti_storico
  sigma = param_sigma
  lambda2 = param_lambda
  gamma = param_gamma
  tau = param_tau
  distanziamento = param_distanziamento
  interazioni = param_interazioni
  deceduti_storico = []
  deceduti_storico.append(0)
  infettivi_storico = []
  infettivi_storico.append(0)
  guariti_storico = []
  guariti_storico.append(0)

"""Restituisce il numero di righe di una matrice (la
matrice è implementata come una lista di liste).
Parametri:
    m: la matrice (una lista)
Valore restituito:
il numero di righe (un numero naturale, la lunghezza della lista)
"""
def RigheMatrice(m):
  global param_sigma, param_lambda, param_gamma, param_tau, param_distanziamento, param_interazioni, c, r, infettivi, durata, output_parziale, mu, sigma, valore, riga, colonna, x1, y1, x2, y2, interazioni, giorno, prima_riga, c_v_a, n_c, tizio, d, i, stringa_matrice, lambda2, segno, argomento_radice, n_r, riga_interessata, guariti, deceduti, gamma, infettivi_storico, stato_individuo, nuova_riga, tau, risultati_giorno, coordinate, j, distanziamento, deceduti_storico, caio, guariti_storico
  return len(m)

"""Imposta una situazione iniziale inserendo alcuni
infettivi in posizioni scelte casualmente.
Parametri:
    m: la matrice (una lista)
    infettivi: il numero di infettivi da inserire
Valore restituito:
la matrice modificata dalla presenza degli infettivi (una lista)
"""
def SituazioneInizialeCasuale(m, infettivi):
  global param_sigma, param_lambda, param_gamma, param_tau, param_distanziamento, param_interazioni, c, r, durata, output_parziale, mu, sigma, valore, riga, colonna, x1, y1, x2, y2, interazioni, giorno, prima_riga, c_v_a, n_c, tizio, d, i, stringa_matrice, lambda2, segno, argomento_radice, n_r, riga_interessata, guariti, deceduti, gamma, infettivi_storico, stato_individuo, nuova_riga, tau, risultati_giorno, coordinate, j, distanziamento, deceduti_storico, caio, guariti_storico
  for count in range(int(infettivi)):
    c = random.randint(1, ColonneMatrice(m))
    r = random.randint(1, ColonneMatrice(m))
    while LeggiElementoMatrice(m, c, r) > 0:
      c = random.randint(1, ColonneMatrice(m))
      r = random.randint(1, ColonneMatrice(m))
    ScriviElementoMatrice(m, c, r, 1)
  infettivi_storico[0] = infettivi
  return m

"""Procedura principale dell'applicazione: avvia la simulazione.
Chiama PropagazioneContagio, FineGiornata e Risultati.
Parametri:
    m: matrice della popolazione
durata: unità di tempo (giorni) della simulazione (un numero naturale)
output_parziale: se True mostra lo stato della simulazione per ogni unità di tempo
Valore resituito:
    la matrice (una lista) con la situazione finale
"""
def SimulazioneSIR(m, durata, output_parziale):
  global param_sigma, param_lambda, param_gamma, param_tau, param_distanziamento, param_interazioni, c, r, infettivi, mu, sigma, valore, riga, colonna, x1, y1, x2, y2, interazioni, giorno, prima_riga, c_v_a, n_c, tizio, d, i, stringa_matrice, lambda2, segno, argomento_radice, n_r, riga_interessata, guariti, deceduti, gamma, infettivi_storico, stato_individuo, nuova_riga, tau, risultati_giorno, coordinate, j, distanziamento, deceduti_storico, caio, guariti_storico
  for giorno in (1 <= float(durata)) and upRange(1, float(durata), 1) or downRange(1, float(durata), 1):
    if not PropagazioneContagio(m, interazioni):
      return False
    FineGiornata(m)
    risultati_giorno = Risultati(m)
    if output_parziale:
      print(risultati_giorno)
  return m

"""Restituisce il numero di colonne di una matrice (implementata come
lista di liste). Ricava un elemento della matrice (una lista)
e ne determina il numero di elementi (il numero di colonne).
Parametri:
    m: la matrice (una lista)
Valore restituito:
il numero di colonne (un numero naturale, la lunghezza di una riga)
"""
def ColonneMatrice(m):
  global param_sigma, param_lambda, param_gamma, param_tau, param_distanziamento, param_interazioni, c, r, infettivi, durata, output_parziale, mu, sigma, valore, riga, colonna, x1, y1, x2, y2, interazioni, giorno, prima_riga, c_v_a, n_c, tizio, d, i, stringa_matrice, lambda2, segno, argomento_radice, n_r, riga_interessata, guariti, deceduti, gamma, infettivi_storico, stato_individuo, nuova_riga, tau, risultati_giorno, coordinate, j, distanziamento, deceduti_storico, caio, guariti_storico
  prima_riga = m[0]
  return len(prima_riga)

"""Determina un valore casuale usando una distribuzione casuale
il cui valore atteso è mu e la varianza sigma. In pratica,
si sceglie un ordinata casuale tra zero e il massimo della
gaussiana e poi si determina un'ascissa compresa tra la media e
l'ascissa corrispondente all'ordinata determinata casualmente.
Parametri:
    mu: valore atteso (un numero reale)
sigma: varianza, tanto più è piccolo tanto più il valore
casuale è "vicino" al valore atteso (un numero reale)
"""
def CasualeValoreAtteso(mu, sigma):
  global param_sigma, param_lambda, param_gamma, param_tau, param_distanziamento, param_interazioni, c, r, m, infettivi, durata, output_parziale, valore, riga, colonna, x1, y1, x2, y2, interazioni, giorno, prima_riga, c_v_a, n_c, tizio, d, i, stringa_matrice, lambda2, segno, argomento_radice, n_r, riga_interessata, guariti, deceduti, gamma, infettivi_storico, stato_individuo, nuova_riga, tau, risultati_giorno, coordinate, j, distanziamento, deceduti_storico, caio, guariti_storico
  if random.random() < 0.5:
    segno = 1
  else:
    segno = -1
  argomento_radice = math.sqrt(-2 * math.log((random.random() * (1 / (sigma * math.sqrt(2 * math.pi)))) * (sigma * math.sqrt(2 * math.pi))))
  c_v_a = math.ceil(mu + segno * (sigma * argomento_radice))
  return c_v_a

"""Scrive un valore in una specifica posizione in una
matrice (la matrice è implementata come lista di liste).
Parametri:
    m: la matrice (una lista)
c, r: la posizione dell'elemento (due numeri
naturali), colonna e riga nella matrice
valore: il valore da inserire nell'elemento della matrice (di tipo qualsiasi)
Valore restituito:
    nessuno
"""
def ScriviElementoMatrice(m, c, r, valore):
  global param_sigma, param_lambda, param_gamma, param_tau, param_distanziamento, param_interazioni, infettivi, durata, output_parziale, mu, sigma, riga, colonna, x1, y1, x2, y2, interazioni, giorno, prima_riga, c_v_a, n_c, tizio, d, i, stringa_matrice, lambda2, segno, argomento_radice, n_r, riga_interessata, guariti, deceduti, gamma, infettivi_storico, stato_individuo, nuova_riga, tau, risultati_giorno, coordinate, j, distanziamento, deceduti_storico, caio, guariti_storico
  if c >= 1 and c <= ColonneMatrice(m):
    if r >= 1 and r <= RigheMatrice(m):
      riga_interessata = m[int(r - 1)]
      riga_interessata[int(c - 1)] = valore

"""Sceglie un individuo tra la popolazione diverso da quello in
esame (è usata dalla funzione Contagio per simulare gli incontri,
le interazioni tra individui). La scelta è fatta calcolando una
riga e una colonna casuali mediante la funzione CasualeValoreAtteso
Parametri:
    m: la matrice della popolazione (una lista)
c, r: posizione nella matrice dell'individuo in esame (due numeri naturali)
Valore restituito:
una lista con colonna e riga dell'elemento scelto (diverso da quello in esame)
"""
def ScegliIndividuo(m, c, r):
  global param_sigma, param_lambda, param_gamma, param_tau, param_distanziamento, param_interazioni, infettivi, durata, output_parziale, mu, sigma, valore, riga, colonna, x1, y1, x2, y2, interazioni, giorno, prima_riga, c_v_a, n_c, tizio, d, i, stringa_matrice, lambda2, segno, argomento_radice, n_r, riga_interessata, guariti, deceduti, gamma, infettivi_storico, stato_individuo, nuova_riga, tau, risultati_giorno, coordinate, j, distanziamento, deceduti_storico, caio, guariti_storico
  n_c = c
  n_r = r
  # Continua a scegliere un individuo (una colonna e una riga
  # della matrice) finché l'individuo scelto è diverso da quello
  # in esame (almeno una tra colonna e riga sono differenti)
  while n_c == c and n_r == r:
    n_c = CasualeValoreAtteso(c, sigma)
    n_r = CasualeValoreAtteso(r, sigma)
    if n_c < 1 or n_c > ColonneMatrice(m):
      n_c = c
      n_r = r
    elif n_r < 1 or n_r > RigheMatrice(m):
      n_c = c
      n_r = r
  return [n_c, n_r]

"""Ricava il valore in una specifica posizione in una
matrice (la matrice è implementata come lista di liste).
Parametri:
    m: la matrice (una lista)
c, r: la posizione dell'elemento (due numeri
naturali), colonna e riga nella matrice
Valore restituito:
il contenuto dell'elemento della matrice (di tipo qualsiasi)
"""
def LeggiElementoMatrice(m, c, r):
  global param_sigma, param_lambda, param_gamma, param_tau, param_distanziamento, param_interazioni, infettivi, durata, output_parziale, mu, sigma, valore, riga, colonna, x1, y1, x2, y2, interazioni, giorno, prima_riga, c_v_a, n_c, tizio, d, i, stringa_matrice, lambda2, segno, argomento_radice, n_r, riga_interessata, guariti, deceduti, gamma, infettivi_storico, stato_individuo, nuova_riga, tau, risultati_giorno, coordinate, j, distanziamento, deceduti_storico, caio, guariti_storico
  if c < 1 or c > ColonneMatrice(m):
    return float('inf')
  if r < 1 or r > RigheMatrice(m):
    return float('inf')
  riga_interessata = m[int(r - 1)]
  return riga_interessata[int(c - 1)]

"""Simula le interazioni, gli incontri, di uno specifico individuo con
un altro scelto casualmente chiamando la funzione ScegliIndividuo.
Parametri:
    m: la matrice della popolazione (una lista)
riga, colonna: la posizione dell'individuo nella matrice (due numeri naturali)
Valore restituito:
    nessuno
"""
def Contagio(m, riga, colonna):
  global param_sigma, param_lambda, param_gamma, param_tau, param_distanziamento, param_interazioni, c, r, infettivi, durata, output_parziale, mu, sigma, valore, x1, y1, x2, y2, interazioni, giorno, prima_riga, c_v_a, n_c, tizio, d, i, stringa_matrice, lambda2, segno, argomento_radice, n_r, riga_interessata, guariti, deceduti, gamma, infettivi_storico, stato_individuo, nuova_riga, tau, risultati_giorno, coordinate, j, distanziamento, deceduti_storico, caio, guariti_storico
  tizio = LeggiElementoMatrice(m, colonna, riga)
  x1 = colonna
  y1 = riga
  # Sceglie un individuo diverso da quello in esame, la variabile
  # coordinate contiene la colonna e la riga dell'individuo scelto.
  coordinate = ScegliIndividuo(m, colonna, riga)
  x2 = coordinate[0]
  y2 = coordinate[-1]
  caio = LeggiElementoMatrice(m, x2, y2)
  # Controlla se avviene il contagio, in altre parole controlla se solo
  # uno dei due individui è infettivo e solo uno dei due è suscettibile.
  if tizio == 0 and caio > 0 or caio == 0 and tizio > 0:
    # Il contagio avviene con una probabilità pari a lambda pesata dal distanziamento.
    if random.random() < lambda2 / distanziamento:
      if tizio == 0:
        ScriviElementoMatrice(m, x1, y1, 1)
      else:
        ScriviElementoMatrice(m, x2, y2, 1)

"""Calcola la distanza tra due elementi della matrice come fossero due
punti su un piano cartesiano e riga e colonna ordinata e ascissa.
Parametri:
    m: la matrice (una lista)
x1, y1: la posizione del primo elemento (due numeri
naturali), colonna e riga nella matrice di questo elemento
x2, y2: la posizione del secondo elemento (due numeri
naturali), colonna e riga nella matrice di questo elemento
Valore restituito:
la distanza cartesiana tra i due elementi (un numero reale)
"""
def DistanzaMatrice(m, x1, y1, x2, y2):
  global param_sigma, param_lambda, param_gamma, param_tau, param_distanziamento, param_interazioni, c, r, infettivi, durata, output_parziale, mu, sigma, valore, riga, colonna, interazioni, giorno, prima_riga, c_v_a, n_c, tizio, d, i, stringa_matrice, lambda2, segno, argomento_radice, n_r, riga_interessata, guariti, deceduti, gamma, infettivi_storico, stato_individuo, nuova_riga, tau, risultati_giorno, coordinate, j, distanziamento, deceduti_storico, caio, guariti_storico
  d = x1 ** 2 + y1 ** 2
  d = (d if isinstance(d, Number) else 0) + (x2 ** 2 + y2 ** 2)
  return math.sqrt(d)

"""Simula l'interazione di ogni individuo della popolazione controlla
lo sviluppo del contagio chiamando la funzione Contagio.
Parametri:
    m: la matrice (una lista)
interazioni: il numero di interazioni, di incontri, di
ciascun individuo della popolazione (un numero naturale)
Valore resituito:
    True se non ci sono errori
"""
def PropagazioneContagio(m, interazioni):
  global param_sigma, param_lambda, param_gamma, param_tau, param_distanziamento, param_interazioni, c, r, infettivi, durata, output_parziale, mu, sigma, valore, riga, colonna, x1, y1, x2, y2, giorno, prima_riga, c_v_a, n_c, tizio, d, i, stringa_matrice, lambda2, segno, argomento_radice, n_r, riga_interessata, guariti, deceduti, gamma, infettivi_storico, stato_individuo, nuova_riga, tau, risultati_giorno, coordinate, j, distanziamento, deceduti_storico, caio, guariti_storico
  riga_end = float(RigheMatrice(m))
  for riga in (1 <= riga_end) and upRange(1, riga_end, 1) or downRange(1, riga_end, 1):
    colonna_end = float(ColonneMatrice(m))
    for colonna in (1 <= colonna_end) and upRange(1, colonna_end, 1) or downRange(1, colonna_end, 1):
      for count2 in range(int(interazioni)):
        Contagio(m, riga, colonna)
  return True

"""Mostra il contenuto della matrice. Funzione di debug.
Parametri:
    m: la matrice (una lista)
Valore restituito:
    nessuno
"""
def StampaMatrice(m):
  global param_sigma, param_lambda, param_gamma, param_tau, param_distanziamento, param_interazioni, c, r, infettivi, durata, output_parziale, mu, sigma, valore, riga, colonna, x1, y1, x2, y2, interazioni, giorno, prima_riga, c_v_a, n_c, tizio, d, i, stringa_matrice, lambda2, segno, argomento_radice, n_r, riga_interessata, guariti, deceduti, gamma, infettivi_storico, stato_individuo, nuova_riga, tau, risultati_giorno, coordinate, j, distanziamento, deceduti_storico, caio, guariti_storico
  i_end = float(RigheMatrice(m))
  for i in (1 <= i_end) and upRange(1, i_end, 1) or downRange(1, i_end, 1):
    print(m[int(i - 1)])

"""Determina lo stato della simulazione alla fine di ogni unità di
tempo (giorno). In sostanza: incrementa di uno lo stato degli
infettivi (per gli infettivi lo stato indica il numero di giorni
dal contagio); controlla se un infettivo guarisce o muore.
Parametri:
    m: la matrice (una lista di liste)
Valore restituito:
    nessuno
"""
def FineGiornata(m):
  global param_sigma, param_lambda, param_gamma, param_tau, param_distanziamento, param_interazioni, c, r, infettivi, durata, output_parziale, mu, sigma, valore, riga, colonna, x1, y1, x2, y2, interazioni, giorno, prima_riga, c_v_a, n_c, tizio, d, i, stringa_matrice, lambda2, segno, argomento_radice, n_r, riga_interessata, guariti, deceduti, gamma, infettivi_storico, stato_individuo, nuova_riga, tau, risultati_giorno, coordinate, j, distanziamento, deceduti_storico, caio, guariti_storico
  riga_end2 = float(RigheMatrice(m))
  for riga in (1 <= riga_end2) and upRange(1, riga_end2, 1) or downRange(1, riga_end2, 1):
    colonna_end2 = float(ColonneMatrice(m))
    for colonna in (1 <= colonna_end2) and upRange(1, colonna_end2, 1) or downRange(1, colonna_end2, 1):
      stato_individuo = LeggiElementoMatrice(m, colonna, riga)
      # Se l'individuo è un infettivo (stato maggiore di
      # zero), lo incrementa: un giorno in più dal contagio
      if stato_individuo > 0:
        stato_individuo = (stato_individuo if isinstance(stato_individuo, Number) else 0) + 1
        # Gamma indica il tasso di guarigione, la
        # probabilità che un infettivo sia rimosso (guarito)
        # Tau indica il tasso di mortaliztà, la probabilità
        # che un infettivo sia rimosso (deceduto)
        if random.random() < tau:
          stato_individuo = -1
        elif random.random() < gamma:
          stato_individuo = -2
        ScriviElementoMatrice(m, colonna, riga, stato_individuo)

"""Trasforma in stringa il contenuto di una matrice
(la matrice è implementata come lista di liste).
Ogni elemento è separato da una virgola, ogni riga è racchiusa
tra [ e ] e separata da una virgola dalla riga seguente.
Parametri:
    m: matrice (una lista di liste)
Valore restiuito:
    una string di testo con gli elementi della matrice
"""
def TestoMatrice(m):
  global param_sigma, param_lambda, param_gamma, param_tau, param_distanziamento, param_interazioni, c, r, infettivi, durata, output_parziale, mu, sigma, valore, riga, colonna, x1, y1, x2, y2, interazioni, giorno, prima_riga, c_v_a, n_c, tizio, d, i, stringa_matrice, lambda2, segno, argomento_radice, n_r, riga_interessata, guariti, deceduti, gamma, infettivi_storico, stato_individuo, nuova_riga, tau, risultati_giorno, coordinate, j, distanziamento, deceduti_storico, caio, guariti_storico
  stringa_matrice = '['
  i_end2 = float(len(m))
  for i in (1 <= i_end2) and upRange(1, i_end2, 1) or downRange(1, i_end2, 1):
    stringa_matrice = str(stringa_matrice) + str(''.join([str(x) for x in ['[', m[int(i - 1)], ']']]))
    if i < len(m):
      stringa_matrice = str(stringa_matrice) + ','
  stringa_matrice = str(stringa_matrice) + ']'
  return stringa_matrice

"""Determina il numero attuale di infettivi, guariti e deceduti.
I valori sono memorizzati in tre liste (variabili globali).
Parametri:
    m: la matrice (una lista)
Valore restituito:
una lista di tre elementi: infettivi, deceduti e guariti attuali
"""
def Risultati(m):
  global param_sigma, param_lambda, param_gamma, param_tau, param_distanziamento, param_interazioni, c, r, infettivi, durata, output_parziale, mu, sigma, valore, riga, colonna, x1, y1, x2, y2, interazioni, giorno, prima_riga, c_v_a, n_c, tizio, d, i, stringa_matrice, lambda2, segno, argomento_radice, n_r, riga_interessata, guariti, deceduti, gamma, infettivi_storico, stato_individuo, nuova_riga, tau, risultati_giorno, coordinate, j, distanziamento, deceduti_storico, caio, guariti_storico
  infettivi = 0
  guariti = 0
  deceduti = 0
  riga_end3 = float(RigheMatrice(m))
  for riga in (1 <= riga_end3) and upRange(1, riga_end3, 1) or downRange(1, riga_end3, 1):
    colonna_end3 = float(ColonneMatrice(m))
    for colonna in (1 <= colonna_end3) and upRange(1, colonna_end3, 1) or downRange(1, colonna_end3, 1):
      stato_individuo = LeggiElementoMatrice(m, colonna, riga)
      if stato_individuo > 0:
        infettivi = (infettivi if isinstance(infettivi, Number) else 0) + 1
      elif stato_individuo == -1:
        deceduti = (deceduti if isinstance(deceduti, Number) else 0) + 1
      elif stato_individuo <= -2:
        guariti = (guariti if isinstance(guariti, Number) else 0) + 1
  infettivi_storico.append(infettivi)
  deceduti_storico.append(deceduti)
  guariti_storico.append(guariti)
  return [infettivi, guariti, deceduti]


Inizializzazione(4, 0.06, 0.05, 0.05, 1, 4)
m = CreaMatrice(5, 5)
m = SituazioneInizialeCasuale(m, 2)
m = SimulazioneSIR(m, 10, False)
print(infettivi_storico)
print(TestoMatrice(m))
