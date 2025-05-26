import analisi as lib
import numpy as np
from iminuit import Minuit
from iminuit.cost import LeastSquares
import matplotlib.pyplot as plt
from scipy.signal import find_peaks



# sferica -> 1/r o 1/r^2 (se abbiamo intensità o campo)
# piana: no dipend. da r
unosur = lambda r,a,b: a/r + b # andamento 1/r
# angolo 30
r = np.array([ 17.5, 17.8, 18, 18.2, 18.8, 19.1, 19.4,  19.7, 19.9, 20.1, 20.3, 20.5, 20.7, 20.9, 21.1,
              21.3, 21.5, 21.7, 21.9, 22.1,  22.3, 22.5, 22.7, 22.9, 23.1, 23.3, 23.5, 23.7, 23.9, 24.1, 24.3, 24.5, 24.7,
               24.9, 25.1,  25.3, 25.5, 25.7, 25.9, 26.1, 26.3, 26.5, 26.7, 26.9, 27.1, 27.3, 27.5, 27.7, 27.9,
               28.1, 28.6, 29.1, 29.6, 30.1, 30.6, 31.1, 31.6, 32.1, 32.6, 33.1, 34.1, 34.9, 36.5, 37.85, 39.2, 40.7,
               42.3, 43.6, 45.1, 46.5, 47.9

               ])

r_fin = np.array([34.9, 36.5, 37.85, 39.2, 40.7,
               42.3, 43.6, 45.1, 46.5, 47.9])
V = np.array([ 4.44, 4.40, 4.05, 3.70,  3.94, 4.38, 4.02 , 3.63,  3.66, 3.70, 3.83, 4.15, 4.31,  4.01, 3.60,
              3.50, 3.59, 3.76, 4.00, 4.21, 4.03, 3.57, 3.37, 3.39, 3.54, 3.75, 3.98, 4.01, 3.62, 3.32, 3.26, 3.31,
               3.54, 3.77, 3.88, 3.58, 3.23,  3.17, 3.23, 3.40, 3.63, 3.73, 3.49, 3.18,  3.07, 3.12, 3.31,  3.48, 3.57,
               3.48, 3.01, 3.37, 3.22, 2.95, 3.25, 3.04, 2.88, 3.14, 2.87, 2.84, 2.73, 2.95, 2.85, 2.81, 2.76, 2.69,
               2.55, 2.48, 2.40, 2.3, 2.23

                ])
V_fin = np.array([ 2.95, 2.85, 2.81, 2.76, 2.69,
               2.55, 2.48, 2.45, 2.36, 2.32])

print(len(r_fin), len(V_fin))

# Lunghezza d'onda/2 = distanza tot / numero di massimi
picchi_idx , _ = find_peaks(V, height=3.1)  # trova i picchi sopra 3.1
picchi_r = r[picchi_idx]  # ottieni le posizioni dei picchi
picchi_V = V[picchi_idx]  # ottieni i valori dei picchi
r_fin = np.concatenate((r_fin, picchi_r))  # aggiungi i picchi alla lista delle posizioni
V_fin = np.concatenate((V_fin, picchi_V))  # aggiungi i picchi ai valori
print("Posizioni dei picchi:", picchi_r)
# Calcola la distanza tra i picchi
distanze = np.diff(picchi_r)
print("Distanze tra i picchi:", distanze)
# Calcola la media delle distanze tra i picchi
media_distanza = np.mean(distanze)
print("Media delle distanze tra i picchi:", media_distanza)
# Calcola la lunghezza d'onda
lunghezza_onda = media_distanza * 2
print("Lunghezza d'onda:", lunghezza_onda)
print(len(r), len(V))
plt.plot(r, V)
plt.plot(picchi_r, picchi_V, 'ro', label='Picchi trovati')
#plt.xlim(0,25)
#plt.ylim(4.0,4.5)
plt.show()
plt.scatter(r, V)
sigma_V = 0.05*np.ones(len(V_fin))  # errore gaussiano su V
lib.plot_fit(r_fin, V_fin, yerr=sigma_V, func=unosur, p0=[1, 1], confidence_intervals=False, prediction_band=False)
plt.show()

def function(x, A, B, C) :
  return A*np.abs(np.cos(x)) + B/(x**2)+C