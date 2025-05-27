import analisi as lib
import numpy as np
from iminuit import Minuit
from iminuit.cost import LeastSquares
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# emettitore 45grafi rispeto al centro cubo
#

# prima misuraz.

angolo_bragg = np.radians(np.array([ 100, 90, 87.5,  85, 82.5, 81, 80, 79 ,  77.5,  75, 72.5,  70, 60, 50  ])) # rispetto theta i
V_bragg = np.array([ 0.040, 0.746, 0.83, 1.8, 2.58, 2.76, 2.74, 2.49, 2.18, 1.96, 1.66, 1, 0.280,  0.800    ]) # V

d = 0.03825 # distanza centri reticolo cristallino TODO: controllare (15.3/4)
sigmad = 0.001 # m

plt.plot(angolo_bragg, V_bragg  )
plt.show()

# calcolo lambda:
n = 1  # primo ordine

# Trova l'angolo di Bragg corrispondente al massimo segnale
i_max = np.argmax(V_bragg)
theta_max = angolo_bragg[i_max]

# Calcola lambda
lambda_ = 2 * d * np.sin(theta_max) / n

# Propagazione dell'errore su lambda (solo errore su d, trascuriamo errore su theta per ora)
# dλ/dd = 2 * sin(theta)
sigma_lambda = 2 * np.sin(theta_max) * sigmad

# Risultato
print(f"Angolo di Bragg massimo: {np.degrees(theta_max):.2f}°")
print(f"Lunghezza d'onda λ: {lambda_*1e9:.2f} nm ± {sigma_lambda*1e9:.2f} nm")
# emettitore 45gradi rispeto al centro cubo
# distanze 100-30
#3x multiplier

# prima misuraz.

angolo_bragg = np.radians(np.array([ 100, 95, 90,  86, 85,84, 83, 82, 81,  80,79,78,77,76, 75, 70, 65 , 60 ,  55, 50,45  ])) # rispetto theta i
V_bragg = np.array([ 0.440, 2.045, 1.08,3.33, 3.50,3.89, 4.04, 4.17, 4.37,  4.31,4.50, 4.36,4.02,3.80, 3.07, 1.589, 1.403 ,1.180, 1.57 , 1.430, 1.280 ]) # V
sigma_V_b = 0.05*np.ones(len(V_bragg))
angolo_bragg_parabola = np.radians(np.array([84, 83, 82, 81,  80,79,78,77,76]))
V_bragg_parabola = np.array([3.89, 4.04, 4.17, 4.37,  4.31,4.50, 4.36,4.02,3.80])
sigma_V_b_p = 0.05*np.ones(len(V_bragg_parabola))
plt.plot(angolo_bragg_parabola, V_bragg_parabola, label='Dati sperimentali')
def parabola(x, a, b, c):
    return a * (x - b)**2 + c
lib.fit(angolo_bragg_parabola, V_bragg_parabola, function=parabola, p0=[-1, 0.5, 4], xerr=None, yerr=sigma_V_b_p)
lib.plot_fit(angolo_bragg_parabola, V_bragg_parabola, func=parabola, p0=[-1, 0.5, 4], xerr=None, yerr=sigma_V_b_p, confidence_intervals=False, prediction_band=False)
plt.show()
d = 0.03825 # distanza centri reticolo cristallino TODO: controllare (15.3/4)
sigmad = 0.001 # m

plt.plot(angolo_bragg, V_bragg  )
plt.show()