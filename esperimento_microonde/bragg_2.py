import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Tutti i dati
angolo_bragg_2 = np.radians(np.array([ 100, 90, 87.5,  85, 82.5, 81, 80, 79 ,  77.5,  75, 72.5,  70, 60, 55, 52.5, 50, 47.5, 45  ])) # rispetto theta i
V_bragg_2 = np.array([ 0.040, 0.746, 0.83, 1.8, 2.58, 2.76, 2.74, 2.49, 2.18, 1.96, 1.66, 1, 0.280,  0.6, 0.87, 0.79, 0.58, 0.57      ]) # V
sigma_V_b = 0.05 * np.ones(len(V_bragg_2))

# Dati selezionati per fit parabolico (zona del massimo)
angolo_bragg_parabola = np.radians(np.array([85, 82.5, 81, 80, 79 ,  77.5]))
V_bragg_parabola = np.array([1.8, 2.58, 2.76, 2.74, 2.49, 2.18])
sigma_V_b_p = 0.05 * np.ones(len(V_bragg_parabola))

# Fit parabolico: y = a*(x - b)^2 + c
def parabola(x, a, b, c):
    return a * (x - b)**2 + c

popt, pcov = curve_fit(parabola, angolo_bragg_parabola, V_bragg_parabola, p0=[-1, np.radians(80), 4], sigma=sigma_V_b_p, absolute_sigma=True)
a_fit, b_fit, c_fit = popt
sigma_b = np.sqrt(pcov[1, 1])  # errore sul massimo

# Calcolo lambda al massimo della parabola
d = 0.03825  # m
sigmad = 0.001  # m
n = 1
theta_max = b_fit
lambda_ = 2 * d * np.sin(theta_max) / n
sigma_lambda = 2 * np.sin(theta_max) * sigmad  # trascurando incertezza su theta

# Output
print(f"Massimo (angolo fit): {np.degrees(theta_max):.2f}° ± {np.degrees(sigma_b):.2f}°")
print(f"Lunghezza d'onda λ: {lambda_*1e9:.2f} nm ± {sigma_lambda*1e9:.2f} nm")

# Plot
theta_fit = np.linspace(np.radians(75), np.radians(85), 500)
V_fit = parabola(theta_fit, *popt)

plt.errorbar(angolo_bragg_2, V_bragg_2, yerr=sigma_V_b, fmt='o', label='Tutti i dati')
plt.errorbar(angolo_bragg_parabola, V_bragg_parabola, yerr=sigma_V_b_p, fmt='o', color='orange', label='Dati per fit')
plt.plot(theta_fit, V_fit, label='Fit parabolico', color='red')
plt.axvline(theta_max, color='green', linestyle='--', label=f'Massimo: {np.degrees(theta_max):.2f}°')

plt.xlabel("Angolo (rad)")
plt.ylabel("Tensione (V)")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
