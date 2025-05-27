import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Tutti i dati
angolo_bragg = np.radians(np.array([100, 95, 90, 86, 85, 84, 83, 82, 81, 80, 79, 78, 77, 76, 75, 70, 65, 60, 55, 50, 45]))
V_bragg = np.array([0.440, 2.045, 1.08, 3.33, 3.50, 3.89, 4.04, 4.17, 4.37, 4.31, 4.50, 4.36, 4.02, 3.80, 3.07, 1.589, 1.403, 1.180, 1.57, 1.430, 1.280])
sigma_V_b = 0.05 * np.ones(len(V_bragg))

# Dati selezionati per fit parabolico (zona del massimo)
angolo_bragg_parabola = np.radians(np.array([84, 83, 82, 81, 80, 79, 78, 77, 76]))
V_bragg_parabola = np.array([3.89, 4.04, 4.17, 4.37, 4.31, 4.50, 4.36, 4.02, 3.80])
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

plt.errorbar(angolo_bragg, V_bragg, yerr=sigma_V_b, fmt='o', label='Tutti i dati')
plt.errorbar(angolo_bragg_parabola, V_bragg_parabola, yerr=sigma_V_b_p, fmt='o', color='orange', label='Dati per fit')
plt.plot(theta_fit, V_fit, label='Fit parabolico', color='red')
plt.axvline(theta_max, color='green', linestyle='--', label=f'Massimo: {np.degrees(theta_max):.2f}°')

plt.xlabel("Angolo (rad)")
plt.ylabel("Tensione (V)")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()