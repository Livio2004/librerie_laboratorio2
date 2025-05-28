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


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# --- Dati forniti ---
angolo_bragg_gradi = np.array([100, 95, 90, 86, 85, 84, 83, 82, 81, 80, 79, 78, 77, 76, 75, 70, 65, 60, 55, 50, 45])
V_bragg = np.array([0.440, 2.045, 1.08, 3.33, 3.50, 3.89, 4.04, 4.17, 4.37, 4.31, 4.50, 4.36, 4.02, 3.80, 3.07, 1.589, 1.403, 1.180, 1.57, 1.430, 1.280])
sigma_V_b = 0.05 * np.ones(len(V_bragg))
angolo_bragg_rad = np.radians(angolo_bragg_gradi)

# Dati selezionati per fit parabolico (zona del massimo)
# È importante che questi dati selezionati centrino bene il picco
# Potresti dover aggiustare questi indici o valori se il picco è altrove
angolo_bragg_parabola_gradi = np.array([84, 83, 82, 81, 80, 79, 78, 77, 76]) # Gradi
V_bragg_parabola = np.array([3.89, 4.04, 4.17, 4.37, 4.31, 4.50, 4.36, 4.02, 3.80])
sigma_V_b_p = 0.05 * np.ones(len(V_bragg_parabola))
angolo_bragg_parabola_rad = np.radians(angolo_bragg_parabola_gradi)

# --- Fit parabolico per trovare il massimo ---
# y = a*(x - b)^2 + c
# b è la posizione x del massimo, c è l'altezza del massimo
def parabola(x, a, b, c):
    return a * (x - b)**2 + c

# Stime iniziali per il fit:
# a: curvatura (negativa per un picco)
# b: posizione approssimativa del picco in radianti
# c: altezza approssimativa del picco
# Controlla che np.radians(80) sia una buona stima per il tuo picco
p0_parabola = [-10, np.radians(80), 4.5] # Ho aggiustato p0 per a e c basandomi sui dati

try:
    popt_parabola, pcov_parabola = curve_fit(
        parabola,
        angolo_bragg_parabola_rad,
        V_bragg_parabola,
        p0=p0_parabola,
        sigma=sigma_V_b_p,
        absolute_sigma=True
    )
    a_fit, b_fit_rad, c_fit = popt_parabola
    sigma_b_rad = np.sqrt(pcov_parabola[1, 1])  # errore sull'angolo del massimo in radianti

    theta_max_rad = b_fit_rad
    sigma_theta_max_rad = sigma_b_rad

    print(f"--- Risultati del Fit Parabolico ---")
    print(f"Parametri del fit (a, b_rad, c): {popt_parabola}")
    print(f"Angolo del massimo (dal fit): {np.degrees(theta_max_rad):.2f}° ± {np.degrees(sigma_theta_max_rad):.2f}°")
    print(f"Tensione al massimo (dal fit): {c_fit:.3f} V")

    # --- Calcolo del passo del reticolo d ---
    lambda_nota = 0.0315  # m (3.1 cm)
    # Assumiamo che lambda_nota sia esatta o con incertezza trascurabile per questo calcolo
    # Se lambda_nota avesse un'incertezza significativa, sigma_lambda_nota, andrebbe propagata.
    sigma_lambda_nota = 0 # Assumiamo esatta

    n_ordine = 1 # Assumiamo primo ordine di diffrazione

    # Legge di Bragg: n * lambda = 2 * d * sin(theta)
    # d = (n * lambda) / (2 * sin(theta))
    d_calcolato = (n_ordine * lambda_nota) / (2 * np.sin(theta_max_rad))

    # Propagazione dell'incertezza per d
    # d = f(theta_max_rad, lambda_nota)
    # sigma_d^2 = (dd/d(theta))^2 * sigma_theta^2 + (dd/d(lambda))^2 * sigma_lambda^2
    # dd/d(theta) = - (n * lambda_nota * cos(theta_max_rad)) / (2 * sin(theta_max_rad)^2)
    # dd/d(lambda) = n / (2 * sin(theta_max_rad))

    parziale_d_theta = - (n_ordine * lambda_nota * np.cos(theta_max_rad)) / (2 * np.sin(theta_max_rad)**2)
    parziale_d_lambda = n_ordine / (2 * np.sin(theta_max_rad))

    sigma_d_calcolato = np.sqrt(
        (parziale_d_theta**2 * sigma_theta_max_rad**2) + \
        (parziale_d_lambda**2 * sigma_lambda_nota**2)
    )
    # Se sigma_lambda_nota è 0, il secondo termine svanisce:
    # sigma_d_calcolato = np.abs(parziale_d_theta) * sigma_theta_max_rad

    print(f"\n--- Calcolo del Passo del Reticolo d ---")
    print(f"Lunghezza d'onda nota (lambda): {lambda_nota*100:.2f} cm")
    print(f"Ordine di diffrazione (n): {n_ordine}")
    print(f"Angolo di Bragg utilizzato (theta_max): {np.degrees(theta_max_rad):.2f}° ± {np.degrees(sigma_theta_max_rad):.2f}°")
    print(f"Passo del reticolo calcolato (d): {d_calcolato*100:.3f} cm ± {sigma_d_calcolato*100:.3f} cm")
    print(f"Passo del reticolo calcolato (d): {d_calcolato*1e3:.2f} mm ± {sigma_d_calcolato*1e3:.2f} mm")


    # --- Plot ---
    # Genera punti per la curva di fit parabolico
    theta_plot_rad = np.linspace(min(angolo_bragg_parabola_rad) * 0.95, max(angolo_bragg_parabola_rad) * 1.05, 500)
    V_plot_fit = parabola(theta_plot_rad, *popt_parabola)

    plt.figure(figsize=(10, 6))
    plt.errorbar(angolo_bragg_rad, V_bragg, yerr=sigma_V_b, fmt='o', capsize=3, label='Tutti i dati sperimentali', alpha=0.6)
    plt.errorbar(angolo_bragg_parabola_rad, V_bragg_parabola, yerr=sigma_V_b_p, fmt='s', color='orange', capsize=3, label='Dati usati per il fit')
    plt.plot(theta_plot_rad, V_plot_fit, label=f'Fit parabolico (per trovare θ_max)', color='red', linestyle='-')
    plt.axvline(theta_max_rad, color='green', linestyle='--', label=f'Massimo θ_max: {np.degrees(theta_max_rad):.2f}°')

    plt.xlabel("Angolo di Bragg (radianti)")
    plt.ylabel("Tensione (V)")
    plt.title(f"Fit Parabolico e Stima di d (λ = {lambda_nota*100:.1f} cm)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

except RuntimeError:
    print("Fit parabolico non riuscito. Controlla i dati selezionati e i parametri iniziali (p0).")
    plt.figure(figsize=(10, 6))
    plt.errorbar(angolo_bragg_rad, V_bragg, yerr=sigma_V_b, fmt='o', capsize=3, label='Tutti i dati sperimentali', alpha=0.6)
    plt.errorbar(angolo_bragg_parabola_rad, V_bragg_parabola, yerr=sigma_V_b_p, fmt='s', color='orange', capsize=3, label='Dati usati per il fit (tentativo)')
    plt.xlabel("Angolo di Bragg (radianti)")
    plt.ylabel("Tensione (V)")
    plt.title("Dati Bragg - Fit Parabolico Fallito")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()