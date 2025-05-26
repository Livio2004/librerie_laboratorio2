import analisi as lib
import numpy as np
import matplotlib.pyplot as plt
def calcola_errore_V (V) :
    sigma = 0.005*V + 2*0.1*10**(-3)
    b = V + sigma
    a = V - sigma
    sigma_gaus = ((b-a))/np.sqrt(12)
    return sigma_gaus

def calcola_errore_V_v2 (V) :
    sigma = 0.05*V + 2*0.1*10**(-3)
    b = V + sigma
    a = V - sigma
    sigma_gaus = ((b-a))/np.sqrt(12)
    return sigma_gaus

R1 = 32.8 # Ohm R1 = R2

Vin = 5.357 # V
sigma_Vin = calcola_errore_V(Vin)


Vout = np.array([ 0.555, 0.960, 1.498, 1.755, 2.025, 2.182, 2.322, 2.436, 2.495, 2.554, 2.586, 2.639, 2.651, 2.663, 2.672, 2.677, 2.683, 2.688    ])
RL = np.array([ 4, 10, 20, 30, 50, 70, 100, 150, 200, 300, 400, 800, 1000 , 1500, 2000, 3000,  5000, 10000   ])   #hannop errori
sigma_Vout = []
for i in Vout :
  sigma_Vout.append(calcola_errore_V_v2(i))

sigma_R = 0.05*np.ones(len(RL))

# errore Vin misura
print('sigma Vin : ', sigma_Vin)

# threshold RL accettabile:
Vin_mezzi_inf = Vin/2 - sigma_Vin/2

def funzione_partitore_resistivo(RL, R_usata):
  v_in = 5.357
  return RL*v_in/(R_usata+2*RL)
result = lib.fit(RL, Vout, yerr = sigma_Vout, function = funzione_partitore_resistivo, p0= [1], parameter_names=['R_usata'])
print(result['residuals'])
lib.plot_fit(RL, Vout, xerr = sigma_R,  yerr = sigma_Vout, func = funzione_partitore_resistivo, p0= [1], xlabel = 'RL',ylabel = 'Vout', title= 'Fit partitore resistivo' ,  parameter_names=['R_usata'], show_fit_params=True, show_chi_squared=True, confidence_intervals=True, error_band=1, prediction_band=True)
plt.axhline(y=Vin_mezzi_inf, color='r', linestyle='--', label='Vin/2-sigma/2')
plt.axhline(y=Vin/2, color='g', linestyle='--', label='Vin/2')
plt.show()





