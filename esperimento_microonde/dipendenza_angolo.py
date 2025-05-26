import analisi as lib
import numpy as np
from iminuit import Minuit
from iminuit.cost import LeastSquares
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def calcola_errore_V(V): #errore gaussiano su V
    sigma = V*0.08
    b = V + sigma
    a = V - sigma
    sigma_gaus = ((b-a))/np.sqrt(12)
    return sigma_gaus



def onda_piana(θ,A):
    return A * np.cos(θ)

def onda_sferica(θ,A):
    return A

def modello(x, I0, a):
    return I0 * (a + (1 - a) * np.cos(x)**2)

#angolo alpha: 0
# dist= 41cm
theta = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90    ])
V_T_1 = np.array([3.53, 3.13, 2.09, 1.048, 0.442, 0.180, 0.035, 0.014, 0.013,  0.0026     ])
theta_2 = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, -10, -20, -30 , -40, -50, -60, -70,-80, -90 ])
V_T_2 = np.array([3.52 ,3.24 ,2.05 ,0.680 ,0.418 ,0.139,0.0217 ,0.0152 ,0.0131 ,0.03, 3.52, 2.617, 1.604, 0.730 , 0.270, 0.101, 0.008, 0.005, 0.002])
theta_1_rad = np.radians(theta)
theta_2_rad = np.radians(theta_2)
sigma_ang_1 = np.ones(len(theta_1_rad))*0.05
sigma_ang_2 = np.ones(len(theta_2_rad))*0.05
sigma_V_T_1 = np.maximum(calcola_errore_V(V_T_1), 0.05)
sigma_V_T_2 = np.maximum(calcola_errore_V(V_T_2), 0.05)
#sigma_V_T_1 = np.ones(len(V_T_1)) * 0.05
#sigma_V_T_2 = np.ones(len(V_T_2)) * 0.05
print(len(theta_2))
print(len(V_T_2))
def modello2(x, A, B):
    return A * np.cos(x)**2 + B * np.abs(np.cos(x))

lib.plot_fit(theta_2_rad, V_T_2,xerr = sigma_ang_2,  yerr = sigma_V_T_2, func= modello, p0= [10, 0.5] , confidence_intervals=False, prediction_band=False)
plt.show()
lib.plot_fit(theta_2_rad, V_T_2, xerr = sigma_ang_2,  yerr = sigma_V_T_2, func= modello2, p0= [10, 0.5] , confidence_intervals=False, prediction_band=False)

def modello3(theta, V0, n):
    return V0 * np.cos(theta)**n

lib.plot_fit(theta_2_rad, V_T_2, xerr = sigma_ang_2, yerr =sigma_V_T_2, func= modello3, p0= [4, 2] , confidence_intervals=False, prediction_band=False)
plt.show()

def modello4(theta, A, sigma):
  return A * np.exp( - (theta)**2 / (2*sigma**2) )


lib.plot_fit(theta_1_rad, V_T_1, yerr = sigma_V_T_1, func= modello3, p0= [4, 2] , confidence_intervals=False, prediction_band=False)
plt.show()

#fascio altamente collimato con pochi gradi di apertura
#oppure può essere collimato con divergenza gaussiana 
#non era molto grande la distanza tra emettitore e ricevitore e quindi magari può influire su come si vede la divergenza