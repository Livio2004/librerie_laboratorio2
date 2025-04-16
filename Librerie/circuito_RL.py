import numpy as np 
import matplotlib.pyplot as plt
import analisi as lib
from scipy.optimize import curve_fit

R_L = 39.0
R = 993
Vin = 4


def mod_HL_f_RL(omega, L): # R noto, RL noto
    num = np.sqrt(R_L**2 + omega**2 * L**2)
    den2 = (R+R_L)**2 + omega**2 * L**2
    return num/np.sqrt(den2)

def arg_HL_f_RL(omega,L):
    return np.arctan(omega*L/R_L) -np.arctan(omega*L / (R + R_L))

def arg_HL_f_RL_2(omega,L):
    return np.arctan(omega*L) -np.arctan(omega*L / (R))


# fit mod HL (Circuito RL)
frequenza_L_2 = np.array([ 200, 380, 715, 1350, 2550, 4800, 9050, 17000, 32300, 61000, 115000  ])

V_L_2 = np.array([ 220*10**(-3), 480*10**(-3), 864*10**(-3), 1.46, 2.40, 3.32, 3.80, 3.90, 3.96, 3.98, 4  ])
sigma_HL_2 = 2*np.array ([2, 1.6, 1.6, 8, 6, 6, 6, 16, 16, 20, 20  ])*10**(-2)/np.sqrt(12) # millivolt

sfase_ang_2_RL = np.array([  90 , 82, 78, 65, 51.3, 35, 21.2,  11, 3 ,1,  0  ])
sigma_ang_2_RL = np.array([ 6, 1.5, 1.6,  1.4, 1.4, 1.4, 1.5, 1.6, 1.8, 1.5, 1.9   ])*2/np.sqrt(12)*(np.pi/180)

omega_L_2 = 2*np.pi*frequenza_L_2
mod_HL_2 = V_L_2 / Vin
arg_HL_2 = sfase_ang_2_RL * np.pi / 180



result = lib.fit(omega_L_2, mod_HL_2, function= mod_HL_f_RL, yerr = sigma_HL_2,  p0= [10**(-2)], parameter_names=['L']    )
L_mod_HL= result['parameters'][0]
sigma_L_mod_HL = result['covariance'][0][0]
print(result['parameters'])
lib.plot_fit(omega_L_2, mod_HL_2, yerr = sigma_HL_2, func= mod_HL_f_RL, p0= [10**(-2)], confidence_intervals=True, prediction_band=True,
             xlabel = '$\omega$ [rad/s]', ylabel = 'modulo $H_L$')
plt.xscale('log')
plt.show()


#fit con resistenza parassita che non viene
result = lib.fit(omega_L_2, arg_HL_2, function= arg_HL_f_RL, yerr =sigma_ang_2_RL ,  p0= [0.04], parameter_names=['L']    )
print(result['parameters'])
#L_arg_HL= result['parameters'][0]
#sigma_L_arg_HL = result['covariance'][0][0]
lib.plot_fit(omega_L_2, arg_HL_2, yerr =sigma_ang_2_RL , func= arg_HL_f_RL, p0= [0.07], prediction_band=True, confidence_intervals=True,
             xlabel = '$\omega$ [rad/s]', ylabel = 'fase $H_L$')
plt.xscale('log')
asse_x = np.linspace(np.min(omega_L_2), np.max(omega_L_2), 1000)
asse_y = arg_HL_f_RL(asse_x, 0.050)
plt.plot(asse_x, asse_y)
plt.show()

L = 0.049
def arg_HL_f_RL(omega,R_L):
    return np.arctan(omega*L/R_L) -np.arctan(omega*L / (R + R_L))

lib.plot_fit(omega_L_2, arg_HL_2, func= arg_HL_f_RL, yerr =sigma_ang_2_RL ,  p0= [1], parameter_names=['R_L']  , prediction_band=True, confidence_intervals=True)
plt.xscale('log')
plt.show()

def arg_HL_f_RL(omega,R):
    return np.arctan(omega*L/R_L) -np.arctan(omega*L / (R + R_L))

lib.plot_fit(omega_L_2, arg_HL_2, func= arg_HL_f_RL, yerr =sigma_ang_2_RL ,  p0= [1000], parameter_names=['R']  , prediction_band=True, confidence_intervals=True)
plt.xscale('log')
plt.show()

def arg_HL_f_RL(omega,R, R_L):
    return np.arctan(omega*L/R_L) -np.arctan(omega*L / (R + R_L))

lib.plot_fit(omega_L_2, arg_HL_2, func= arg_HL_f_RL, yerr =sigma_ang_2_RL ,  p0= [1000, 10], parameter_names=['R', 'R_L']  , prediction_band=True, confidence_intervals=True)
plt.xscale('log')
plt.show()







