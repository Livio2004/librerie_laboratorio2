from iminuit import Minuit
from iminuit.cost import LeastSquares
import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
Vin = 4

C_fit4= 9.70914278e-08
L_fit4 = 0.048
R = 993
R_L = 39
def mod_HCL_f_RLC(omega, C, L):
    d = (1 - omega**2 * L * C)**2 + omega**2 * C**2 * R_L**2
    numerator = np.sqrt(d)
    denominator = np.sqrt((1 - omega**2 * L * C)**2 + omega**2 * C**2 * (R + R_L)**2)
    return numerator / denominator


frequenze_RLC_RL = np.array([  200, 500, 1000, 1300, 1600, 1800, 1950, 2050, 2150, 2250, 2350, 2400,  2450,  2690, 3000, 3300, 3600,  4600
 , 6100, 9300, 12300, 16500, 20100,   45100 ]) #

V_RLC_RL = np.array([ 4, 3.90, 3.24, 2.60, 1.92, 1.36, 1,  720*10**(-3), 520*10**(-3), 340*10**(-3), 310*10**(-3),
                     400*10**(-3), 480*10**(-3), 880*10**(-3), 1.44, 1.84,
                          2.16,  2.96, 3.52, 3.92, 4,  4,  4,  4 ]) # 4

#V_RLC_RL_norm = V_RLC_RL / np.max(V_RLC_RL) #oooooooooooooooooooooooooooooooooooo

sfase_ang_RLC_RL = np.array([  -4.32, -16.2,  -33.5,  -43.5,  -53, -59.6, -59.8, -70, -76 , 79 , 83 , 87 , 90, 63, 59.5 , 57.1 ,
                          52.9, 40.2, 30.3,  19.2, 15.5, 11,  8.95, 4.55 ])

# - 56 al posto di 70
omega_RLC_RL = 2*np.pi*frequenze_RLC_RL
#omega_RLC_RL = omega_RLC_RL - min(omega_RLC_RL)

mod_HCL = V_RLC_RL / Vin

sigma = np.ones(len(frequenze_RLC_RL))

sigma_RLC_HCL =  (2*0.1/np.sqrt(12))*sigma



#V_RLC_RL_norm = V_RLC_RL / np.max(V_RLC_RL) #oooooooooooooooooooooooooooooooooooo
#sigma_RLC_HCL_norm =  (2*0.1/np.sqrt(12))*sigma #ooooooooooooooooooooooooooo


arg_HCL = sfase_ang_RLC_RL*2*np.pi
sigma_ang_RLC_RL = 2*0.5/np.sqrt(12)*sigma*np.pi/180

least_squares = LeastSquares (omega_RLC_RL, mod_HCL, sigma_RLC_HCL, mod_HCL_f_RLC)
my_minuit = Minuit (least_squares, C = C_fit4, L = L_fit4)  # starting values for m and q
my_minuit.migrad () #trova il minimo della funzione least squared
my_minuit.hesse ()
C_fit = my_minuit.values ['C']
L_fit = my_minuit.values ['L']

fig, ax = plt.subplots ()

x_fit = np.linspace(min(omega_RLC_RL), max(omega_RLC_RL), 1000)
y_fit = np.array ([mod_HCL_f_RLC(x, C_fit, L_fit) for x in x_fit])

ax.plot(x_fit, y_fit)
ax.set_ylabel ('y')
ax.errorbar (omega_RLC_RL, mod_HCL, xerr = 0.0, yerr = sigma_RLC_HCL, linestyle = 'None', marker = 'o')
plt.xscale('log')
plt.show()



#con scipy

from scipy.optimize import curve_fit


popt, pcov = curve_fit(mod_HCL_f_RLC, omega_RLC_RL, mod_HCL, p0=[C_fit4, L_fit4], sigma=sigma_RLC_HCL, absolute_sigma=True)
C_fit, L_fit = popt
C_err, L_err = np.sqrt(np.diag(pcov))

fig, ax = plt.subplots()

x_fit = np.linspace(min(omega_RLC_RL), max(omega_RLC_RL), 1000)
y_fit = mod_HCL_f_RLC(x_fit, C_fit, L_fit)

ax.plot(x_fit, y_fit, label=f'Fit: C = {C_fit:.2e} F, L = {L_fit:.2e} H')
ax.errorbar(omega_RLC_RL, mod_HCL, yerr=sigma_RLC_HCL, fmt='o', label='Dati')
ax.set_xlabel('ω [rad/s]')
ax.set_ylabel('Modulo normalizzato')
ax.set_xscale('log')
ax.legend()
plt.show()
