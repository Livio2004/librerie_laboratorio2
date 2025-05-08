from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import analisi as lib
from scipy.optimize import curve_fit
from scipy.stats import norm
from iminuit import Minuit
from iminuit.cost import LeastSquares
R_L = 39.0
R = 993
Vin = 4

def calcola_errore_Ang (fase):
  sigma = 0.05*fase
  a = fase - sigma
  b = fase + sigma
  sigma_gaus = ((b-a))/np.sqrt(12)
  return sigma_gaus

def arg_HCL_f_RLC(omega, C, L):
    term1 = np.arctan(omega * C * R_L / (1 - omega**2 * L * C))
    term2 = np.arctan(omega * C * (R + R_L) / (1 - omega**2 * L * C))
    return term1 - term2

#DATI
from scipy import stats
frequenze_RLC_RL = np.array([  200, 500, 1000, 1300, 1600, 1800, 1950, 2050, 2150, 2250, 2350, 2400,  2450,  2690, 3000, 3300, 3600,  4600
 , 6100, 9300, 12300, 16500, 20100,   45100 ]) #
omega_RLC_HCL = 2*np.pi*frequenze_RLC_RL

'''sfase_ang_RLC_RL = np.array([  -4.32, -16.2,  -33.5,  -43.5,  -53, -59.6, -59.8, -70, -76 , 79 , 83 , 87 , 90, 63, 59.5 , 57.1 ,
                          52.9, 40.2, 30.3,  19.2, 15.5, 11,  8.95, 4.55 ])'''

sfase_grad_RLC_HCL = np.array([  -4.32, -16.2,  -33.5,  -43.5,  -53, -59.6, -63.8, -70, -76 , 79 , 73 , 71 , 69, 63, 59.2 , 54.1 ,
                          49.9, 40.0, 30.3,  19.2, 15.2, 11,  8.95, 4.55 ])

sfase_rad_RLC_HCL = sfase_grad_RLC_HCL*np.pi/180

#INCERTEZZE
sigma_grad_RLC_HCL = np.abs(calcola_errore_Ang (sfase_grad_RLC_HCL))

sigma_rad_RLC_HCL = sigma_grad_RLC_HCL*np.pi/180 * 2
popt, pcov = curve_fit(arg_HCL_f_RLC, omega_RLC_HCL, sfase_rad_RLC_HCL, p0=[10**(-7), 10**(-2)], sigma=sigma_rad_RLC_HCL, absolute_sigma=True)
C_fit, L_fit = popt
print(popt)
lib.plot_fit(omega_RLC_HCL, sfase_rad_RLC_HCL, yerr = sigma_rad_RLC_HCL, func= arg_HCL_f_RLC, p0= [10**(-7), 10**(-2)], prediction_band=True, confidence_intervals=True,
             xlabel = '$\omega$ [rad/s]', ylabel = 'fase $H_{CL}$', title = 'Circuito RLC in serie')
plt.show()
#######################################################################################################################################################

#attnzione il grande gemini



# Define the cost function (least squares)
least_squares = LeastSquares(omega_RLC_HCL, sfase_rad_RLC_HCL, sigma_rad_RLC_HCL, arg_HCL_f_RLC)

# Create a Minuit object with initial values for C and L
my_minuit = Minuit(least_squares, C=1e-7, L=1e-2)  # Adjust initial values if needed

# Perform the fit
my_minuit.migrad()  # Minimization
my_minuit.hesse()   # Error estimation

# Print the fit results
print("C =", my_minuit.values['C'])
print("L =", my_minuit.values['L'])

# Access fit parameters and errors
C_fit_phase = my_minuit.values['C']
L_fit_phase = my_minuit.values['L']
C_err = my_minuit.errors['C']
L_err = my_minuit.errors['L']

#Qsquared e altro
Q2 = my_minuit.fval
ndof = my_minuit.ndof
chi2_red = Q2 / ndof
p_value = 1 - stats.chi2.cdf(chi2_red, ndof)

# Print the results
print(f"C = ({C_fit_phase:.3e} +/- {C_err:.3e})")
print(f"L = ({L_fit_phase:.3e} +/- {L_err:.3e})")

'''asse_x = np.linspace(np.min(omega_RLC_RL), np.max(omega_RLC_RL), 1000)
asse_y = arg_HCL_f_RLC(asse_x, C_fit_phase, L_fit_phase)
plt.plot(asse_x, asse_y, label=f'fit, $\\chi^2_{{red}}$ = {chi2_red:.2f}, p-value = {p_value:.3f}', color='red')
plt.errorbar (omega_RLC_RL, sfase_rad_RLC_RL, sigma_rad_RLC_RL, linestyle = 'None', marker = 'o', capsize = 4, label = 'dati')
plt.legend()
plt.xscale('log')
plt.xlabel('$\omega$ [rad/s]')
plt.ylabel('fase $H_CL$')
plt.show()'''



##################################################################################################################################

x_fit = np.linspace(min(omega_RLC_HCL), max(omega_RLC_HCL), 1000)
y_fit = arg_HCL_f_RLC(x_fit, C_fit_phase, L_fit_phase)


z = stats.norm.ppf(0.8413)
#z = stats.norm.ppf(0.975)

C_conf_int = [C_fit_phase - z * C_err, C_fit_phase + z * C_err]
L_conf_int = [L_fit_phase - z * L_err, L_fit_phase + z * L_err]

cov_matrix = my_minuit.covariance
cov_CL = cov_matrix['C', 'L']

eps_C = C_err / 10
eps_L = L_err / 10

y_fit_C_plus = arg_HCL_f_RLC(x_fit, C_fit_phase + eps_C, L_fit_phase)
y_fit_C_minus = arg_HCL_f_RLC(x_fit, C_fit_phase - eps_C, L_fit_phase)
dy_dC = (y_fit_C_plus - y_fit_C_minus) / (2 * eps_C)* 5

y_fit_L_plus = arg_HCL_f_RLC(x_fit, C_fit_phase, L_fit_phase + eps_L)
y_fit_L_minus = arg_HCL_f_RLC(x_fit, C_fit_phase, L_fit_phase - eps_L)
dy_dL = (y_fit_L_plus - y_fit_L_minus) / (2 * eps_L) * 5

sigma_model = np.sqrt(
    (dy_dC * C_err)**2 +
    (dy_dL * L_err)**2 +
    2 * dy_dC * dy_dL * cov_CL
)

y_pred = y_fit
y_pred_upper = y_pred + z * sigma_model
y_pred_lower = y_pred - z * sigma_model

fig, ax = plt.subplots()

#plt.plot(x_fit, y_fit, label=f'fit, $\\chi^2_{{red}}$ = {chi2_red:.2f}, p-value = {p_value:.3f}', color='red')

plt.plot(
    x_fit,
    y_fit,
    label=(
        r"Fit" "\n"
        r"$\chi^2/n_{{\mathrm{{dof}}}}$ = "
        f"{Q2:.1f} / {ndof} = {chi2_red:.2f}\n"
        r"$p$-value = "
        f"{p_value:.3f}\n"
        r"$C$ = "
        f"{C_fit_phase:.3e} $\pm$ {C_err:.3e}\n"
        r"$L$ = "
        f"{L_fit_phase:.3f} $\pm$ {L_err:.3f}"
    ),
    color = 'red'
)


ax.fill_between(x_fit, y_pred_lower, y_pred_upper, color='green', alpha=0.3
              , label='Prediction band, $ 1(\sigma)$')
ax.errorbar(omega_RLC_HCL, sfase_rad_RLC_HCL, yerr=sigma_rad_RLC_HCL, fmt='o', label='Dati', capsize = 3, markersize = 2)
ax.set_ylabel('Fase [rad]')
ax.set_xlabel('ω [rad/s]')
ax.set_xscale('log')
ax.legend()
ax.legend(fontsize=9)


plt.show()


###################################################################################################################################

fig, ax = plt.subplots()

#plt.plot(x_fit, y_fit, label=f'fit, $\\chi^2_{{red}}$ = {chi2_red:.2f}, p-value = {p_value:.3f}', color='red')

plt.plot(
    x_fit,
    y_fit,
    label=(
        r"Fit" "\n"
        r"$\chi^2/n_{{\mathrm{{dof}}}}$ = "
        f"{Q2:.1f} / {ndof} = {chi2_red:.2f}\n"
        r"$p$-value = "
        f"{p_value:.3f}\n"
        r"$C$ = "
        f"{C_fit_phase:.3e} $\pm$ {C_err:.3e}\n"
        r"$L$ = "
        f"{L_fit_phase:.3f} $\pm$ {L_err:.3f}"
    ),
    color = 'red'
)


ax.fill_between(x_fit, y_pred_lower, y_pred_upper, color='green', alpha=0.3
              , label='Prediction band, $ 1(\sigma)$')
ax.errorbar(omega_RLC_HCL, sfase_rad_RLC_HCL, yerr=sigma_rad_RLC_HCL, fmt='o', label='Dati', capsize = 3, markersize = 2)
ax.set_ylabel('Fase [rad]')
ax.set_xlabel('ω [rad/s]')
ax.legend()
ax.legend(fontsize=9)

plt.show()


###################################################################################################################################



residui = sfase_rad_RLC_HCL - arg_HCL_f_RLC(C_fit_phase, L_fit_phase, omega_RLC_HCL)

plt.errorbar(
    omega_RLC_HCL, residui, yerr=sigma_rad_RLC_HCL,
    fmt='o', capsize=2, color='black', markersize = 2
)
plt.axhline(0, color='black', linestyle='--')
plt.xlabel('$\omega $[rad/s]')
plt.ylabel('Residui [rad]')
plt.grid(True)
plt.tight_layout()
plt.show()