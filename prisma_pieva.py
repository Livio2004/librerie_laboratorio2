# ...
import numpy as np
import matplotlib.pyplot as plt

import sys
print(sys.path)

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from librerie_laboratorio2.Librerie import analisi as lib

from iminuit import Minuit
from iminuit.cost import LeastSquares
from IPython.display import display, Latex
from librerie_laboratorio2.spettrometro import spectra_reppresentation as spettro

def indice_rifrazione(δm,alpha): # α = 60.18, δm = angolo minima deviazione
    return np.sin((δm + alpha)/2 * np.pi/180) / np.sin(alpha/2 *np.pi/180)

def sigma_indice_rifrazione(δm, alpha, sigma_δm, sigma_alpha):
    rad = np.pi / 180  # conversione gradi → radianti

    # Parte comune
    num = np.sin((δm + alpha)/2 * rad)
    den = np.sin(alpha/2 * rad)
    n = num / den

    # Derivate parziali
    dn_dδm = (0.5 * rad * np.cos((δm + alpha)/2 * rad)) / den
    dn_dalpha = (
        (0.5 * rad * np.cos((δm + alpha)/2 * rad)) / den -
        (num * 0.5 * rad * np.cos(alpha/2 * rad)) / (den**2)
    )

    # Propagazione dell'errore
    sigma_n = np.sqrt(
        (dn_dδm * sigma_δm)**2 +
        (dn_dalpha * sigma_alpha)**2
    )
    
    return sigma_n

# PRIMA GIORNATA

# stima del valore di alpha

alpha = 60
sigma_alpha = 0

# SECONDA GIORNATA

angolo_centr = 73

array_angoli = []
array_sigma = []

# orange 1 forte

orange = 123.8 - angolo_centr
array_angoli.append(orange)
array_sigma.append(0.10)

# green

green = 124.3  - angolo_centr
array_angoli.append(green)
array_sigma.append(0.1)

# azzurro meno forte quasi verde

azver = 125.4 - angolo_centr
array_angoli.append(azver)
array_sigma.append(0.12)

# azzurro più forte (cmq nn ben visibile)

#azurr = 125.7 - angolo_centr
#array_angoli.append(azurr)

# blu molto bene

blu = 127 - angolo_centr
array_angoli.append(blu )
array_sigma.append(0.17)

# viola meno visibile

viola = 128.5 - angolo_centr
array_angoli.append(viola )
array_sigma.append(0.09)

# viola più visibile

viola2 = 128.7 - angolo_centr
array_angoli.append(viola2 )
array_sigma.append(0.08)

array_n3 = []

for j in (array_angoli) :
  array_n3.append(indice_rifrazione(j, 60))
print(array_n3)

# formula di Cauchy
def n_cauchy (lambd, a , b):
  return a + b/np.pow(lambd, 2)

# Lunghezze d'onda NIST per Hg (in nm) BASATE SULLA TABELLA E INTENSITÀ
# Questa associazione è CRUCIALE e va verificata con gli appunti di lab
array_lambda_calcolate_nm = np.array([
    579.0,  # Per "orange": la più lunga delle gialle intense (o 576.9 o una media). 
             # Se hai visto distintamente arancione e poi giallo, questa tabella è incompleta per le tue osservazioni.
             # Se "orange" era la tua percezione del giallo più deviato (lambda maggiore)
    546.0,  # Per "green": Verde Hg
    491.6,  # Per "azver": Blu/Verde Hg (segnata come debole, ma la tua successiva è più deviata) [La scarto dal fit perchè deboel]
             # Se questa non fosse visibile, allora azver e blu dovrebbero essere molto vicine alla 435.8 nm
    435.8,  # Per "blu": Azzurro/blu Hg (Intensa)
    407.7,  # Per "viola": Viola Hg (Media)
    404.6   # Per "viola2": Viola Hg (Media, più deviata)
])

lam_3 = np.array(array_lambda_calcolate_nm) # *10**(-9)
n_3 =  array_n3


sigma_n = []
for i in range (len(array_angoli)):
    sigma_n.append(sigma_indice_rifrazione(array_angoli[i], alpha, array_sigma[i], sigma_alpha))


print('n: ',np.round(array_n3,5),' ± ',np.round(sigma_n,5))

results_r2 = lib.fit(lam_3, n_3,  yerr = sigma_n,
                  function= n_cauchy, p0 = [1.61665, 10000.8],  )
lib.plot_fit(lam_3, n_3,  yerr = sigma_n,
                  func= n_cauchy, p0 = [1.61665, 10000.8], confidence_intervals = False,  
                  xlabel = ' $\lambda$ (m)', ylabel= 'indice rifrazione $n$')
plt.title('fit finale')
plt.show()

a_fit2 = results_r2['parameters'][0] 
b_fit2 = results_r2['parameters'][1] 
print('a_fit2, b_fit2: ', a_fit2, b_fit2)

sigma_a = np.sqrt(results_r2['covariance'][0][0])
sigma_b = np.sqrt(results_r2['covariance'][1][1])
rho_ab = (results_r2['covariance'][1][0])
print('rho ab: ', rho_ab)
# spectrum

spettro.draw_spectrum(array_lambda_calcolate_nm,ytitle='Hg')



# inverto formula

def ricava_lambda (a, b, n):
  return np.sqrt(b/np.abs(n-a))

def sigma_lambda(n, a, b, sigma_n, sigma_a, sigma_b, rho_ab):
    delta = n - a
    abs_delta = np.abs(delta)
    sgn = np.sign(delta)

    # lambda
    lam = np.sqrt(b / abs_delta)

    # derivate
    dlam_dn = -0.5 * b / (abs_delta**2 * lam) * sgn
    dlam_da =  0.5 * b / (abs_delta**2 * lam) * sgn
    dlam_db =  0.5 / np.sqrt(b * abs_delta)

    # varianza con termine di covarianza
    var_lam = (dlam_dn * sigma_n)**2 + (dlam_da * sigma_a)**2 + (dlam_db * sigma_b)**2
    var_lam += 2 * dlam_da * dlam_db * rho_ab * sigma_a * sigma_b

    return np.sqrt(var_lam)*10**(-9)


# GAS IGNOTO CON PRISMA

array_gas2 = []
array_deltam = []

indef = 83
indef = 80.1
# rosso secondario (meno intenso)

rosso =  120

n_rosso = indice_rifrazione(np.mean(rosso)-indef,  60 )
array_gas2.append(n_rosso)

#

rosso2 = 120.5
delta_rosso2 = 162.5-162 # oscillazione max min riga di colore

n_rosso2 = indice_rifrazione(np.mean(rosso2)-indef,  60 )
array_gas2.append(n_rosso2)

print(n_rosso2)

orange= 121

n_orange = indice_rifrazione(np.mean(orange)-indef,  60 )
array_gas2.append(n_orange)

print(n_orange)


verde_acqua = 121.8

n_verde_acqua= indice_rifrazione(np.mean(verde_acqua)-indef,  60 )
array_gas2.append(n_verde_acqua)

print(n_verde_acqua)

# verde acqua meno forte
azure = 122

n_azure= indice_rifrazione(np.mean(azure)-indef,  60 )
array_gas2.append(n_azure)

print(n_azure)

# blue meno forte

blu = 122.3
n_blue = indice_rifrazione(np.mean(blu )-indef ,  60 )

array_gas2.append(n_blue )
print(n_blue )

violet  = 122.7

n_violet = indice_rifrazione(np.mean(violet )-indef ,  60 )
array_gas2.append(n_violet )

print(n_violet )

array_deltam = np.array([ 118.5, 119.4, 120.4, 121.2, 122.7, 123, 123.2  ])-80.1 # 80.1

print('array angoli  gas ignoto: ', array_deltam)
print('array n  gas ignoto: ', array_gas2)


n_gas = [  ]        
sigman_gas = []      
λ_ignoto = []
sigmaλ_ignoto = []

for i in range (len(array_gas2)):
   n_gas.append(indice_rifrazione(array_deltam[i], alpha))
   sigman_gas.append(sigma_indice_rifrazione(array_deltam[i], 0.01, alpha, sigma_alpha))
   λ_ignoto.append(ricava_lambda(a=a_fit2, b=b_fit2, n=n_gas[i])) 
   sigmaλ_ignoto.append(sigma_lambda(n_gas[i], a_fit2, b_fit2, sigman_gas[i], sigma_a, sigma_b, rho_ab))

print('λ ignoto [nm]: ',  λ_ignoto, '+/-', sigmaλ_ignoto)
      
lambda_gas = []
for j in array_gas2 :
  lambda_gas.append(ricava_lambda(a_fit2, b_fit2, j))
  print('lambda :', ricava_lambda(a_fit2, b_fit2, j))
spettro.draw_spectrum(lambda_gas,ytitle='ign')


# --- CORREZIONE DELLA FUNZIONE sigma_lambda ---
def sigma_lambda_propagated(n_val, a_val, b_val, sigma_n, sigma_a, sigma_b, cov_ab):
    """
    Calcola l'incertezza su lambda = sqrt(b / |n-a|) propagando gli errori.
    Tutti i valori di input (n, a, b e le loro sigma, cov) devono essere in unità consistenti.
    Se b è in nm^2, lambda sarà in nm.
    """
    delta = n_val - a_val
    
    if np.isclose(delta, 0) and b_val !=0 :
        print(f"sigma_lambda: n-a ({delta}) è vicino a zero. Impossibile calcolare l'errore in modo stabile.")
        return np.inf # O un valore molto grande per indicare divergenza
    if np.isclose(delta,0) and np.isclose(b_val,0):
        print(f"sigma_lambda: n-a ({delta}) e b ({b_val}) sono vicini a zero.")
        return np.nan

    abs_delta = np.abs(delta)
    sgn_delta = np.sign(delta) # segno di (n-a)

    # Lambda calcolato
    if b_val / abs_delta < 0: # Dovrebbe essere già gestito da ricava_lambda
        print("sigma_lambda: b/|n-a| < 0, lambda non reale.")
        return np.nan
        
    lam = np.sqrt(b_val / abs_delta)
    if np.isclose(lam,0): # Se b è zero e delta no
        print("sigma_lambda: lambda calcolato è zero, le derivate potrebbero divergere.")
        # Per ora, se b_val è piccolo, dlam_db sarà grande.
        if np.isclose(b_val,0): # Se b è effettivamente zero o molto vicino
             if not np.isclose(delta,0): # e n != a, allora lambda = 0
                 # dlam_db in questo caso è problematico.
                 # Se b=0 (esatto), lambda=0. Se b = 0 +/- sigma_b, lambda = sqrt(sigma_b/abs_delta)
                 return np.sqrt(sigma_b / abs_delta) if abs_delta > 1e-9 else np.inf


    # Derivate parziali (controlla che lam non sia zero se b!=0)

    dlam_dn = -0.5 * np.sqrt(b_val) * np.power(abs_delta, -1.5) * sgn_delta
    
    dlam_da = -dlam_dn 
    dlam_db = 0.5 * lam / b_val if not np.isclose(b_val, 0) else np.inf # Gestisce b_val = 0

    # Varianza
    var_lam = (dlam_dn * sigma_n)**2 + \
              (dlam_da * sigma_a)**2 + \
              (dlam_db * sigma_b)**2 + \
              2 * dlam_da * dlam_db * cov_ab # Termine di covarianza per a,b

    # Nota: non c'è covarianza tra n e (a,b) perché n è misurato indipendentemente

    return np.sqrt(var_lam)*10**(-9)


for i in range (len(array_gas2)):
   n_gas.append(indice_rifrazione(array_deltam[i], alpha))
   sigman_gas.append(sigma_indice_rifrazione(array_deltam[i], 0.01, alpha, sigma_alpha))
   λ_ignoto.append(ricava_lambda(a=a_fit2, b=b_fit2, n=n_gas[i])) 
   sigmaλ_ignoto.append(sigma_lambda_propagated(n_gas[i], a_fit2, b_fit2, sigman_gas[i], sigma_a, sigma_b, rho_ab))

print('λ ignoto 2 [nm]: ',  λ_ignoto, '+/-', sigmaλ_ignoto)
  

spettro.draw_spectrum(λ_ignoto,ytitle='ign2')

