import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import scipy.stats as stats
from scipy.optimize import curve_fit
from scipy.odr import RealData, Model, ODR
from typing import Union, List, Dict, Callable, Optional, Tuple
from scipy.integrate import quad
from scipy.stats import chi2, t, norm
from jacobi import propagate



#------------propagazione errore----------------
'''
    formula = "log(x)+exp(y)"
    variables = ["x", "y"]
    errors = [1, 2]
    values = [10,4] 
'''

def error_propagation(formula, variables, errors, values):
    sym_variables = [sp.symbols(var) for var in variables]
    expr = sp.sympify(formula)
    partial_derivates = [sp.diff(expr, var) for var in sym_variables]
    partial_derivates_num = [deriv.subs(zip(sym_variables, values)) for deriv in partial_derivates]

    error_squared_formula= sum((deriv*error)**2 for deriv, error in zip(partial_derivates, errors))
    errore_propagato_formula = sp.sqrt(error_squared_formula)

    error_squared_num = sum((deriv*error)**2 for deriv, error in zip(partial_derivates_num, errors))
    errore_propagato_num = sp.sqrt(error_squared_num)
    return errore_propagato_formula, errore_propagato_num.evalf(4)


#------------fit----------------

def fit(x: np.ndarray, 
        y: np.ndarray, 
        function: Callable, 
        xerr: Optional[Union[List[float], np.ndarray]] = None, 
        yerr: Optional[Union[List[float], np.ndarray]] = None, 
        p0: Optional[Union[float, List[float], np.ndarray]] = None, 
        chi_square: bool = True, 
        method: str = 'auto', 
        parameter_names: Optional[List[str]] = None) -> Dict:
    
    if p0 is not None:
        p0 = np.atleast_1d(p0)

    x = np.asarray(x)
    y = np.asarray(y)
    if xerr is not None:
        xerr = np.asarray(xerr)
    if yerr is not None:
        yerr = np.asarray(yerr)
    else:
        yerr = np.zeros_like(y)

    if method == 'auto':
        if xerr is not None and yerr is not None:
            method = 'odr'
        else:
            method = 'curve_fit'

    if method == 'curve_fit':
        # Approssima la derivata per propagare l'errore di x in y
        if xerr is not None:
            def df_dx(x_val, *params):
                dx = 1e-8
                return (function(x_val + dx, *params) - function(x_val - dx, *params)) / (2 * dx)

            deriv = df_dx(x, *p0)
            sigma_eff = np.sqrt(yerr**2 + (deriv * xerr)**2)
        else:
            sigma_eff = yerr

        popt, pcov = curve_fit(function, x, y, p0=p0, sigma=sigma_eff, absolute_sigma=True)
        residuals = y - function(x, *popt)
        chi2 = np.sum((residuals / sigma_eff) ** 2) if chi_square else None

    elif method == 'odr':
        model = Model(lambda beta, x: function(x, *beta))
        data = RealData(x, y, sx=xerr, sy=yerr)
        odr = ODR(data, model, beta0=p0)
        output = odr.run()
        popt = output.beta
        pcov = output.cov_beta
        residuals = y - function(x, *popt)
        chi2 = output.res_var * len(x) if chi_square else None

    else:
        raise ValueError("Invalid method. Choose 'curve_fit', 'odr', or 'auto'.")

    result = {
        'parameters': popt,
        'covariance': pcov,
        'residuals': residuals,
        'chi2': chi2
    }

    if parameter_names:
        result['parameters_named'] = dict(zip(parameter_names, popt))

    return result

#------------grafico fit----------------

def plot_fit(
    x: np.ndarray,
    y: np.ndarray,
    xerr: Optional[np.ndarray] = None,
    yerr: Optional[np.ndarray] = None,
    func: Callable = None,
    p0: Optional[List[float]] = None,
    xlabel: str = '',
    ylabel: str = '',
    title: str = '',
    fit_line: bool = True,
    figsize: Tuple[int, int] = (10, 6),
    save_path: Optional[str] = None,
    dpi: int = 300,
    grid: bool = False,
    legend_loc: str = 'best',
    fit_points: int = 1000,
    fit_range: Optional[Tuple[float, float]] = None,
    residuals: bool = False,
    fmt: str = '+',
    markersize: int = 6,
    linewidth: float = 1.5,
    show_fit_params: bool = True,
    show_chi_squared: bool = True,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    capsize: float = 3,
    axis_fontsize: int = 12,
    title_fontsize: int = 14,
    parameter_names: Optional[List[str]] = None,
    method: str = 'auto',
    confidence_intervals: bool = False,
    prediction_band: bool = False,
    error_band: int = 1
):
    result = fit(x, y, func, xerr=xerr, yerr=yerr, p0=p0, method=method, parameter_names=parameter_names)

    fig, ax1 = plt.subplots(figsize=figsize)

    ax1.errorbar(x, y, xerr=xerr, yerr=yerr, fmt=fmt, markersize=markersize, capsize=capsize, label="Dati")

    fit_x = np.linspace(fit_range[0], fit_range[1], fit_points) if fit_range else np.linspace(min(x), max(x), fit_points)
    fit_y = func(fit_x, *result['parameters'])

    if fit_line:
        ax1.plot(fit_x, fit_y, color='red', linewidth=linewidth, label='Fit')
        if confidence_intervals:
            def f_wrapped(p):
                return func(fit_x, *p)
    
            _ , cov_y = propagate(f_wrapped, result['parameters'], result['covariance'])
            sigma_y = np.sqrt(np.diag(cov_y))

            ax1.fill_between(
                fit_x,
                fit_y - error_band * sigma_y,
                fit_y + error_band * sigma_y,
                color='red',
                alpha=0.2,
                label=f"Banda {error_band}σ"
            )
        if prediction_band:
            def f_wrapped(p):
                return func(fit_x, *p)
            
            _ , cov_y = propagate(f_wrapped, result['parameters'], result['covariance'])
            sigma_fit = np.sqrt(np.diag(cov_y))
            sigma_residuo = np.sqrt(np.sum(result['residuals']**2)/ (len(x) - len(result['parameters'])))
            sigma_pred = np.sqrt(sigma_fit**2 + sigma_residuo**2)
            

            ax1.fill_between(
                fit_x,
                fit_y - error_band * sigma_pred,
                fit_y + error_band * sigma_pred,
                color='green',
                alpha=0.2,
                label=f'Prediction band ({error_band}σ)'
            )


    # Calcolo gradi di libertà e p-value
    ndof = len(x) - len(result['parameters'])
    chi2_red = result['chi2'] / ndof if result['chi2'] is not None else None
    p_val = p_value(result['chi2'], ndof) if result['chi2'] is not None else None

    # Costruzione legenda con i parametri
    fit_info = []
    if show_chi_squared and result['chi2'] is not None:
        fit_info.append(f"$\\chi^2$/$n_\\mathrm{{dof}}$ = {result['chi2']:.1f} / {ndof} = {chi2_red:.2f}")
        fit_info.append(f"$p$-value = {p_val:.3f}")
    
    if show_fit_params:
        for i, val in enumerate(result['parameters']):
            err = np.sqrt(result['covariance'][i][i]) if result['covariance'] is not None else 0
            name = parameter_names[i] if parameter_names else f"p{i}"
            if abs(val) < 1e-3 or abs(val) >= 1e+3:
                fit_info.append(f"{name} = ${val:.3e} \\pm {err:.3e}$")
            else:
                fit_info.append(f"{name} = ${val:.3f} \\pm {err:.3f}$")
            

    ax1.legend(title="\n".join(fit_info), frameon=False, loc=legend_loc)

    ax1.set_xlabel(xlabel, fontsize=axis_fontsize)
    ax1.set_ylabel(ylabel, fontsize=axis_fontsize)
    ax1.set_title(title, fontsize=title_fontsize)
    ax1.grid(grid)

    if xlim:
        ax1.set_xlim(xlim)
    if ylim:
        ax1.set_ylim(ylim)

    if residuals:
        residual = result['residuals']
        fig_res, ax2 = plt.subplots(figsize=(figsize[0], figsize[1] / 2))
        ax2.errorbar(x, residual, xerr=xerr, yerr=yerr, fmt=fmt, color='black', markersize=markersize, capsize=capsize)
        ax2.axhline(0, color='gray', linestyle='--')
        ax2.set_xlabel(xlabel, fontsize=axis_fontsize)
        ax2.set_ylabel("Residui", fontsize=axis_fontsize)
        ax2.grid(grid)
        if save_path:
            fig_res.savefig(f"{save_path}_residui.png", dpi=dpi)
        

    if save_path:
        fig.savefig(f"{save_path}.png", dpi=dpi)

    







#------------integrali----------------
# Calcola integrale [Scipy]
def integral_scipy(f, a, b) : 
  integral = quad(f, a,b)
  return integral[0], integral[1]


#----------media pesata----------------
def media_pesata(x, sigma) :
    x = np.array(x)
    sigma = np.array(sigma)
    
    pesi = 1 / sigma**2
    m = np.sum(x * pesi) / np.sum(pesi)
    sigma_m = 1 / np.sqrt(np.sum(pesi))
    return m, sigma_m



#------------PDF's & CDF's----------------
def Gaussian(x, mu = 0, sigma = 1) :
	return (1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-((x - mu)**2) / (2 * sigma**2))

def Gaussian_standard(z):
  return (1/np.sqrt(2*np.pi))*np.exp((-z**2)/2) 

def Gaussian_cdf_ext(bin_edges, s, mu, sigma) :
  return s*norm.cdf(bin_edges, mu, sigma)

def Gaussian_cdf(bin_edges, mu, sigma) :
  return norm.cdf(bin_edges, mu, sigma)


#----------------HYPOTHESIS TESTING----------------

# Calcolo del p-value
def p_value(chi_square, ndof) :
  s = 1-chi2.cdf(chi_square, ndof)
  r = s*100
  return r

# z_test double sided
def z_test1(x1,x2,s1,s2) : 
  z = np.absolute(x1-x2)/np.sqrt(s1**2+s2**2)  #t di confronto
  R = quad(Gaussian_standard,-z,z) #calcolo del rapporto con l'integrale
  p_value = (1 - R[0])
  return p_value

# z test di ipotesi con un valore calcolato
def z_test2(x1,X,s) :  
  z = np.absolute(x1-X)/s  #t di confronto
  R = quad(Gaussian_standard,-z,z) #calcolo del rapporto con l'integrale
  p_value = (1 - R[0])
  return p_value

# t test con 1 vincolo
def t_test1(x1, X, N, err_media) : # N = N_dati 
	t_test = np.absolute(x1-X)/err_media
	R = t.cdf(-t_test, df=N-1)
	p_value = R*2
	return p_value
	
# t test con 2 vincoli
def t_test2(x1, x2, N, err1, err2) : 
	t_test = np.absolute(x1-x2)/np.sqrt(err1**2+err2**2)
	R = t.cdf(-t_test, df = N-1)
	p_value = R*2
	return p_value


