import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


from cherab.openadas import OpenADAS
from cherab.core.atomic import carbon, hydrogen
from scipy.optimize import lsq_linear


plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
#-------------------------------------------------------------------------------
# Definir funções (taxa de recomb, taxa de ionização, balanço-ionico)
#-------------------------------------------------------------------------------

def get_rates_recombination(element):

    coef_rec = {}
    for i in np.arange(1, element.atomic_number + 1):
        coef_rec[i] = adas.recombination_rate(element,int(i))

    return coef_rec

def get_rates_ionisation(element):
    coef_ionis = {}
    for i in np.arange(0, element.atomic_number):
        coef_ionis[i] = adas.ionisation_rate(element, int(i))

    return coef_ionis

def get_cx_rates(element):
    coef_ions = {}
    for i in np.arange(1, element.atomic_number + 1):
        coef_ions[i] = adas.thermal_cx_rate(hydrogen, 0,
                             elemento, int(i))
    return coef_ions

def temperature_profile(pho,nu,t0 = None):
    if t0 is None:
        return print("You need your temperature on reference system.")
    return t0*(1 - pho)**nu

def solve_ion_balance(element, n_e, t_e, coef_ion, coef_recom, nh0=None, coef_tcx=None): #Peguei da documentação

    atomic_number = element.atomic_number

    # construct the fractional abundance matrix
    matbal = np.zeros((atomic_number + 1, atomic_number + 1))

    matbal[0, 0] -= coef_ion[0](n_e, t_e)
    matbal[0, 1] += coef_recom[1](n_e, t_e)
    matbal[-1, -1] -= coef_recom[atomic_number](n_e, t_e)
    matbal[-1, -2] += coef_ion[atomic_number - 1](n_e, t_e)

    if nh0 is not None:
        matbal[0, 1] += nh0 / n_e * coef_tcx[1](n_e, t_e)
        matbal[-1, -1] -= nh0 / n_e * coef_tcx[atomic_number](n_e, t_e)

    for i in range(1, atomic_number):
        matbal[i, i - 1] += coef_ion[i - 1](n_e, t_e)
        matbal[i, i] -= (coef_ion[i](n_e, t_e) + coef_recom[i](n_e, t_e))
        matbal[i, i + 1] += coef_recom[i + 1](n_e, t_e)
        if nh0 is not None:
            matbal[i, i] -= nh0 / n_e * coef_tcx[i](n_e, t_e)
            matbal[i, i + 1] += nh0 / n_e * coef_tcx[i + 1](n_e, t_e)

    # for some reason calculation of stage abundance seems to yield better results than calculation of fractional abun.
    matbal = matbal * ne  # multiply by ne to calculate abundance instead of fractional abundance

    # add sum constraints. Sum of all stages should be equal to electron density
    matbal = np.concatenate((matbal, np.ones((1, matbal.shape[1]))), axis=0)

    # construct RHS of the balance steady-state equation
    rhs = np.zeros((matbal.shape[0]))
    rhs[-1] = ne

    abundance = lsq_linear(matbal, rhs, bounds=(0, ne))["x"]

    # normalize to ne to get fractional abundance
    frac_abundance = abundance / ne

    return frac_abundance


#-------------------------------------------------------------------------------
# Define a base de dados atômica
#-------------------------------------------------------------------------------
adas = OpenADAS(permit_extrapolation=True)

#-------------------------------------------------------------------------------
# parametros
#-------------------------------------------------------------------------------


elemento = carbon #elemento desejado
ne = 2.05e19 # m^-3
n0 = 1e15
temperature_steps = 100
a_coef = get_rates_recombination(elemento) # Coeficientes de recombinação para C I - VI
S_coef = get_rates_ionisation(elemento) # Coeficientes de ionização para C I - VI
cx_pec = get_cx_rates(elemento) # Coeficientes de emissão via CX para CVIII
te = [10 ** x for x in np.linspace(np.log10(S_coef[1].raw_data["te"].min()), np.log10(S_coef[1].raw_data["te"].max()),num=temperature_steps)] # temperaturas eletrônicas diversas [eV]
t0 = 400.0 # Temperatura eletrônica inicial no TCABR [eV] (SEVERO 2003)
a = 61
r = np.linspace(0, a, temperature_steps)


#-------------------------------------------------------------------------------
# Linhas importantes do carbono
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Criando matriz de abundancias parciais
#-------------------------------------------------------------------------------

partial_fract = np.zeros((len(te), elemento.atomic_number + 1))
pec_carbon_cx = np.zeros((len(te), elemento.atomic_number + 1))

#-------------------------------------------------------------------------------
# Cálculo do balanço iônico (Seguindo Severo 2004)
#-------------------------------------------------------------------------------

for idx, t_val in enumerate(te):
    partial_fract[idx, :] = solve_ion_balance(elemento, ne,t_val, S_coef, a_coef)
    for z in range(1, elemento.atomic_number + 1):
        pec_carbon_cx[idx, z] = cx_pec[z](ne, t_val)

#-------------------------------------------------------------------------------
# Calculo do perfil de simulado primário (Será trocado no futuro) 
#-------------------------------------------------------------------------------

pho = r/a # raio normalizado [0]
nu = 5/3 # Indice de perfilização

Temperatura = temperature_profile(pho,nu,t0) # Perfil toroidal tirado do Wesson - Apenas para testes

partial_fract_radius = np.zeros((len(Temperatura), elemento.atomic_number + 1)) # 

for idx, t_val in enumerate(Temperatura):
   partial_fract_radius[idx,:] = solve_ion_balance(elemento, ne, t_val, S_coef, a_coef, nh0=n0, coef_tcx=cx_pec)



#-------------------------------------------------------------------------------
# Plotando os resultados das abundâncias parciais + resultado dos coeficientes de CX
#-------------------------------------------------------------------------------

labels = ['CI', 'CII', 'CIII', 'CIV', 'CV', 'CVI', 'CVII']
plt.ion()

plt.figure(figsize=(9, 6))
for z in range(elemento.atomic_number + 1):
    plt.semilogx(te, partial_fract[:, z], label=labels[z], linewidth=2)

plt.xscale('log')
plt.xlabel(r'\(\\rho\) [eV]', fontsize=12)
plt.ylabel(r'Abundância Parcial \(f_Z\)', fontsize=12)
plt.title(f'Abundância Parcial de Íons de Carbono (ADAS) - \(n_e = {ne:.2e}\) \(m^{{-3}}\)', fontsize=13)
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.legend(fontsize=10)
plt.ylim(0, 1.05)
plt.xlim(1, 1000)

plt.figure(figsize=(9, 6))
for z in range(1,elemento.atomic_number +1 ):
    plt.loglog(te, pec_carbon_cx[:, z], label=labels[z], linewidth=2)
plt.xlabel(r'\(T_e\) [eV]', fontsize=12)
plt.ylabel(r'Coeficientes de Recombinação [fótons \(\mathrm{m^{-3} s^{-1}}\)]', fontsize=12)
plt.title('Emissividade por Troca de Carga', fontsize=13)
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.legend(fontsize=10)
plt.tight_layout()


plt.figure(figsize=(9,6))
plt.plot(pho,Temperatura, c = "r", label = f"\(\\nu = {nu:.2f}\)")
plt.xlabel(f"\(\\rho\) (r/a)", fontsize = 12)
plt.ylabel(r"\(T_e\) [eV]", fontsize = 12)
plt.title(f'Perfil de temperatura toroidal',fontsize=13)
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.legend(fontsize = 10)
plt.tight_layout()

plt.figure(figsize=(9,6))
for z in range(1, elemento.atomic_number+1):
    plt.plot(pho,partial_fract_radius[:, z], label = labels[z])
plt.xlabel(f"\(\\rho\) (r/a)", fontsize = 12)
plt.ylabel(r"Frações [eV]", fontsize = 12)
plt.title(f'Perfil das frações parciais',fontsize=13)
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.legend(fontsize = 10)
plt.tight_layout()


plt.ioff()

plt.show()
