import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

from cherab.openadas import OpenADAS
from cherab.core.atomic import carbon


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

def get_ion_balance(n_e,t_e, element, S_coef, a_coef):
    # Usando o resultado do balanço iônico 
    r = np.zeros(elemento.atomic_number + 1)
    r[0] = 1.0

    for z in range(elemento.atomic_number):
        S_val = S_coef[z](n_e, t_e)
        a_val = a_coef[z + 1](n_e, t_e)
        r[z + 1] = r[z] * (S_val/a_val)

    partial_fraction =  r / np.sum(r)

    return partial_fraction

#-------------------------------------------------------------------------------
# Define a base de dados atômica
#-------------------------------------------------------------------------------
adas = OpenADAS()

#-------------------------------------------------------------------------------
# parametros
#-------------------------------------------------------------------------------

elemento = carbon #elemento desejado
ne = 2.05e19 # m^-3
te = np.linspace(1,1000, int(1e4)) # eV
a_coef = get_rates_recombination(elemento) 
S_coef = get_rates_ionisation(elemento)


#-------------------------------------------------------------------------------
# Criando matriz de abundancias parciais
#-------------------------------------------------------------------------------

partial_fract = np.zeros((len(te), elemento.atomic_number + 1))

#-------------------------------------------------------------------------------
# Cálculo do balanço iônico (Seguindo Severo 2004)
#-------------------------------------------------------------------------------

for idx, t_val in enumerate(te):
    partial_fract[idx, :] = get_ion_balance(ne,t_val, elemento, S_coef, a_coef)


#-------------------------------------------------------------------------------
# Plotando os resultados das abundâncias parciais
#-------------------------------------------------------------------------------
plt.figure(figsize=(9, 6))

labels = ['CI', 'CII', 'CIII', 'CIV', 'CV', 'CVI', 'CVII']

for z in range(elemento.atomic_number + 1):
    plt.plot(te, partial_fract[:, z], label=labels[z], linewidth=2)

plt.xscale('log')
plt.xlabel(r'$\(T_e\) [eV]$', fontsize=12)
plt.ylabel(r'Abundância Parcial $\(f_Z\)$', fontsize=12)
plt.title(f'Abundância Parcial de Íons de Carbono (ADAS) - $\(n_e = {ne:.2e}\)$ $\(m^{{-3}}\)$', fontsize=13)
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.legend(fontsize=10)
plt.ylim(0, 1.05)
plt.xlim(1, 1000)

plt.show()


