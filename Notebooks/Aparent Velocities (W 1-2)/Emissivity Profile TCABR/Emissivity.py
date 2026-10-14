import pandas as pd
import numpy as np
import scipy as sp
from scipy.constants import h, c
import matplotlib.pyplot as plt
import re 


from cherab.core.atomic import carbon, hydrogen
from cherab.openadas import OpenADAS 

### Ler o documento
def read_file(path):
    with open(path, 'r') as archive:
        lines = archive.read().splitlines()
        lines.pop(0)
        return lines
        
### Encontrar as linhas de cada linha atômica
def new_specline(self):
    return "/TYPE" in self and "A" in self

def get_parameters(self):
    temp = []
    densit = []

    idx = 0
    length_data = len(self)
   
    while idx < length_data:
        row = self[idx]
        
        if new_specline(row):
            match = re.search(r"(\d+\.\d+)A",row)
            match_type = re.search(r'/TYPE\s=\s(\S+)',row)
            
            if match and match.group(1) == "33.8" and match_type.group(1) == "EXCIT":
                idx += 1
                
                for _ in range(3):
                    row_data = self[idx]
                    # Separa os valores tratadas pela notação científica do FORTRAN (D -> E):
                    for num in row_data.replace('D', 'E').split():
                        densit.append(float(num)*1e6)
        
                    idx += 1
            if match and match.group(1) == "33.8" and match_type.group(1) == "RECOM":  
                idx += 4
                
                for _ in range(3):
                    row_data = self[idx]
                    # Separa os valores tratadas pela notação científica do FORTRAN (D -> E):
                    for num in row_data.replace('D', 'E').split():
                        temp.append(float(num))
        
                    idx += 1
        idx += 1
    
    return temp, densit
    
### Gerar o dicionário com as PECS, TRANSIÇÕES e COMPRIMENTOS DE ONDA para cada comprimento de onda
def pec_matrix_creation(self):
    dic_excit = {} 
    dic_recom = {}
    temp, density = get_parameters(self)
    # Indice no qual iremos iterar
    idx = 0
    length_data = len(self)

    ### Iterar sobre todas as linhas
    while idx < length_data:
        
        row = self[idx]
        
        ### Checar se linha é comprimento de onda novo (True se sim e False se não)
        if new_specline(row): 
            match_type = re.search(r'/TYPE\s=\s(\S+)',row) # Encontra o tipo de transição
            
            ### Se o comprimento de onda estiver ali pula as linhas até atingir as pecs
            if match_type.group(1) == "EXCIT":
                match = re.search(r"(\d+\.\d+)A",row) ### Encontra o comprimento de onda
                wave_value = float(match.group(1)) # Pega o comprimento de onda em [A]
                
                ###Pular as lihas
                idx += 7
                
                ### Define as listas de valores e a lista que servirá como matriz
                values = []
                
                ### Itera sobre os blocos de 3 em 3 linhas para criar as linhas:
                for _ in range(int(3 * 24)):
                    if idx < length_data:
                        row_data = self[idx]
                        
                        # Separa os valores tratadas pela notação científica do FORTRAN (D -> E):
                        for num in row_data.replace('D', 'E').split():
                            values.append(float(num))
                
                        idx += 1
                        
                pec_for_wavelength = np.array(values).reshape(24, 21)   
                
                ### Adiciona as matrizes para cada comprimento de onda dentro da matriz pec
                dic_excit[wave_value] = {
                    "PEC": pec_for_wavelength,
                    "temp": temp,
                    "density": density
                }
                ### Continue apenas para não bugar o loop? (gepeto ajudou)
                continue
                
            if match_type.group(1) == "RECOM":
                match = re.search(r"(\d+\.\d+)A",row) ### Encontra o comprimento de onda
                wave_value = float(match.group(1)) # Pega o comprimento de onda em [A]
                
                ###Pular as lihas
                idx += 7
                
                ### Define as listas de valores e a lista que servirá como matriz
                values = []
                
                ### Itera sobre os blocos de 3 em 3 linhas para criar as linhas:
                for _ in range(int(3 * 24)):
                    if idx < length_data:
                        row_data = self[idx]
                        
                        # Separa os valores tratadas pela notação científica do FORTRAN (D -> E):
                        for num in row_data.replace('D', 'E').split():
                            values.append(float(num))
                
                        idx += 1
                        
                pec_for_wavelength = np.array(values).reshape(24, 21)   
                
                ### Adiciona as matrizes para cada comprimento de onda dentro da matriz pec
                dic_recom[wave_value] = {
                    "PEC": pec_for_wavelength,
                    "temp": temp,
                    "density": density
                }
                                             
                ### Continue apenas para não bugar o loop? (gepeto ajudou)
                continue
                
        ### Dicionário contendo os comprimentos de onda, o tipo de transição e as matrizes correspondentes)                       
        idx += 1
    
    ### Retorna um dicionário com matrizes de dimensão (24x21) - (N_densidades x N_temperaturas)
    return dic_excit, dic_recom


def get_pec(excit, recom,  wavelength, ne, te, kind = None, n0 = None, permit_extrapolation = False):
  
    # Pegando as listas de temperaturas e densidades
    temp_e = np.log(excit[wavelength]["temp"])
    dens_e = np.log(excit[wavelength]["density"])

    ## Calculando  
    if (min(dens_e) <= np.log(ne) <= max(dens_e)) and (min(temp_e)<= np.log(te) <= max(temp_e)):
        if kind == "EXCIT":
            pec_discret = excit[wavelength]["PEC"]
            interp_pec = sp.interpolate.RegularGridInterpolator((dens_e, temp_e), pec_discret, method = "linear", bounds_error = True, fill_value = permit_extrapolation)
            return interp_pec((np.log(ne), np.log(te)))
        if kind == "RECOM":
            pec_discret = recom[wavelength]["PEC"]
            interp_pec = sp.interpolate.RegularGridInterpolator((dens_e, temp_e), pec_discret, method = "linear", bounds_error = True, fill_value = permit_extrapolation)
            return interp_pec((np.log(ne), np.log(te)))
        else:
            return print("Please specify the kind of transition: 'EXCIT' or 'RECOM'")
    else:
        return print("The density and temperature values are out of the range of the data. Please provide values within the range of the data.")

def get_emissivity(ne, te, n_c5_eff, n_c6_eff, n0=None, energy = None, Zeff = None, B0 = None, ti = None):
    
    pec_excit_adas = adas.impact_excitation_pec(carbon, 5, (8, 7)) # Wm^3
    pec_recomb_adas = adas.recombination_pec(carbon, 5, (8, 7)) # Wm^3
    
    wavelength = 529.97e-9  # meters
    conv_factor = (h * c) / wavelength # SI 
    
    emissivity_excit = n_c5_eff * (pec_excit_adas(ne, te)/conv_factor) # Wm^3
    emissivity_recomb = n_c6_eff * (pec_recomb_adas(ne, te)/conv_factor)  #Wm^3
    
    total_rate = emissivity_excit + emissivity_recomb
    
    if n0 is not None:
        pec_cx_adas = adas.beam_cx_pec(hydrogen, carbon, 6, (8, 7))
        cxr1, cxr2 = pec_cx_adas
        emissivity_cx = (n0 / ne) * n_c6_eff * cxr1(energy, ti, n0, 1, B0)
        total_rate += emissivity_cx/conv_factor

    return ne * total_rate
    
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
        coef_ions[i] = adas.thermal_cx_pec(hydrogen, 0,
                             elemento, int(i))
    return coef_ions


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
    matbal = matbal * n_e  # multiply by ne to calculate abundance instead of fractional abundance

    # add sum constraints. Sum of all stages should be equal to electron density
    matbal = np.concatenate((matbal, np.ones((1, matbal.shape[1]))), axis=0)

    # construct RHS of the balance steady-state equation
    rhs = np.zeros((matbal.shape[0]))
    rhs[-1] = n_e

    abundance = sp.optimize.lsq_linear(matbal, rhs, bounds=(0, n_e))["x"]

    # normalize to ne to get fractional abundance
    frac_abundance = abundance / n_e

    return frac_abundance
