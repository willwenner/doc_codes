import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re 

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
                        densit.append(float(num))
        
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
        
    ### Combinar a temperatura e densidade para cada PEC.
    
    ### Retorna um dicionário com matrizes de dimensão (24x21) - (N_densidades x N_temperaturas)
    return dic_excit, dic_recom
