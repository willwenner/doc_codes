# PEC - Coeficientes de emissividade de photons.

Os coeficientes de emissividade de fótons (PEC) são utilizados para o cálculo da emissividade de uma transição atômica, eles são obtidos através da biblioteca 
[Open-ADAS](https://open.adas.ac.uk/) que armazena dados atômicos obtidos por modelos atômicos colisional-radioativos. 

## Versão 0.1:

- Leitura do arquivo tipo ADF15 que contÊm os valores de PEC para as transições do CARBONO.
- O programa lê e retorna dois dicionários: excit, recom. Os dicionários são escritos da seguinte forma:

    - Temperatura Eletrônica (eV)
    - Densidade Eletrônica ($cm^(-3)$)
    - PEC: Matrix N_densidades x N_temperaturas (24x21)

    
### BUGS:
- Problemas de valores repetidos (???) 
- Não lê a pec para CRX

### Arquivos utilizados:

Arquivo: "pec96#c_bnd#c5.dat"  
Data: 29 de setembro de 2026 - 16h - UTC-3 (BRT)  
URL: [OPEN-ADAS Coeficientes de emisividade de photons do Carbono](https://open.adas.ac.uk/adf15?element=Carbon&charge=&wave_min=&wave_max=&resolveby=transition&searching=1#searchbutton)  
Checksum: f4cb15270cde2935c520eb30aa5b737b1c6e1e23468171bf5d841122368bf16d

