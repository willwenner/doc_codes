import numpy as np

def onion_peeling_method(n_cordas, y_i):
    """
    Função implementa o método de Pearse, descrito por Dasch (1992), num conjunto de cordas y_i.

    Parâmetros:
    n_cordas (int): número de cordas a ser utilizado no experimento
    y_i (arr): malha criada com o número de cordas.

    Retorna:
    W_ij (arr): Matriz de pesos para o caso direto
    D_ij (arr): Matriz de deconvolução para a inversão através da inversão de W_ij
    """
    # Criação da malha 
    dr = y_i[1]-y_i[0]
    
    # Criação dos anéis utilizando o passo
    r_externo = y_i[:]+dr/2
    r_interno = y_i[:]-dr/2
    
    # Tratamento dos valores de borda
    r_externo[-1] = n_cordas*dr
    r_interno[0] = 0
    
    # Criando a matriz W_ij
    W_ij = np.zeros((n_cordas,n_cordas))
    
    # Preenchendo a matriz W_ij
    for i in range(0,n_cordas):
        for j in range(0,n_cordas):
            if r_externo[j] > abs(y_i[i]):
                W_ij[i,j] = 2*(np.sqrt(max(0,r_externo[j]**2-y_i[i]**2))-np.sqrt(max(0,r_interno[j]**2-y_i[i]**2)))
            else:
                W_ij[i,j] = 0
    
    D_ij = np.linalg.inv(W_ij)
    return W_ij, D_ij