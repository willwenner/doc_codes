import numpy as np


def Vphi_model(x_i, a, DeltaV, gama):
    """Calcula o modelo de sino da velocidade toroidal eliminando os valores imaginários e divisões por zero

    Parâmetros:
        x_i (float): Intervalos de valores de raio 
        a (float): raio maior do tokamak 
        DeltaV (float): Amplitude da velocidade 
        Gama (float):  Expoente efetivo
    Retorna:
        V_phi (arr): V_phi(x_i) simulada seguindo modelo descrito em Wesson J., Tokamaks (2004).
    """
    # Recusa valores de a negativos e iguais a 0
    if a <= 0:
        return print("Please insert the original Radius of your Tokamak")

    # Adapta o código para listas e números
    x_i = np.asarray(x_i)

    # Tratamento de possíveis valores negativos
    power = gama - 1
    base = np.maximum(1e-18, 1 - (x_i/a)**2)

    resultado = DeltaV * gama * (-2 * (x_i / a**2)) * (base**power)
    
    return resultado


def onion_peeling_method(n_cordas, y_i):
    """
    Função implementa o método de Pearse, descrito por Dasch (1992), num conjunto de cordas y_i.

    Parâmetros:
    n_cordas (int): número de cordas a ser utilizado no experimento
    r_max (int, float): valor do raio máximo a ser utilizado 

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

def Vapp(x, eps, Vphi, W_ij):
    """Aplicação do modelo direto propagando a velocidade de plasma.

    Parâmetros:
        x (list, float): Valores de emissividade ao longo do raio.
        Vphi (float): Valores de velocidade toroidal ao longo do raio.
        W_ij (arr): Matriz de peso calculada utilizando o método de Onion Peeling.
    Retorna: 
        Vapp (arr): Curvas de valores de velocidade aparente.    
    """
    x = np.asarray(x)

    W_ij, D_ij = onion_peeling_method(len(x),x)
    
    U = W_ij @ (eps*Vphi)
    S = W_ij @ eps
    return U/S