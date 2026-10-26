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
    base = 1 - (x_i/a)**2

    resultado = DeltaV * gama * (-2 * (x_i / a**2)) * (base**power)
    
    if power < 0:
        invalido = (base<=0)
    else:
        invalido = base<0
    
    return np.where(invalido,0.0,resultado)

def Vapp(eps, Vphi, W_ij):
    """Aplicação do modelo direto propagando a velocidade de plasma.

    Parâmetros:
        x_i (float): Valores de emissividade ao longo do raio.
        Vphi (float): Valores de velocidade toroidal ao longo do raio.
        W_ij (arr): Matriz de peso calculada utilizando o método de Onion Peeling.
    Retorna: 
        Vapp (arr): Curvas de valores de velocidade aparente.    
    """
    U = W_ij @ (eps*Vphi)
    S = W_ij @ eps
    return U/S