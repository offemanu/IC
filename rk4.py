import numpy as np

def metodo_rk4(f, y0, t):
    n = len(t) # Número de pontos
    y = np.zeros(n) # Vetor para armazenar a solução
    y[0] = y0
    for i in range(n - 1):
        h = t[i + 1] - t[i] # Tamanho do passo 
        ti = t[i] # Tempo atual
        yi = y[i] # Solução atual
        # Estágios:
        k1 = f(ti, yi)
        k2 = f(ti + h / 2, yi + (h/2)*k1)
        k3 = f(ti + h / 2, yi + (h/2)*k2)
        k4 = f(ti + h, yi + h*k3)
        y[i + 1] = yi + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
    return y # Vetor 1D contendo a solução aproximada nos tempos especificados

def metodo_rk4_vector(f, y0, t):
    n = len(t) # Número de pontos
    y = np.zeros((n, len(y0)))  # Vetor 2D: n linhas (tempos) × m colunas (variáveis do sistema)
    y[0] = y0 # Vetor de condições iniciais [y1_0, y2_0, ..., yn_0] no tempo t[0]
    
    for i in range(n - 1):
        h = t[i+1] - t[i] # Tamanho do passo 
        yi = np.asarray(y[i])
        # Estágios do método RK4 - VERSÃO VETORIAL:
        k1 = f(t[i], yi)
        k2 = f(t[i] + h/2, yi + k1/2)
        k3 = f(t[i] + h/2, yi + k2/2)
        k4 = f(t[i] + h, yi + k3)
        
        y[i+1] = yi + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
        
    return y # Vetor 2D onde:- Cada LINHA representa o estado completo do sistema em um tempo
#                            - Cada COLUNA representa a evolução temporal de uma variável
#                            - Formato: [[y1(t0), y2(t0), ...], [y1(t1), y2(t1), ...],...]