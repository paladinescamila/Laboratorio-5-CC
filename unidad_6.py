# UNIDAD 6: ECUACIONES DIFERENCIALES ORDINARIAS

import time
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt

t, y = sym.symbols("t y")


# Método de Euler
def euler(ODE, t0, y0, tf, h):
    """
    Entrada: una Ecuación Diferencial Ordinaria de primer orden ODE(t, y), dos reales
            t0 & y0 que representan las condiciones iniciales y(t0) = y0, un real
            tf que es el límite superior del intervalo de análisis, y un real h que 
            representa el incremento de t.
    Salida: los puntos (tk, yk) que corresponden a la solución de la ODE, para
            t0 <= tk <= tf.
    """

    f = sym.lambdify([t, y], ODE)
    pasos = [(t0, y0)]
    tk, yk = pasos[-1]

    while (tk < tf):
        yk += f(tk, yk)*h
        tk += h
        pasos.append((tk, yk))
    
    return pasos


# Método de Serie de Taylor
def taylor(ODE, t0, y0, tf, h):
    """
    Entrada: una Ecuación Diferencial Ordinaria de primer orden ODE(t, y), dos reales
            t0 & y0 que representan las condiciones iniciales y(t0) = y0, un real
            tf que es el límite superior del intervalo de análisis, y un real h que 
            representa el incremento de t.
    Salida: los puntos (tk, yk) que corresponden a la solución de la ODE, para
            t0 <= tk <= tf.
    """

    f = sym.lambdify([t, y], ODE)
    pasos = [(t0, y0)]
    tk, yk = pasos[-1]

    ypp = sym.diff(ODE, t) + sym.diff(ODE, y)*ODE
    fpp = sym.lambdify([t, y], ypp)

    while (tk < tf):
        yk += f(tk, yk)*h + (fpp(tk, yk)/2)*(h**2)
        tk += h
        pasos.append((tk, yk))
    
    return pasos


# Método de Runge-Kutta (orden 2)
def runge_kutta_2(ODE, t0, y0, tf, h):
    """
    Entrada: una Ecuación Diferencial Ordinaria de primer orden ODE(t, y), dos reales
            t0 & y0 que representan las condiciones iniciales y(t0) = y0, un real
            tf que es el límite superior del intervalo de análisis, y un real h que 
            representa el incremento de t.
    Salida: los puntos (tk, yk) que corresponden a la solución de la ODE, para
            t0 <= tk <= tf.
    """

    f = sym.lambdify([t, y], ODE)
    pasos = [(t0, y0)]
    tk, yk = pasos[-1]

    while (tk < tf):
        k1 = f(tk, yk)*h
        k2 = f(tk + h, yk + k1)*h
        yk += (1/2) * (k1 + k2)
        tk += h
        pasos.append((tk, yk))
    
    return pasos


# Método de Runge-Kutta (orden 4)
def runge_kutta_4(ODE, t0, y0, tf, h):
    """
    Entrada: una Ecuación Diferencial Ordinaria de primer orden ODE(t, y), dos reales
            t0 & y0 que representan las condiciones iniciales y(t0) = y0, un real
            tf que es el límite superior del intervalo de análisis, y un real h que 
            representa el incremento de t.
    Salida: los puntos (tk, yk) que corresponden a la solución de la ODE, para
            t0 <= tk <= tf.
    """

    f = sym.lambdify([t, y], ODE)
    pasos = [(t0, y0)]
    tk, yk = pasos[-1]

    while (tk < tf):
        k1 = f(tk, yk)*h
        k2 = f(tk + h/2, yk + k1/2)*h
        k3 = f(tk + h/2, yk + k2/2)*h
        k4 = f(tk + h, yk + k3)*h
        yk += (1/6) * (k1 + 2*k2 + 2*k3 + k4)
        tk += h
        pasos.append((tk, yk))
    
    return pasos


# Método Multipaso (2 pasos)
def multipaso_2(ODE, t0, y0, tf, h):
    """
    Entrada: una Ecuación Diferencial Ordinaria de primer orden ODE(t, y), dos reales
            t0 & y0 que representan las condiciones iniciales y(t0) = y0, un real
            tf que es el límite superior del intervalo de análisis, y un real h que 
            representa el incremento de t.
    Salida: los puntos (tk, yk) que corresponden a la solución de la ODE, para
            t0 <= tk <= tf.
    """

    f = sym.lambdify([t, y], ODE)
    t1, y1 = runge_kutta_4(ODE, t0, y0, t0+h, h)[1]
    pasos = [(t0, y0), (t1, y1)]
    tk, yk = pasos[-1]

    while (tk < tf):
        ypk_1 = f(pasos[-2][0], pasos[-2][1])
        yk += (1/2) * (3*f(tk, yk) - ypk_1) * h
        tk += h
        pasos.append((tk, yk))
    
    return pasos


# Método Multipaso (4 pasos)
def multipaso_4(ODE, t0, y0, tf, h):
    """
    Entrada: una Ecuación Diferencial Ordinaria de primer orden ODE(t, y), dos reales
            t0 & y0 que representan las condiciones iniciales y(t0) = y0, un real
            tf que es el límite superior del intervalo de análisis, y un real h que 
            representa el incremento de t.
    Salida: los puntos (tk, yk) que corresponden a la solución de la ODE, para
            t0 <= tk <= tf.
    """

    f = sym.lambdify([t, y], ODE)
    [(t1, y1), (t2, y2), (t3, y3)] = runge_kutta_4(ODE, t0, y0, t0+3*h, h)[1:4]
    pasos = [(t0, y0), (t1, y1), (t2, y2), (t3, y3)]
    tk, yk = pasos[-1]

    while (tk < tf):
        ypk_1 = f(pasos[-2][0], pasos[-2][1])
        ypk_2 = f(pasos[-3][0], pasos[-3][1])
        ypk_3 = f(pasos[-4][0], pasos[-4][1])
        yk += (1/24) * (55*f(tk, yk) - 59*ypk_1 + 37*ypk_2 - 9*ypk_3) * h
        tk += h
        pasos.append((tk, yk))
    
    return pasos


# Resolución de Ecuaciones de Orden Superior
def ODEs_superior(ODE, t0, y0s, tf, h, n):
    """
    Entrada: una Ecuación Diferencial Ordinaria ODE de orden n, un real t0 y una lista
            de reales y0s que representan las condiciones iniciales de cada orden, un 
            real tf que es el límite superior del intervalo de análisis, un real h
            que representa el incremento de t, y un entero n que es el orden de la ODE.
    Salida: los puntos (tk, yk) que corresponden a la solución de la ODE, para
            t0 <= tk <= tf.
    """

    pasos = [(t0, y0s[0])]
    tk, yk = t0, list (y0s)

    while (tk < tf):
        f = sym.lambdify([t, y], ODE)
        for i in range(n-1, -1, -1):
            yk[i] += f(tk, yk[i])*h
            f = sym.lambdify([t, y], yk[i])
        tk += h
        pasos.append((tk, yk[0]))
        
    return pasos


# Método de Diferencias Finitas
def diferencias_finitas(ODE, t0, y0, tf, yf, n): # ¿VARIAMOS N O H?
    """
    """

    ti = np.linspace(t0, tf, n)
    f = sym.lambdify(t, ODE) # ¿DEPENDE SÓLO DE t? (ES OPCIONAL SI SE QUIERE METER UNA Y')
    h = (tf - t0) / (n - 1)

    A = [[0 for _ in range(n - 2)] for _ in range(n - 2)]
    b = [f(ti[i])*h**2 for i in range(1, n - 1)]
    
    for i in range(1, n - 1):

        j = i - 1
        A[j][j] = -2

        if (i == 1): b[j] -= y0
        if (i == n - 2): b[j] -= yf
        if (i > 1): A[j][j-1] = 1
        if (i < n-2): A[j][j+1] = 1
    
    yi = [y0] + list (np.linalg.solve(A, b)) + [yf]
    pasos = [(ti[i], yi[i]) for i in range(n)] # ¿SERÁ NECESARIO HACERLO AQUÍ?
    
    return pasos


# Método de Elementos Finitos
def elementos_finitos(ODE, t0, y0, tf, yf, n):
    """
    """

    ti = np.linspace(t0, tf, n)
    f = sym.lambdify(t, ODE)
    A, b = [], []
    coheficientes_ypp = [i for i in range(n)]

    A.append([t0**i for i in range(n)])
    # PARTE DONDE SE AGREGAN LOS DEMÁS COHEFICIENTES
    A.append([tf**i for i in range(n)])

    b.append(y0)
    # PARTE DONDE SE AGREGAN LOS DEMÁS VALORES DE B
    b.append(yf)


ODE = 6*t # EN EL LABORATORIO NO SE PUEDEN PONER POLINOMIOS
analitica = t**3
t0, y0, tf, yf = 0, 0, 1, 1
elementos_finitos(ODE, t0, y0, tf, yf, 5)