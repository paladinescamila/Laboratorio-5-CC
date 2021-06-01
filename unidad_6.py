# UNIDAD 6: ECUACIONES DIFERENCIALES ORDINARIAS

import time
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt

t, y = sym.symbols("t y")


# Método de Euler
def euler(ODE, t0, y0, h, tn):
    """
    Entrada: una Ecuación Diferencial Ordinaria de primer orden ODE(t, y), dos reales
            t0 & y0 que representan las condiciones iniciales y(t0) = y0, un real
            h que representa el incremento de t, y un real tn que es el límite
            superior del intervalo de análisis.
    Salida: los puntos (tk, yk) que corresponden a la solución de la ODE, para
            t0 <= tk <= tn.
    """

    f = sym.lambdify([t, y], ODE)
    pasos = [(t0, y0)]
    tk, yk = pasos[-1]

    while (tk < tn):
        yk += f(tk, yk)*h
        tk += h
        pasos.append((tk, yk))
    
    return pasos


# Método de Serie de Taylor
def taylor(ODE, t0, y0, h, tn):
    """
    Entrada: una Ecuación Diferencial Ordinaria de primer orden ODE(t, y), dos reales
            t0 & y0 que representan las condiciones iniciales y(t0) = y0, un real
            h que representa el incremento de t, y un real tn que es el límite
            superior del intervalo de análisis.
    Salida: los puntos (tk, yk) que corresponden a la solución de la ODE, para
            t0 <= tk <= tn.
    """

    f = sym.lambdify([t, y], ODE)
    pasos = [(t0, y0)]
    tk, yk = pasos[-1]

    ypp = sym.diff(ODE, t) + sym.diff(ODE, y)*ODE
    fpp = sym.lambdify([t, y], ypp)

    while (tk < tn):
        yk += f(tk, yk)*h + (fpp(tk, yk)/2)*(h**2)
        tk += h
        pasos.append((tk, yk))
    
    return pasos


# Método de Runge-Kutta (orden 2)
def runge_kutta_2(ODE, t0, y0, h, tn):
    """
    Entrada: una Ecuación Diferencial Ordinaria de primer orden ODE(t, y), dos reales
            t0 & y0 que representan las condiciones iniciales y(t0) = y0, un real
            h que representa el incremento de t, y un real tn que es el límite
            superior del intervalo de análisis.
    Salida: los puntos (tk, yk) que corresponden a la solución de la ODE, para
            t0 <= tk <= tn.
    """

    f = sym.lambdify([t, y], ODE)
    pasos = [(t0, y0)]
    tk, yk = pasos[-1]

    while (tk < tn):
        k1 = f(tk, yk)*h
        k2 = f(tk + h, yk + k1)*h
        yk += (1/2) * (k1 + k2)
        tk += h
        pasos.append((tk, yk))
    
    return pasos


# Método de Runge-Kutta (orden 4)
def runge_kutta_4(ODE, t0, y0, h, tn):
    """
    Entrada: una Ecuación Diferencial Ordinaria de primer orden ODE(t, y), dos reales
            t0 & y0 que representan las condiciones iniciales y(t0) = y0, un real
            h que representa el incremento de t, y un real tn que es el límite
            superior del intervalo de análisis.
    Salida: los puntos (tk, yk) que corresponden a la solución de la ODE, para
            t0 <= tk <= tn.
    """

    f = sym.lambdify([t, y], ODE)
    pasos = [(t0, y0)]
    tk, yk = pasos[-1]

    while (tk < tn):
        k1 = f(tk, yk)*h
        k2 = f(tk + h/2, yk + k1/2)*h
        k3 = f(tk + h/2, yk + k2/2)*h
        k4 = f(tk + h, yk + k3)*h
        yk += (1/6) * (k1 + 2*k2 + 2*k3 + k4)
        tk += h
        pasos.append((tk, yk))
    
    return pasos


# Método Multipaso (2 pasos)
def multipaso_2(ODE, t0, y0, h, tn):
    """
    Entrada: una Ecuación Diferencial Ordinaria de primer orden ODE(t, y), dos reales
            t0 & y0 que representan las condiciones iniciales y(t0) = y0, un real
            h que representa el incremento de t, y un real tn que es el límite
            superior del intervalo de análisis.
    Salida: los puntos (tk, yk) que corresponden a la solución de la ODE, para
            t0 <= tk <= tn.
    """

    f = sym.lambdify([t, y], ODE)
    t1, y1 = runge_kutta_4(ODE, t0, y0, h, t0+h)[1]
    pasos = [(t0, y0), (t1, y1)]
    tk, yk = pasos[-1]

    while (tk < tn):
        ypk_1 = f(pasos[-2][0], pasos[-2][1])
        yk += (1/2) * (3*f(tk, yk) - ypk_1) * h
        tk += h
        pasos.append((tk, yk))
    
    return pasos


# Método Multipaso (4 pasos)
def multipaso_4(ODE, t0, y0, h, tn):
    """
    Entrada: una Ecuación Diferencial Ordinaria de primer orden ODE(t, y), dos reales
            t0 & y0 que representan las condiciones iniciales y(t0) = y0, un real
            h que representa el incremento de t, y un real tn que es el límite
            superior del intervalo de análisis.
    Salida: los puntos (tk, yk) que corresponden a la solución de la ODE, para
            t0 <= tk <= tn.
    """

    f = sym.lambdify([t, y], ODE)
    [(t1, y1), (t2, y2), (t3, y3)] = runge_kutta_4(ODE, t0, y0, h, t0+3*h)[1:4]
    pasos = [(t0, y0), (t1, y1), (t2, y2), (t3, y3)]
    tk, yk = pasos[-1]

    while (tk < tn):
        ypk_1 = f(pasos[-2][0], pasos[-2][1])
        ypk_2 = f(pasos[-3][0], pasos[-3][1])
        ypk_3 = f(pasos[-4][0], pasos[-4][1])
        yk += (1/24) * (55*f(tk, yk) - 59*ypk_1 + 37*ypk_2 - 9*ypk_3) * h
        tk += h
        pasos.append((tk, yk))
    
    return pasos


# Resolución de Ecuaciones de Orden Superior
def ODEs_superior(ODE, t0, y0s, h, tn, orden):
    """
    """

    pasos = [(t0, y0s[0])]
    tk, yk = t0, list (y0s)

    while (tk < tn):
        f = sym.lambdify([t, y], ODE)
        for i in range(orden-1, -1, -1):
            yk[i] += f(tk, yk[i])*h
            f = sym.lambdify([t, y], yk[i])
        tk += h
        pasos.append((tk, yk[0]))
        
    return pasos