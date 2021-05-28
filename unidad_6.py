# UNIDAD 6: ECUACIONES DIFERENCIALES ORDINARIAS

import time
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt

t = sym.Symbol("t")
y = sym.Symbol("y")


# Método de Euler
def euler(ODE, h, t0, y0, a, b):
    """
    """

    f = sym.lambdify([t, y], ODE)
    tk, yk, pasos = t0, y0, [(t0, y0)]

    while (tk <= b):
        yk += f(tk, yk)*h
        tk += h
        pasos += [(tk, yk)]
    
    return pasos


# Método de Serie de Taylor
def taylor(ODE, h, t0, y0, a, b):
    """
    """

    f = sym.lambdify([t, y], ODE)
    tk, yk, pasos = t0, y0, [(t0, y0)]
    ypp = sym.diff(ODE, t) + sym.diff(ODE, y)*ODE
    fpp = sym.lambdify([t, y], ypp)

    while (tk <= b):
        yk += f(tk, yk)*h + (fpp(tk, yk)/2)*h**2
        tk += h
        pasos += [(tk, yk)]
    
    return pasos


# Método de Runge-Kutta de orden 2
def runge_kutta_2(ODE, h, t0, y0, a, b):
    """
    """

    f = sym.lambdify([t, y], ODE)
    tk, yk, pasos = t0, y0, [(t0, y0)]

    while (tk <= b):
        k1 = f(tk, yk)*h
        k2 = f(tk + h, yk + k1)*h
        yk += (1/2) * (k1 + k2)
        tk += h
        pasos += [(tk, yk)]
    
    return pasos


# Método de Runge-Kutta de orden 4
def runge_kutta_4(ODE, h, t0, y0, a, b):
    """
    """

    f = sym.lambdify([t, y], ODE)
    tk, yk, pasos = t0, y0, [(t0, y0)]

    while (tk <= b):
        k1 = f(tk, yk)*h
        k2 = f(tk + h/2, yk + k1/2)*h
        k3 = f(tk + h/2, yk + k2/2)*h
        k4 = f(tk + h, yk + k3)*h
        yk += (1/6) * (k1 + 2*k2 + 2*k3 + k4)
        tk += h
        pasos += [(tk, yk)]
    
    return pasos


# Método Multipaso de 2 pasos
def multipaso_2(ODE, h, t0, y0, a, b):
    """
    """

    f = sym.lambdify([t, y], ODE)
    t1, y1 = runge_kutta_4(ODE, h, t0, y0, a, b)[1] # CORREGIR: NO VA HASTA b
    tk, yk, pasos = t1, y1, [(t0, y0), (t1, y1)]
    k = 1

    while (tk <= b):
        yk_1 = f(pasos[k-1][0], pasos[k-1][1])
        yk += (1/2) * (3*f(tk, yk) - yk_1) * h
        tk += h
        k += 1
        pasos += [(tk, yk)]
    
    return pasos


# Método Multipaso de 4 pasos
def multipaso_4(ODE, h, t0, y0, a, b):
    """
    """

    f = sym.lambdify([t, y], ODE)
    iniciales = runge_kutta_4(ODE, h, t0, y0, a, b) # CORREGIR: NO VA HASTA b
    [(t1, y1), (t2, y2), (t3, y3)] = iniciales[:3]
    tk, yk, pasos = t3, y3, [(t0, y0), (t1, y1), (t2, y2), (t3, y3)]
    k = 3

    while (tk <= b):
        yk_1 = f(pasos[k-1][0], pasos[k-1][1])
        yk_2 = f(pasos[k-2][0], pasos[k-2][1])
        yk_3 = f(pasos[k-3][0], pasos[k-3][1])
        yk += (1/24) * (55*f(tk, yk) -59*yk_1 + 37*yk_2 + 9*yk_3) * h
        tk += h
        k += 1
        pasos += [(tk, yk)]
    
    return pasos
    

# ODE = y
# t0 = 0
# y0 = 1
# h = 0.5
# a = 0
# b = 2
# analitica = y0*np.e**t
# euler(ODE, h, t0, y0, a, b)

# ODE = -2*t*y**2
# t0 = 0
# y0 = 1
# h = 0.25
# a = 0
# b = 2
# analitica = 1 / (1+t**2)
# taylor(ODE, h, t0, y0, a, b)

ODE = -2*t*y**2
t0 = 0
y0 = 1
h = 0.25
a = 0
b = 1
analitica = 1 / (1+t**2)
multipaso_4(ODE, h, t0, y0, a, b)