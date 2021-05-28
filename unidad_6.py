# UNIDAD 6: ECUACIONES DIFERENCIALES ORDINARIAS

import time
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt

t = sym.Symbol("t")
y = sym.Symbol("y")


# Método de Euler
def euler(ODE, t0, y0, h, n):
    """
    """

    f = sym.lambdify([t, y], ODE)
    pasos = [(t0, y0)]

    for k in range(n):
        tk, yk = pasos[-1]
        yk += f(tk, yk) * h
        tk += h
        pasos.append((tk, yk))
    
    return pasos


# Método de Serie de Taylor
def taylor(ODE, t0, y0, h, n):
    """
    """

    f = sym.lambdify([t, y], ODE)
    pasos = [(t0, y0)]

    ypp = sym.diff(ODE, t) + sym.diff(ODE, y)*ODE
    fpp = sym.lambdify([t, y], ypp)

    for k in range(n):
        tk, yk = pasos[-1]
        yk += f(tk, yk) * h + (fpp(tk, yk)/2) * h**2
        tk += h
        pasos.append((tk, yk))
    
    return pasos


# Método de Runge-Kutta (orden 2)
def runge_kutta_2(ODE, t0, y0, h, n):
    """
    """

    f = sym.lambdify([t, y], ODE)
    pasos = [(t0, y0)]

    for k in range(n):
        tk, yk = pasos[-1]
        k1 = f(tk, yk)*h
        k2 = f(tk + h, yk + k1)*h
        yk += (1/2) * (k1 + k2)
        tk += h
        pasos.append((tk, yk))
    
    return pasos


# Método de Runge-Kutta (orden 4)
def runge_kutta_4(ODE, t0, y0, h, n):
    """
    """

    f = sym.lambdify([t, y], ODE)
    pasos = [(t0, y0)]

    for k in range(n):
        tk, yk = pasos[-1]
        k1 = f(tk, yk)*h
        k2 = f(tk + h/2, yk + k1/2)*h
        k3 = f(tk + h/2, yk + k2/2)*h
        k4 = f(tk + h, yk + k3)*h
        yk += (1/6) * (k1 + 2*k2 + 2*k3 + k4)
        tk += h
        pasos.append((tk, yk))
    
    return pasos


# Método Multipaso (2 pasos)
def multipaso_2(ODE, t0, y0, h, n):
    """
    """

    f = sym.lambdify([t, y], ODE)
    t1, y1 = runge_kutta_4(ODE, t0, y0, h, 1)[1]
    pasos = [(t0, y0), (t1, y1)]

    for k in range(1, n):
        tk, yk = pasos[-1]
        tk_1, yk_1 = pasos[-2]
        yk += (1/2) * (3*f(tk, yk) - yk_1) * h
        tk += h
        pasos.append((tk, yk))
    
    return pasos


# Método Multipaso (4 pasos)
def multipaso_4(ODE, t0, y0, h, n):
    """
    """

    f = sym.lambdify([t, y], ODE)
    [(t1, y1), (t2, y2), (t3, y3)] = runge_kutta_4(ODE, t0, y0, h, 3)[1:4]
    pasos = [(t0, y0), (t1, y1), (t2, y2), (t3, y3)]

    for k in range(3, n):
        tk, yk = pasos[-1]
        yk_1 = f(pasos[-2][0], pasos[-2][1])
        yk_2 = f(pasos[-3][0], pasos[-3][1])
        yk_3 = f(pasos[-4][0], pasos[-4][1])
        yk += (1/24) * (55*f(tk, yk) -59*yk_1 + 37*yk_2 + 9*yk_3) * h
        tk += h
        pasos.append((tk, yk))
    
    return pasos


# Pintar ejemplo de ODE con cada método (problema de valor inicial)
def ejemplo(ODE, analitica, t0, y0, hs, n, mostrar):
    """
    """

    x_funcion = np.linspace(t0, t0 + max(hs)*n, 1000)
    f_analitica = sym.lambdify([t, y], analitica)
    y_analitica = [f_analitica(i,0) for i in x_funcion]

    nh = len(hs)
    resultados = [[None for _ in range(nh)] for _ in range(6)]
    tiempos = [[None for _ in range(nh)] for _ in range(6)]
    promedios = [[None for _ in range(nh)] for _ in range(6)]
    desviaciones = [[None for _ in range(nh)] for _ in range(6)]

    if (mostrar):
        
        print("ODE = {}\ny(t) = {}".format(ODE, analitica))
        print("t0 = {}\ny0 = {}\nhs = {}\nn = {}".format(t0, y0, hs, n))

        t_metodos, p_metodos, d_metodos, f_metodos, y_metodos = [], [], [], [], []
        colores = ["purple", "red", "green", "orange", "gray"]
        metodos = ["Euler", "Taylor", "Runge-Kutta 2"]
        if (nh > 5): colores = ["purple" for _ in range(nh)]


def main():

    ODE = y
    analitica = 0*np.e**t
    print(euler(ODE, 0, 1, 0.5, 5))
    # ejemplo(ODE, analitica, 0, 1, [0.5], 5, True)

    ODE = -2*t*y**2
    analitica = 1 / (1+t**2)
    print(taylor(ODE, 0, 1, 0.25, 5))

    ODE = -2*t*y**2
    analitica = 1 / (1+t**2)
    print(multipaso_4(ODE, 0, 1, 0.25, 5))


main()