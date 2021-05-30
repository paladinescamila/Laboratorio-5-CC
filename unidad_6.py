# UNIDAD 6: ECUACIONES DIFERENCIALES ORDINARIAS

import time
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt

t = sym.Symbol("t")
y = sym.Symbol("y")


# Método de Euler
def euler(ODE, t0, y0, h, tn):
    """
    Entrada: una Ecuación Diferencial Ordinaria de primer orden ODE(t, y), dos reales
            t0 & y0 que representan las condiciones iniciales y(t0) = y0, un real
            h que representa el incremento de t, y un real tn que es el límite
            superior del intervalo de aplicación.
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
            superior del intervalo de aplicación.
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
            superior del intervalo de aplicación.
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
            superior del intervalo de aplicación.
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
            superior del intervalo de aplicación.
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
            superior del intervalo de aplicación.
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
def ODEs_superior(ODE, t0s, y0s, hs, ns, orden):
    """
    """
    pasos = []
    for i in range(orden):
        pasos = runge_kutta_4(ODE, t0s[i], y0s[i], hs[i], ns[i])

    print("aiuda")


# Pintar ejemplo de ODE con cada método (problema de valor inicial)
def ejemplo(ODE, analitica, t0, y0, hs, tn, mostrar):
    """
    """

    x_funcion = np.linspace(t0, tn, 1000)
    f_analitica = sym.lambdify(t, analitica)
    y_analitica = [f_analitica(i) for i in x_funcion]

    nh = len(hs)
    resultados = [[None for _ in range(nh)] for _ in range(6)]
    tiempos = [[None for _ in range(nh)] for _ in range(6)]
    promedios = [[None for _ in range(nh)] for _ in range(6)]
    desviaciones = [[None for _ in range(nh)] for _ in range(6)]

    if (mostrar):
        
        print("ODE = {}\ny(t) = {}".format(ODE, analitica))
        print("y({}) = {}\ntn = {}\nhs = {}".format(t0, y0, tn, hs))

        t_metodos, p_metodos, d_metodos, y_metodos = [], [], [], []
        colores = ["red", "purple", "green", "blue", "orange", "gray"]
        metodos = ["Euler", "Taylor", "Runge-Kutta 2", "Runge-Kutta 4", 
                   "Multipaso 2", "Multipaso 4"]
        if (nh > 6): colores = ["purple" for _ in range(nh)]

    for i in range(6):

        if (mostrar):

            plt.title("Método de " + metodos[i])
            print("----------------------------------------------")
            print(" MÉTODO DE {}".format(metodos[i].upper()))
            print("----------------------------------------------")
            print(" h\tTiempo\t\tError (Prom)\tError (Desv)")
            print("----------------------------------------------")

        for j in range(nh):

            if (i == 0):
                inicio = time.time()
                pasos = euler(ODE, t0, y0, hs[j], tn)
                tiempo = time.time() - inicio

            elif (i == 1):
                inicio = time.time()
                pasos = taylor(ODE, t0, y0, hs[j], tn)
                tiempo = time.time() - inicio

            elif (i == 2):
                inicio = time.time()
                pasos = runge_kutta_2(ODE, t0, y0, hs[j], tn)
                tiempo = time.time() - inicio
            
            elif (i == 3):
                inicio = time.time()
                pasos = runge_kutta_4(ODE, t0, y0, hs[j], tn)
                tiempo = time.time() - inicio
            
            elif (i == 4):
                inicio = time.time()
                pasos = multipaso_2(ODE, t0, y0, hs[j], tn)
                tiempo = time.time() - inicio
            
            else:
                inicio = time.time()
                pasos = multipaso_4(ODE, t0, y0, hs[j], tn)
                tiempo = time.time() - inicio

            errores = [np.abs(yi - f_analitica(ti)) for ti, yi in pasos]
            errores.pop(0)
            promedio, desviacion = np.mean(errores), np.std(errores)

            if (j == 0 and mostrar):
                t_metodos.append(tiempo)
                p_metodos.append(promedio)
                d_metodos.append(desviacion)
                y_metodos.append(pasos)
            
            resultados[i][j], tiempos[i][j] = pasos, tiempo
            promedios[i][j], desviaciones[i][j] = promedio, desviacion

            if (mostrar):
                print(" {}\t{:.10f}\t{:.10f}\t{:.10f}"
                .format(hs[j], tiempo, promedio, desviacion))
                # plt.plot(t0, y0, color=colores[j], label="h = "+str(hs[j]), marker="o", markersize=4)
                # for ti, yi in pasos:
                #     plt.plot(ti, yi, color=colores[j], marker="o", markersize=4)
                ts = [ti for ti, yi in pasos]
                ys = [yi for ti, yi in pasos]
                plt.plot(ts, ys, color=colores[j], label="h = "+str(hs[j]), marker="o", markersize=4)

        if (mostrar):

            print("----------------------------------------------\n")   
            plt.plot(x_funcion, y_analitica, color="black", label="Analítica")
            plt.legend()
            plt.xlabel('x')
            plt.ylabel('y')
            plt.grid()
            plt.show()

    if (mostrar):

        print()
        plt.title("Métodos con h = "+str(hs[0]))
        print("-----------------------------------------------------")
        print("             PROBLEMAS DE VALOR INICIAL              ")
        print("-----------------------------------------------------")
        print(" Método\t\tTiempo\tError (Prom)\tError (Desv)")
        print("-----------------------------------------------------")
        for i in range(6):
            print(" {}\t{:.10f}\t{:.10f}\t{:.10f}"
            .format(metodos[i], t_metodos[i], p_metodos[i], d_metodos[i]))
            # plt.plot(t0, y0, color=colores[i], label=str(metodos[i]), marker="o", markersize=4)
            # for ti, yi in y_metodos[i]:
            #     plt.plot(ti, yi, color=colores[i], marker="o", markersize=4)
            ts = [ti for ti, yi in y_metodos[i]]
            ys = [yi for ti, yi in y_metodos[i]]
            plt.plot(ts, ys, color=colores[i], label=str(metodos[i]), marker="o", markersize=4)
        print("-----------------------------------------------------")
        
        plt.plot(x_funcion, y_analitica, color="black", label="Analítica", linewidth=2)
        plt.legend()
        plt.xlabel('x')
        plt.ylabel('y')
        plt.grid()
        plt.show()
        print()

    return analitica, resultados, tiempos, promedios, desviaciones


def main():

    ODE = y
    analitica = np.e**t
    hs = [0.1, 0.2, 0.3, 0.4]
    ejemplo(ODE, analitica, 0, 1, hs, 5, True)

    # ODE = -2*t*y**2
    # analitica = 1 / (1+t**2)
    # print(taylor(ODE, 0, 1, 0.25, 0.5))

    # ODE = -2*t*y**2
    # analitica = 1 / (1+t**2)
    # print(multipaso_4(ODE, 0, 1, 0.25, 1))


main()