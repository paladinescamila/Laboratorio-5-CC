# EJEMPLOS DE ODEs DE ORDEN 1 (PROBLEMAS DE VALOR INICIAL)

from unidad_6 import *


def ejemplo_ODE_1(ODE, analitica, t0, y0, tf, hs, mostrar):
    """
    """

    x_funcion = np.linspace(t0, tf, 1000)
    f_analitica = sym.lambdify(t, analitica)
    y_analitica = [f_analitica(i) for i in x_funcion]

    nh = len(hs)
    resultados = [[None for _ in range(nh)] for _ in range(6)]
    tiempos = [[None for _ in range(nh)] for _ in range(6)]
    promedios = [[None for _ in range(nh)] for _ in range(6)]
    desviaciones = [[None for _ in range(nh)] for _ in range(6)]

    if (mostrar):
        
        print("ODE = {}\ny(t) = {}".format(ODE, analitica))
        print("y({}) = {}\ntf = {}\nhs = {}".format(t0, y0, tf, hs))

        t_metodos, p_metodos, d_metodos, y_metodos = [], [], [], []
        colores = ["red", "purple", "green", "blue", "orange", "gray"]
        metodos = ["Euler", "Taylor", "Runge-Kutta 2", "Runge-Kutta 4", 
                   "Multipaso 2", "Multipaso 4"]
        if (nh > 6): colores = ["purple" for _ in range(nh)]

    for i in range(6):

        if (mostrar):

            plt.title("Método de " + metodos[i])
            print("------------------------------------------------------")
            print(" MÉTODO DE {}".format(metodos[i].upper()))
            print("------------------------------------------------------")
            print(" h\tTiempo\t\tError (Prom)\tError (Desv)")
            print("------------------------------------------------------")

        for j in range(nh):

            if (i == 0):
                inicio = time.time()
                pasos = euler(ODE, t0, y0, tf, hs[j])
                tiempo = time.time() - inicio

            elif (i == 1):
                inicio = time.time()
                pasos = taylor(ODE, t0, y0, tf, hs[j])
                tiempo = time.time() - inicio

            elif (i == 2):
                inicio = time.time()
                pasos = runge_kutta_2(ODE, t0, y0, tf, hs[j])
                tiempo = time.time() - inicio
            
            elif (i == 3):
                inicio = time.time()
                pasos = runge_kutta_4(ODE, t0, y0, tf, hs[j])
                tiempo = time.time() - inicio
            
            elif (i == 4):
                inicio = time.time()
                pasos = multipaso_2(ODE, t0, y0, tf, hs[j])
                tiempo = time.time() - inicio
            
            else:
                inicio = time.time()
                pasos = multipaso_4(ODE, t0, y0, tf, hs[j])
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

            print("------------------------------------------------------\n")   
            plt.plot(x_funcion, y_analitica, color="black", label="Analítica")
            plt.legend()
            plt.xlabel('x')
            plt.ylabel('y')
            plt.grid()
            plt.show()

    if (mostrar):

        print()
        plt.title("Métodos con h = "+str(hs[0]))
        print("------------------------------------------------------")
        print("             PROBLEMAS DE VALOR INICIAL              ")
        print("------------------------------------------------------")
        print(" Método\t\tTiempo\tError (Prom)\tError (Desv)")
        print("------------------------------------------------------")
        for i in range(6):
            print(" {}\t{:.10f}\t{:.10f}\t{:.10f}"
            .format(metodos[i], t_metodos[i], p_metodos[i], d_metodos[i]))
            # plt.plot(t0, y0, color=colores[i], label=str(metodos[i]), marker="o", markersize=4)
            # for ti, yi in y_metodos[i]:
            #     plt.plot(ti, yi, color=colores[i], marker="o", markersize=4)
            ts = [ti for ti, yi in y_metodos[i]]
            ys = [yi for ti, yi in y_metodos[i]]
            plt.plot(ts, ys, color=colores[i], label=str(metodos[i]), marker="o", markersize=4)
        print("------------------------------------------------------")
        
        plt.plot(x_funcion, y_analitica, color="black", label="Analítica", linewidth=2)
        plt.legend()
        plt.xlabel('x')
        plt.ylabel('y')
        plt.grid()
        plt.show()
        print()

    return analitica, resultados, tiempos, promedios, desviaciones


def main():

    # EJEMPLO 1
    ODE = 2*t*sym.cos(t**2) - 2*sym.cos(t)
    analitica = sym.sin(t**2) - 2*sym.sin(t)
    hs = [0.2, 0.5, 0.1, 0.4]
    ejemplo_ODE_1(ODE, analitica, 0, 0, 6, hs, True)

    # EJEMPLO 2
    ODE = 2**(sym.sin(t**2)/t)*(2*sym.cos(t**2) - sym.sin(t**2)/t**2)*sym.log(2)
    analitica = 2**(sym.sin(t**2)/t)
    hs = [0.2, 0.5, 0.1, 0.4]
    ejemplo_ODE_1(ODE, analitica, 1, 1.8, 4, hs, True)

    # EJEMPLO 3
    ODE = (-t*sym.sin(t)*sym.cos(sym.cos(t)) + sym.sin(sym.sin(t))*sym.cos(t) + 
           sym.sin(sym.cos(t))) / ((t*sym.sin(sym.cos(t)) - sym.cos(sym.sin(t)))**2 + 1)
    analitica = sym.atan(t*sym.sin(sym.cos(t)) - sym.cos(sym.sin(t)))
    hs = [0.25, 0.5, 0.15, 0.05]
    ejemplo_ODE_1(ODE, analitica, 7.5, 1.01, 9, hs, True)


main()