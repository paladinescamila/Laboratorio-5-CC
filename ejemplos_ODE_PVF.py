# EJEMPLOS DE ODEs  (PROBLEMAS DE VALOR DE FRONTERA)

from unidad_6 import *


def ejemplo_ODE_PVF(ODE, analitica, t0, y0, tf, yf, ns, mostrar):
    """
    """

    x_funcion = np.linspace(t0, tf, 1000)
    f_analitica = sym.lambdify(t, analitica)
    y_analitica = [f_analitica(i) for i in x_funcion]

    nn = len(ns)
    resultados = [[None for _ in range(nn)] for _ in range(2)]
    tiempos = [[None for _ in range(nn)] for _ in range(2)]
    promedios = [[None for _ in range(nn)] for _ in range(2)]
    desviaciones = [[None for _ in range(nn)] for _ in range(2)]

    if (mostrar):
        
        print("ODE = {}\ny(t) = {}".format(ODE, analitica))
        print("y({}) = {}\ny({}) = {}\nns = {}".format(t0, y0, tf, yf, ns))

        t_metodos, p_metodos, d_metodos, y_metodos = [], [], [], []
        colores = ["red", "purple", "green", "blue", "orange", "gray"]
        metodos = ["Diferencias Finitas", "Elementos Finitos"]
        if (nn > 6): colores = ["purple" for _ in range(nn)]

    for i in range(2):

        if (mostrar):

            plt.title("Método de " + metodos[i])
            print("------------------------------------------------------")
            print(" MÉTODO DE {}".format(metodos[i].upper()))
            print("------------------------------------------------------")
            print(" n\tTiempo\t\tError (Prom)\tError (Desv)")
            print("------------------------------------------------------")

        for j in range(nn):

            if (i == 0):
                inicio = time.time()
                puntos = diferencias_finitas(ODE, t0, y0, tf, yf, ns[j])
                tiempo = time.time() - inicio

            else:
                inicio = time.time()
                polinomio = elementos_finitos(ODE, t0, y0, tf, yf, ns[j])
                tiempo = time.time() - inicio
                f_ef = sym.lambdify(t, polinomio)
                puntos = [(i, f_ef(i)) for i in np.linspace(t0, tf, ns[j])]

            errores = [np.abs(yi - f_analitica(ti)) for ti, yi in puntos]
            errores.pop(0)
            errores.pop()
            promedio, desviacion = np.mean(errores), np.std(errores)

            if (j == nn - 1 and mostrar):
                t_metodos.append(tiempo)
                p_metodos.append(promedio)
                d_metodos.append(desviacion)
                y_metodos.append(puntos)
            
            resultados[i][j], tiempos[i][j] = puntos, tiempo
            promedios[i][j], desviaciones[i][j] = promedio, desviacion

            if (mostrar):
                print(" {}\t{:.10f}\t{:.10f}\t{:.10f}"
                .format(ns[j], tiempo, promedio, desviacion))
                # plt.plot(t0, y0, color=colores[j], label="h = "+str(ns[j]), marker="o", markersize=4)
                # for ti, yi in puntos:
                #     plt.plot(ti, yi, color=colores[j], marker="o", markersize=4)
                ts = [ti for ti, yi in puntos]
                ys = [yi for ti, yi in puntos]
                plt.plot(ts, ys, color=colores[j], label="h = "+str(ns[j]), marker="o", markersize=4)

        if (mostrar):

            print("------------------------------------------------------\n")   
            plt.plot(x_funcion, y_analitica, color="black", label="Analítica")
            plt.legend()
            plt.xlabel('t')
            plt.ylabel('y')
            plt.grid()
            plt.show()

    if (mostrar):

        print()
        plt.title("Métodos con n = "+str(ns[-1]))
        print("------------------------------------------------------")
        print("            PROBLEMAS DE VALOR DE FRONTERA            ")
        print("------------------------------------------------------")
        print(" Método\t\tTiempo\tError (Prom)\tError (Desv)")
        print("------------------------------------------------------")
        for i in range(2):
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
        plt.xlabel('t')
        plt.ylabel('y')
        plt.grid()
        plt.show()
        print()

    return resultados, tiempos, promedios, desviaciones


def main():

    ODE = 42*t**5 + 2
    analitica = t**7 + t**2 + 10*t + 2
    ns = [4, 6, 8, 10]
    ejemplo_ODE_PVF(ODE, analitica, 0, 2, 2, 162, ns, True)


main()