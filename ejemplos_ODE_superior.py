# EJEMPLOS DE ODEs DE ORDEN SUPERIOR (PROBLEMAS DE VALOR INICIAL)

from unidad_6 import *


def ejemplo_ODE_superior(ODE, analitica, t0, y0s, tf, hs, n, mostrar):
    """
    """

    x_funcion = np.linspace(t0, tf, 1000)
    f_analitica = sym.lambdify(t, analitica)
    y_analitica = [f_analitica(i) for i in x_funcion]

    nh = len(hs)
    resultados, tiempos, promedios, desviaciones = [], [], [], []

    if (mostrar):
        
        print("ODE = {}\ny".format(ODE)+"'"*n+" = {}".format(analitica))
        for i in range(n):
            print("y"+"'"*i+"({}) = {}".format(t0, y0s[i]))
        print("tf = {}\nhs = {}".format(tf, hs))

        colores = ["red", "purple", "green", "blue", "orange", "gray"]
        if (nh > 6): colores = ["purple" for _ in range(nh)]
    
        plt.title("ODEs de orden superior")
        print("------------------------------------------------------")
        print("              ODEs DE ORDEN SUPERIOR {}              ".format(n))
        print("------------------------------------------------------")
        print(" h\tTiempo\t\tError (Prom)\tError (Desv)")
        print("------------------------------------------------------")

    for i in range(nh):
        
        inicio = time.time()
        pasos = ODEs_superior(ODE, t0, y0s, tf, hs[i], n)
        tiempo = time.time() - inicio

        errores = [np.abs(yi - f_analitica(ti)) for ti, yi in pasos]
        errores.pop(0)
        promedio, desviacion = np.mean(errores), np.std(errores)

        resultados.append(pasos)
        tiempos.append(tiempo)
        promedios.append(promedio)
        desviaciones.append(desviacion)

        if (mostrar):
                print(" {}\t{:.10f}\t{:.10f}\t{:.10f}"
                .format(hs[i], tiempo, promedio, desviacion))
                ts = [ti for ti, yi in pasos]
                ys = [yi for ti, yi in pasos]
                plt.plot(ts, ys, color=colores[i], label="h = "+str(hs[i]), marker="o", markersize=4)

    if (mostrar):

        print("------------------------------------------------------\n")
        plt.plot(x_funcion, y_analitica, color="black", label="Analítica")
        plt.legend()
        plt.xlabel('x')
        plt.ylabel('y')
        plt.grid()
        plt.show()
    
    return resultados, tiempos, promedios, desviaciones


def main():

    ODE = 180*t**2 + 30
    analitica = 3*t**5 + 5*t**3 + 10
    y0s = [0, 0, 30]
    hs = [0.5, 1, 2, 2.5]
    ejemplo_ODE_superior(ODE, analitica, 0, y0s, 10, hs, 3, True)


main()