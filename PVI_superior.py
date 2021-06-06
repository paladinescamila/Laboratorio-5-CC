# ODEs DE ORDEN SUPERIOR (PROBLEMAS DE VALOR INICIAL)

from metodos import *


# Muestra ejemplos de ODEs de Orden Superior
def ejemplo_PVI_superior(ODE, analitica, t0, y0s, tf, hs, n, mostrar):
    """
    Entrada: una Ecuación Diferencial Ordinaria ODE de orden superior, la solución 
            analítica de la ODE, un real t0 y una lista de reales y0s, que 
            representan las condiciones iniciales para cada ODE de orden 1, un real 
            tf que es el límite superior del intervalo de análisis, una lista de 
            reales hs que son el incremento de t, un entero n que representa el 
            orden de la ODE, y un booleano mostrar.
    Salida: los puntos (tk, yk) que corresponden a la solución de la ODE 
            (t0 <= tk <= tf) para cada h, sus tiempos de ejecución, el promedio del 
            error absoluto y la desviación del error.
    """

    x_funcion = np.linspace(t0, tf, 1000)
    f_analitica = sym.lambdify(t, analitica)
    y_analitica = [f_analitica(i) for i in x_funcion]

    nh = len(hs)
    resultados, tiempos, promedios, desviaciones = [], [], [], []

    if (mostrar):
        
        print("ODE = {}\ny(t) = {}".format(ODE, analitica))
        for i in range(n):
            print("y"+"'"*i+"({}) = {}".format(t0, y0s[i]))
        print("tf = {}\nhs = {}".format(tf, hs))

        colores = ["red", "blue", "green", "purple", "orange", "dodgerblue"]
        if (nh > 6): colores = ["blue" for _ in range(nh)]
    
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
        plt.xlabel('t')
        plt.ylabel('y')
        plt.grid()
        plt.show()
    
    return resultados, tiempos, promedios, desviaciones


# Muestra análisis de ejemplos de ODEs de Orden Superior
def analisis_PVI_superior(ODE, analitica, t0, y0s, tf, orden):
    
    print("ODE = {}".format(ODE))

    hs = [round(0.1*(i+1), 1) for i in range(10)]
    tiempo, promedio, desviacion = [], [], []

    for i in hs:
        _, t, p, d = ejemplo_PVI_superior(ODE, analitica, t0, y0s, tf, [i], orden, False)
        tiempo.append(t[0])
        promedio.append(p[0])
        desviacion.append(d[0])

    print("----------------------------------------------------")
    print("               ODEs DE ORDEN SUPERIOR               ")
    print("----------------------------------------------------")
    print(" h\tTiempo\t\tPromedio\tDesviación")
    print("----------------------------------------------------")

    for i in range(10):
        print(" {}\t{:.5f}\t\t{:.5f}\t{:.5f}"
        .format(hs[i], tiempo[i], promedio[i], desviacion[i]))
            
    print("----------------------------------------------------")
    
    graficar_PVI_superior(hs, tiempo, "blue", "Tiempo", "h", "Tiempo", "Tiempo")
    graficar_PVI_superior(hs, promedio, "blue", "Error", "h", "Error", "Error")


def graficar_PVI_superior(x, y, color, title, xlabel, ylabel, label):

    plt.plot(x, y, color=color, label=label, marker="o")
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.grid()
    plt.show()