# ODEs DE ORDEN 1 (PROBLEMAS DE VALOR INICIAL)

from metodos import *


# Muestra ejemplos de ODEs de Primer Orden
def ejemplo_PVI_1(ODE, analitica, t0, y0, tf, hs, mostrar):
    """
    Entrada: una Ecuación Diferencial Ordinaria ODE de orden 1, la solución 
            analítica de la ODE, dos reales t0, y0, que representan las 
            condiciones iniciales y(t0) = y0, un real tf que es el límite 
            superior del intervalo de análisis, una lista de reales hs que son 
            el incremento de t, y un booleano mostrar.
    Salida: los puntos (tk, yk) que corresponden a la solución de la ODE para 
            cada método y cada h, sus tiempos de ejecución, el promedio del 
            error absoluto y la desviación del error.
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

        y_metodos, t_metodos, p_metodos, d_metodos = [], [], [], []
        colores = ["red", "blue", "green", "purple", "orange", "dodgerblue"]
        metodos = ["Euler", "Taylor", "Runge-Kutta 2", "Runge-Kutta 4", 
                   "Multipaso 2", "Multipaso 4"]
        if (nh > 6): colores = ["blue" for _ in range(nh)]

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

            if (j == nh - 1 and mostrar):
                t_metodos.append(tiempo)
                p_metodos.append(promedio)
                d_metodos.append(desviacion)
                y_metodos.append(pasos)
            
            resultados[i][j], tiempos[i][j] = pasos, tiempo
            promedios[i][j], desviaciones[i][j] = promedio, desviacion

            if (mostrar):
                print(" {}\t{:.10f}\t{:.10f}\t{:.10f}"
                .format(hs[j], tiempo, promedio, desviacion))
                ts = [ti for ti, yi in pasos]
                ys = [yi for ti, yi in pasos]
                plt.plot(ts, ys, color=colores[j], label="h = "+str(hs[j]), marker="o", markersize=4)

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
        plt.title("Métodos con h = "+str(hs[-1]))
        print("------------------------------------------------------")
        print("             PROBLEMAS DE VALOR INICIAL              ")
        print("------------------------------------------------------")
        print(" Método\t\tTiempo\tError (Prom)\tError (Desv)")
        print("------------------------------------------------------")
        for i in range(6):
            print(" {}\t{:.10f}\t{:.10f}\t{:.10f}"
            .format(metodos[i], t_metodos[i], p_metodos[i], d_metodos[i]))
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


# Muestra análisis de ejemplos de ODEs de Primer Orden
def analisis_PVI_1(ODE, analitica, t0, y0, tf):
    
    print("ODE = {}".format(ODE))

    colores = ["red", "blue", "green", "purple", "orange", "dodgerblue"]
    metodos = ["Euler", "Taylor", "Runge-Kutta 2", "Runge-Kutta 4", 
                "Multipaso 2", "Multipaso 4"]

    hs = [round(0.1*(i+1), 1) for i in range(10)]
    tiempo = [[] for _ in range(6)]
    promedio = [[] for _ in range(6)]
    desviacion = [[] for _ in range(6)]

    for i in hs:
        _, t, p, d = ejemplo_PVI_1(ODE, analitica, t0, y0, tf, [i], False)
        for j in range(6):
            tiempo[j].append(t[j][0])
            promedio[j].append(p[j][0])
            desviacion[j].append(d[j][0])

    imprimir_PVI_1("Tiempo", hs, tiempo, ["h"] + metodos)
    graficar_PVI_1(hs, tiempo, colores, "Tiempo", "h", "Tiempo", metodos)

    imprimir_PVI_1("Error (Promedio)", hs, promedio, ["h"] + metodos)
    imprimir_PVI_1("Error (Desviación)", hs, desviacion, ["h"] + metodos)
    graficar_PVI_1(hs, promedio, colores, "Error", "h", "Error", metodos)


def graficar_PVI_1(x, y, color, title, xlabel, ylabel, label):

    for i in range(len(y)): 
        plt.plot(x, y[i], color=color[i], label=label[i], marker="o")

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.grid()
    plt.show()


def imprimir_PVI_1(titulo, x, y, c):

    print("------------------------------------------------------")
    print(" " + titulo)
    print("------------------------------------------------------")
    print(" {}\t{}\t{}\t{}\t{}\t{}\t{}".format(c[0], c[1], c[2], c[3], c[4], c[5], c[6]))
    print("------------------------------------------------------")

    for i in range(len(x)):
        y1, y2, y3, y4, y5, y6 = y[0][i], y[1][i], y[2][i], y[3][i], y[4][i], y[5][i]
        print(" {}\t{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}"
        .format(x[i], y1, y2, y3, y4, y5, y6))
            
    print("------------------------------------------------------")