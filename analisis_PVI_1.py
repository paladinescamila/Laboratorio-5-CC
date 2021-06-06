# ANÁLISIS DE COMPLEJIDAD Y EXACTITUD (ODEs de orden 1)

from unidad_6 import *
from ejemplos_PVI_1 import ejemplo_PVI_1


def analisis_PVI_1(ODE, analitica, t0, y0, tf):
    
    print("ODE = {}".format(ODE))

    colores = ["red", "blue", "green", "purple", "orange", "dodgerblue"]
    metodos = ["Euler", "Taylor", "Runge-Kutta 2", "Runge-Kutta 4", 
                "Multipaso 2", "Multipaso 4"]

    hs = [round(0.1*(i+1), 1) for i in range(10)]
    promedio = [[] for _ in range(6)]
    desviacion = [[] for _ in range(6)]
    tiempo = [[] for _ in range(6)]

    for i in hs:
        _, t, p, d = ejemplo_PVI_1(ODE, analitica, t0, y0, tf, [i], False)
        for j in range(6):
            promedio[j].append(p[j][0])
            desviacion[j].append(d[j][0])
            tiempo[j].append(t[j][0])

    imprimir("Tiempo", hs, tiempo, ["h"] + metodos)
    graficar(hs, tiempo, colores, "Tiempo", "h", "Tiempo", metodos)

    imprimir("Error (Promedio)", hs, promedio, ["h"] + metodos)
    imprimir("Error (Desviación)", hs, desviacion, ["h"] + metodos)
    graficar(hs, promedio, colores, "Error", "h", "Error", metodos)


def graficar(x, y, color, title, xlabel, ylabel, label):

    for i in range(len(y)): 
        plt.plot(x, y[i], color=color[i], label=label[i], marker="o")

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.grid()
    plt.show()


def imprimir(titulo, x, y, c):

    print("------------------------------------------------------")
    print(" " + titulo)
    print("------------------------------------------------------")
    print(" {}\t{}\t{}\t{}\t{}\t{}\t{}".format(c[0], c[1], c[2], c[3], c[4], c[5], c[6]))
    print("------------------------------------------------------")

    for i in range(len(x)):
        y1, y2, y3, y4, y5, y6 = y[0][i], y[1][i], y[2][i], y[3][i], y[4][i], y[5][i]
        # print(" {}\t{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}".format(x[i], y1, y2, y3, y4, y5, y6))
        print("\t\t{} & {:.5f} & {:.5f} & {:.5f} & {:.5f} & {:.5f} & {:.5f} \\\\"
        .format(x[i], y1, y2, y3, y4, y5, y6))
            
    print("------------------------------------------------------")


def main():

    print("EJEMPLO 1")
    ODE = -2*y + 3*sym.cos(t)
    analitica = 3*sym.sin(t)/5 + 6*sym.cos(t)/5
    analisis_PVI_1(ODE, analitica, 5, -0.23495994224201155, 8)

    print("EJEMPLO 2")
    ODE = 2*t**3 - 5*t**2
    analitica = t**4/2 - 5*t**3/3
    analisis_PVI_1(ODE, analitica, -15, 30937.5, 20)

    print("EJEMPLO 3")
    ODE = 2*t*sym.cos(t**2) - 3*sym.cos(t)
    analitica = -3*sym.sin(t) + sym.sin(t**2)
    analisis_PVI_1(ODE, analitica, 2, -3.4846947757849733, 3)


# main()