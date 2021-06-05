# ANÁLISIS DE COMPLEJIDAD Y EXACTITUD (ODEs DE UN PVF)

from unidad_6 import *
from ejemplos_PVF import ejemplo_PVF


def analisis_PVF(ODE, analitica, t0, y0, tf, yf):
    
    print("ODE = {}".format(ODE))

    colores = ["blue", "red"]
    metodos = ["Diferencias Finitas", "Elementos Finitos"]

    ns = [(i+3) for i in range(10)]
    promedio = [[] for _ in range(2)]
    desviacion = [[] for _ in range(2)]
    tiempo = [[] for _ in range(2)]

    for i in ns:
        _, t, p, d = ejemplo_PVF(ODE, analitica, t0, y0, tf, yf, [i], False)
        for j in range(2):
            promedio[j].append(p[j][0])
            desviacion[j].append(d[j][0])
            tiempo[j].append(t[j][0])

    imprimir("Error (Promedio)", ns, promedio, ["n"] + metodos)
    imprimir("Error (Desviación)", ns, desviacion, ["n"] + metodos)
    graficar(ns, promedio, colores, "Error", "n", "Error", metodos)

    imprimir("Tiempo", ns, tiempo, ["n"] + metodos)
    graficar(ns, tiempo, colores, "Tiempo", "n", "Tiempo", metodos)


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
    print(" {}\t{}\t{}".format(c[0], c[1], c[2]))
    print("------------------------------------------------------")

    for i in range(len(x)):
        y1, y2 = y[0][i], y[1][i]
        print(" {}\t{:.10f}\t{:.10f}".format(x[i], y1, y2))
            
    print("------------------------------------------------------")


def main():

    ODE = 42*t**5 + 2
    analitica = t**7 + t**2 + 10*t + 2
    analisis_PVF(ODE, analitica, 0, 2, 2, 162)


main()