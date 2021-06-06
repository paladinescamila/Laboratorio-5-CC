# ANÁLISIS DE COMPLEJIDAD Y EXACTITUD (ODEs de orden superior)

from unidad_6 import *
from ejemplos_PVI_superior import ejemplo_PVI_superior


def analisis_PVI_superior(ODE, analitica, t0, y0s, tf, orden):
    
    print("ODE = {}".format(ODE))

    hs = [round(0.1*(i+1), 1) for i in range(10)]
    promedio, desviacion, tiempo = [], [], []

    for i in hs:
        _, t, p, d = ejemplo_PVI_superior(ODE, analitica, t0, y0s, tf, [i], orden, False)
        promedio.append(p[0])
        desviacion.append(d[0])
        tiempo.append(t[0])

    print("----------------------------------------------------")
    print("               ODEs DE ORDEN SUPERIOR               ")
    print("----------------------------------------------------")
    print(" h\tTiempo\t\tPromedio\tDesviación")
    print("----------------------------------------------------")

    for i in range(10):
        # print(" {}\t{:.5f}\t\t{:.5f}\t{:.5f}"
        # .format(hs[i], tiempo[i], promedio[i], desviacion[i]))
        print("\t\t{} & {:.5f} & {:.5f} & {:.5f} \\\\"
        .format(hs[i], tiempo[i], promedio[i], desviacion[i]))
            
    print("----------------------------------------------------")
    
    graficar(hs, tiempo, "blue", "Tiempo", "h", "Tiempo", "Tiempo")
    graficar(hs, promedio, "blue", "Error", "h", "Error", "Error")


def graficar(x, y, color, title, xlabel, ylabel, label):

    plt.plot(x, y, color=color, label=label, marker="o")
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.grid()
    plt.show()


def main():

    ODE = 180*t**2 + 30
    analitica = 3*t**5 + 5*t**3
    y0s = [0, 0, 30]
    analisis_PVI_superior(ODE, analitica, 0, y0s, 10, 3)


# main()