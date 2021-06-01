# EJEMPLOS DE ODEs  (PROBLEMAS DE VALOR DE FRONTERA)

from unidad_6 import *


def ejemplo_ODE_PVF(ODE, analitica, t0, y0, tf, yf, ns, mostrar):
    """
    """

    x_analitica = np.linspace(t0, tf, 1000)
    f_analitica = sym.lambdify(t, analitica)
    y_analitica = [f_analitica(i) for i in x_analitica]

    colores = ["red", "purple", "green", "blue", "orange", "gray"]
    pasos = diferencias_finitas(ODE, t0, y0, tf, yf, ns[0])
    if (mostrar):
        ti = [i for i,j in pasos]
        yi = [j for i,j in pasos]
        plt.plot(ti, yi, color="blue", marker="o", markersize=4)
        plt.plot(x_analitica, y_analitica, color="black")
        for i, j in pasos: print(i, j)


def main():

    # EJEMPLO DE LA CLASE
    ODE = 6*t
    analitica = t**3
    t0, y0, tf, yf = 0, 0, 1, 1
    ejemplo_ODE_PVF(ODE, analitica, t0, y0, tf, yf, [5], True)

main()