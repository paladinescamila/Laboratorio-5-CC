# EJEMPLOS DE ODEs  (PROBLEMAS DE VALOR DE FRONTERA)

from unidad_6 import *


def ejemplo_ODE_PVF(ODE, analitica, t0, y0, tf, yf, n, mostrar):
    solucion = diferencias_finitas(ODE, t0, y0, tf, yf, n)
    print(solucion)


def main():

    # EJEMPLO DE LA CLASE
    ODE = 6*t
    analitica = t**3
    t0, y0, tf, yf = 0, 0, 1, 1
    ejemplo_ODE_PVF(ODE, analitica, t0, y0, tf, yf, 5, True)

main()