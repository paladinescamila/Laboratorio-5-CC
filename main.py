# EJEMPLOS DE APLICACIÓN DE LOS MÉTODOS

from metodos import *
from PVI_1 import ejemplo_PVI_1, analisis_PVI_1
from PVI_superior import ejemplo_PVI_superior, analisis_PVI_superior
from PVF import ejemplo_PVF, analisis_PVF


def main():

    # ODEs DE PRIMER ORDEN (PROBLEMA DE VALOR INICIAL)

    print("EJEMPLO 1")
    ODE = -2*y + 3*cos(t)
    analitica = 3*sin(t)/5 + 6*cos(t)/5
    hs = [0.09375, 0.1875, 0.375, 0.75]
    ejemplo_PVI_1(ODE, analitica, 5, -0.23495994224201155, 8, hs, True)
    analisis_PVI_1(ODE, analitica, 5, -0.23495994224201155, 8)

    print("EJEMPLO 2")
    ODE = 2*t**3 - 5*t**2
    analitica = t**4/2 - 5*t**3/3
    hs = [1.09375, 2.1875, 4.375, 8.75]
    ejemplo_PVI_1(ODE, analitica, -15, 30937.5, 20, hs, True)
    analisis_PVI_1(ODE, analitica, -15, 30937.5, 20)

    print("EJEMPLO 3")
    ODE = 2*t*cos(t**2) - 3*cos(t)
    analitica = -3*sin(t) + sin(t**2)
    hs = [0.03125, 0.0625, 0.125, 0.25]
    ejemplo_PVI_1(ODE, analitica, 2, -3.4846947757849733, 3, hs, True)
    analisis_PVI_1(ODE, analitica, 2, -3.4846947757849733, 3)


    # ODEs de ORDEN SUPERIOR (PROBLEMA DE VALOR INICIAL)
    
    ODE = 180*t**2 + 30
    analitica = 3*t**5 + 5*t**3
    y0s = [0, 0, 30]
    hs = [0.3125, 0.625, 1.25, 2.5]
    ejemplo_PVI_superior(ODE, analitica, 0, y0s, 10, hs, 3, True)
    analisis_PVI_superior(ODE, analitica, 0, y0s, 10, 3)


    # ODEs bajo un Problema de Valor de Frontera

    ODE = 42*t**5 + 2
    analitica = t**7 + t**2
    ns = [4, 6, 8, 10]
    ejemplo_PVF(ODE, analitica, 0, 0, 15, 170859600, ns, True)
    analisis_PVF(ODE, analitica, 0, 0, 15, 170859600)


main()