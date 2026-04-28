# Algoritmo de descomposicion dual para el despacho economico de plantas termicas

## Resumen

Resuelve el problema de despacho de plantas termicas usando el algoritmo de dual descomposition

## Parametros

Los parametros del sistema se presentan a continuacion:


| Uni   | $a$     |  $b$    | $p_{max} |
| ----: | ------: | ------: | --------: | 
| 1     | 0.30494 | 38.539  |   100     | 
| 2     | 0.21174 | 46.1591 |   100     | 
| 3     | 0.07092 | 38.3055 |   200     | 
| 4     | 0.05606 | 40.3965 |   200     | 
| 5     | 0.03598 | 38.2704 |   300     | 
| 6     | 0.04222 | 36.3278 |   500     | 


Cada modelo de optimizacion tiene la siguiente forma:

$p_i \leftarrow \underset{p_\text{min}\leq p_i \leq p_\text{max}}{\text{argmin}} \left\\{ \frac{a_i}{2}p_i^2 + b_ip_i + \lambda p_i \right\\}$



## Resultados

Iteracion 1 p = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

Iteracion 2 p = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

Iteracion 3 p = [70.38, 65.37, 200.0, 200.0, 300.0, 500.0]

Iteracion 4 p = [59.25, 49.34, 200.0, 200.0, 300.0, 480.31]

Iteracion 5 p = [51.96, 38.84, 200.0, 200.0, 300.0, 427.67]

Iteracion 6 p = [50.45, 36.66, 200.0, 200.0, 300.0, 416.73]

Iteracion 7 p = [50.13, 36.21, 200.0, 200.0, 300.0, 414.46]

Iteracion 8 p = [50.07, 36.12, 200.0, 200.0, 300.0, 413.98]

Iteracion 9 p = [50.05, 36.1, 200.0, 200.0, 300.0, 413.89]

Iteracion 10 p = [50.05, 36.09, 200.0, 200.0, 300.0, 413.87]



---
## Contacto

Alejandro Garcés Ruiz
(https://github.com/alejandrogarces)

## Licencia

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC_BY--NC--SA_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
