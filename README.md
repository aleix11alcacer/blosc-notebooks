# Trabajo desarrollado

Este repositorio contiene todo el código que he desarrollado durante mi estancia en prácticas.

Todo el código que genera resultados está escrito en notebooks de jupyter para poder visualizar también los resultados sin necesidad de ejecutarse.

He organizado todos los archivos por quincenas para facilitar la redacción de los informes quincenales y la memoria final.


## Primera quincena

- `pycblosc2/`. Direcotrio que contiene los ficheros de código del paquete *pycblosc2*.

- `shuffle-iterativo/`. Directorio que contiene un análisis de como afecta a los datos el shuffle iterativo de *Blosc*.

- `orden-permutaciones/`. Directorio que contiene un análisis del orden del shuffle iterativo, pues se puede considerar como una permutación de datos.


## Segunda quincena

- `analisis-imagenes/`. Directorio que contiene un análisis de los resultados de la compresión de imágenes, cuando se aplica una transformación previa a la imagen.

- `analisis-datos-reales/`. Directorio que contiene un análisis de los resultados de la compresión de datos reales, cuando se aplica una transformación previa al conjunto de datos.

- `eliminar-bucles/`. Directorio que contiene un algoritmo que realiza transformaciones a conjuntos de datos en *n* dimensiones sin usar bucles anidados.

## Tercera quincena

- `eliminar-bucles_C/`. Directorio que contiene un algoritmo que realiza transformaciones a conjuntos de datos en *n* dimensiones sin usar bucles anidados y escrito en el lenguaje de programación *C*.

- `arreglar-algoritmo/`. Directorio que contiene los resultados finales, junto con unos tests, del algoritmo creado. Se ha creado para detectar y corregir los posibles fallos del algoritmo.

- `algortimo-general/`. Directorio que contiene el resultado final del algoritmo. Contiene el algoritmo implementado en *numpy*, el algortimo implementado en *Python* y el algortimo implementado en *C*.

## Cuarta quincena

-`algoritmo-bucles-for/`. Directorio que contiene el algoritmo implementado, tanto la versión simple como la general, con bucles *for* para ver si la implementación resulta más sencilla o no. NO SE VA A REALIZAR.

-`resultados/`. Directorio que contiene un notebook que utiliza todo lo que se ha desarrollado hasta la fecha.

-`version-completa/`. Directorio que contiene la libreria que junta el algoritmo simple y el general en una sola función.
