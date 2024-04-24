# Script para generar un PWM a partir de datos de nucleotidos  

Script creado por Sebastian Quiroz para Pontificia Universidad Católica de Chile  
Contacto. <squirozpa@uc.cl>

## Utilizacion

El script tiene diferentes funciones, pero para todas es necesario estar en la carpeta inicial del script (Carpeta total del repositorio)
El script contiene 4 funciones principales, las cuales son:
    - PWM(genera un PFM, PRM y PWM a partir de un archivo fasta o raw de un msa)
    - Fuzzy(genera un sequencia fuzzy a partir de la tabla PPM o PRM, de la funcion anterior y un threshold)
    - Alignment(busca un aliniamiento a partir de un fuzzy (o sequencia normal), una secuencia target y un maximo mismatch, y entrega un arhivo con los aliniamientos encontrados y su mismatch)
    - Energy profile (genera un perfil de energia a partir de un PWM y una secuencia target, entrega un archivo y un grafico con la energia de cada posicion de la secuencia target)

Estos scripts toman los archivos de la carpeta input_files y los entrega en output_files, si se desea utilizar un archivo desde el output_files, se debe agregar = antes de el "path" que deberia tomar como input.

Guias de comandos:
() = argumento absoluto
[] = argumento relativo. ej paths de archivos, nombres del output

### PWM

Para generar un PWM se debe correr el script con el siguiente comando:

Functions_general.py (PWM) [path del archivo fasta o raw] [nombre del output]

### Fuzzy

Para generar un fuzzy se debe correr el script con el siguiente comando:

Functions_general.py (Fuzzy) [path del archivo PRM] [threshold] [nombre del output]

### Alignment

Para generar un alineamiento se debe correr el script con el siguiente comando:

Functions_general.py (Alignment) [path del archivo fasta o raw] [path del archivo fasta o raw target] [nombre del output] [maximo mismatch]

### Energy profile

Para generar un perfil de energia se debe correr el script con el siguiente comando:

Functions_general.py (Energy_profile) [path del archivo PWM o PRM] [PRM/PWM] [path del archivo fasta o raw target] [nombre del output]

#### Change Log

V4.2.3 (23/4)

- Se consolidan las funciones, se encapsulan un par de errores y donde solucionarlos
- La funcion de aliniamineto entrega un archivo ordenado por cantidad de mismatches
- Se modifica el algoritmo de fuzzy donde el treshold ahora reprenta el porcentaje entre que el nyucleotido no sea menos que aleatorio y bsu maximo para cada condición

V4.2.2 (10/4)

- Se modofica la funcion de energia para poder entregar un grafico

V4.2.1 (9/4)

- Se agrega la funcion para generar un perfil de energia a partir de un PWM y una secuencia target

V4.2.0 (5/4)

- Se cambian las funciones para entregar archivos en .tab
- Se agrega el score de la PRM y PWM
- Se agrega un alinimianto para la reversa de la secuencia target

V4.1.0 (2/4)

- Se agregan las funciones condensadas para cada funcion del script
- Se agrega la funcion para PWM
- Se agregan las funciones para el scoring de PWM y rel

V4.0.1 (27/3)

- Se agregó Función para aliniamiento
- Se modifico las funciones para fuzzy y para guardarlo como archivo en dos funciones separadas
- Se agregaron varias carpetas para modificar el flujo y logica del problema principal, esto para ordenar el acceso a las funciones para limitar la cantidad de imports necesarios
si no se van a utilzar ciertos modulos del script. Esto permitiria al usuario solo ver los requirements de los archivos que necesita. Solo se generó la arquitectura para hacer estos cambios todavia el flujo del progarama no ha cambiado.
- Se agrearon carpetas de input y output donde se deberian poner los files de input y autopmaticamente se crean los files en el output

V3.2 (12/3)

- Se agregó docstring a todos las funciones

V3.1 (10/3)  

- Modificaciónes de formato y de readme  

V3 (10/3)

- Se agregaron las funciones asociadas a la creacion de un archico tipo fuzzy  
- Se modificó el flujo del programa para poder aceptar las dos opciones, y poder recibir todo en una linea del terminal a diferencia de requerir inputs posteriores a correr el script.
- Queda a analisis de su uso, correción de posibles errores de programación y encapsulamiento de errores de usuario.  

V2 (8/3)  

- Se agrega la función para guardar el archivo como un txt  
- Se agregan la funcionalidad de correr el el script agregando el file a usar por consola, si bien ahora el programa tiene instrucciones posteriores a correrlo, todas estas opciones se pueden cambiar a argumentos por consola.  
- El programa esta completamente funcional, queda a revisión de su funcionalidad práctica.

V1.1 (7/3)  

- Se agrega la función que recibe la información de las otras funciones y escribe en un string en el formato PWM, tanto horizontal como vertical.

V1.0 (6/3)  

- Se crean las funciones principales para procesar la información.
- Abre un archivo (falta poder agregarle el argumento por consola). Siempre y cuando los archivos sean fasta, o solo contengan las lineas de nucleotidos, el script puede generar una lista con diccionarios que tienen como llave el nucleotido, y el conteo de ellos en cada posición.
- Falta generar una función para poder generar el PWM, y agrerar los argumentos por consola
