# Script para generar un PWM a partir de datos de nucleotidos  

Script creado por Sebastian Quiroz para Pontificia Universidad Católica de Chile  
Contacto. <squirozpa@uc.cl>

Utilizacion:  
Correr el archivo principal Nuc-PWM.py en una carpeta con el archivo de nucleotidos a procesar, puede estar en formato fasta o raw en la forma:  
El script recibe 4 argumentos  

1. Primero el nombre del arhivo de nucleotidos tanto en fasta como raw  
2. Segundo el tipo de analisis que se busca hacer: -f para fuzzy y -p para PWM  
3. Tercero, dependiendo de si se requiere fuzzy o PWM. Para fuzzy se requiere entregar el valor the threshold para el analisis en formato de frequencia relativa (0-0.99)
si es PWM se requiere entregar la orientacion del output file -h para horizontal y -v para vertical  
4. Cuarto el nombre del arhico a guardar, se entrega en un formato txt, no agregar extensión ni espacios.  

Ejemplos:  

- Fuzzy con threshold de 0.6  
(pyhton) Nuc-PWM.py [nombre_archivo.fasta] -f 0.6 [output]  

- PWM con orientacion vertical  
(python3) Nuc_PWM.py [nombre_arhivo.raw] -p -h [output]  

Change Log  
V1.0 (6/3)  
Se crean las funciones principales para procesar la información.  
Abre un archivo (falta poder agregarle el argumento por consola). Siempre y cuando los archivos sean fasta, o solo contengan las lineas de nucleotidos, el script puede generar una lista con diccionarios que tienen como llave el nucleotido, y el conteo de ellos en cada posición.
Falta generar una función para poder generar el PWM, y agrerar los argumentos por consola

V1.1 (7/3)  
Se agrega la función que recibe la información de las otras funciones y escribe en un string en el formato PWM, tanto horizontal como vertical.

V2 (8/3)  
Se agrega la función para guardar el archivo como un txt  
Se agregan la funcionalidad de correr el el script agregando el file a usar por consola, si bien ahora el programa tiene instrucciones posteriores a correrlo, todas estas opciones se pueden cambiar a argumentos por consola.  
El programa esta completamente funcional, queda a revisión de su funcionalidad práctica.

V3 (10/3)
Se agregaron las funciones asociadas a la creacion de un archico tipo fuzzy  
Se modificó el flujo del programa para poder aceptar las dos opciones, y poder recibir todo en una linea del terminal a diferencia de requerir inputs posteriores a correr el script.
Queda a analisis de su uso, correción de posibles errores de programación y encapsulamiento de errores de usuario.  

V3.1 (10/3)  
Modificaciónes de formato y de readme  

V3.2 (12/3)
Se agregó docstring a todos las funciones
