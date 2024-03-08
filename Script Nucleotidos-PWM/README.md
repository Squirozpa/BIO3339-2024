# Script para generar un PWM a partir de datos de nucleotidos  

Script creado por Sebastian Quiroz para Pontificia Universidad Católica de Chile  
Contacto. <squirozpa@uc.cl>

Utilizacion:  
Correr el archivo principal Nuc-PWM.py en una carpeta con el archivo de nucleotidos a procesar, puede estar en formato fasta o raw en la forma: (pyhton) Nuc-PWM.py [nombre_archivo]
Responder si requiere que el output sea en horizontal o vertical  
Por último poner el nombre con el cuál se quiere guardar el archivo

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
