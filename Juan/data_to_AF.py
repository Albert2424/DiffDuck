import csv

# Nombre del archivo de entrada y las columnas que deseas extraer
archivo_entrada = 'clean_data.csv'
columnas_a_extraer = [0,6,1]  # Índices de las columnas que deseas extraer (0-indexed)

# Nombre del archivo de salida
archivo_salida = 'AF_input.csv'

# Leer el archivo de entrada y escribir las columnas seleccionadas en el archivo de salida
with open(archivo_entrada, 'r', newline='') as f_input, open(archivo_salida, 'w', newline='') as f_output:
    csv_reader = csv.reader(f_input)
    csv_writer = csv.writer(f_output)

    for row in csv_reader:
        # Seleccionar solo las columnas especificadas
        nueva_fila = [row[i] for i in columnas_a_extraer]
        nueva_fila[1] = "'" + nueva_fila[1] +"'"
        # Escribir la fila en el archivo de salida
        csv_writer.writerow(nueva_fila)

print("Se han extraído y guardado las columnas seleccionadas en", archivo_salida)
