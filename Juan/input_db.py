import csv

# Archivos de entrada y salida
archivo_a = 'data_A.csv'
archivo_b = 'data_B.csv'
archivo_resultado = 'input_DD.csv'
nombres_columnas = ['complex_name', 'chainA', 'chainB','SMILES']
path = '/home/ramon/juan/alphafold/working_structures/'

# Leer archivos de entrada y escribir las columnas 2 y 3 en el archivo de salida
with open(archivo_a, 'r') as file_a, open(archivo_b, 'r') as file_b, open(archivo_resultado, 'w', newline='') as file_resultado:
    reader_a = csv.reader(file_a)
    reader_b = csv.reader(file_b)
    writer = csv.writer(file_resultado)

    writer.writerow(nombres_columnas)
    next(reader_a)
    next(reader_b)

    # Iterar sobre las filas de ambos archivos simultáneamente y escribir en el archivo de salida
    for row_a, row_b in zip(reader_a, reader_b):
        # Combinar las columnas 2 y 3 de ambos archivos
        row = [row_a[1] + '_' + row_b[1] + '_'+ row_a[6], path + row_a[1] + '.pdb', path + row_b[1] + '.pdb', row_b[3]]
        # Escribir la fila combinada en el archivo de salida
        writer.writerow(row)
