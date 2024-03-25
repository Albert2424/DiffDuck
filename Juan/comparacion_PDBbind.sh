#!/bin/bash

#una vez descargado los resultados del blastp en csv format y los index de la base de datos PDBbind
blastp=$1

#me quedo con los resultados con identity de 100%
grep '100.000' $blastp > identity.csv

#me quedo solo con el id y el nombre del pdb y le quito si tiene alguna terminacion el PDB
cut -d ',' -f 1,2 identity.csv > columnas.csv
sed 's/_.$//' columnas.csv > blastp_pdb.csv

awk -F ',' '{print $2}' blastp_pdb.csv | sort > palabras_csv.txt
sort pdbbind/index/2020_index.lst > palabras_lst.txt
cut -d ' ' -f1 palabras_lst.txt > lst.txt
tr '[:upper:]' '[:lower:]' < palabras_csv.txt > csv.txt

join csv.txt lst.txt | awk -F ',' '{print $0}' > coincidencias.csv


# Crear un archivo temporal para almacenar las líneas coincidentes
touch match.csv

# Leer cada línea del archivo de palabras y buscar coincidencias en el archivo CSV
while IFS= read -r palabra; do
    grep -i "$palabra" identity.csv >> match.csv
done < coincidencias.csv

cat match.csv | cut -d ',' -f1 |sort -n |uniq >> pdbbind_match.csv

rm identity.csv columnas.csv blastp_pdb.csv palabras_lst.txt palabras_csv.txt lst.txt csv.txt coincidencias.csv match.csv 