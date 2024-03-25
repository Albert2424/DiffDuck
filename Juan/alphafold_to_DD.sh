#!/bin/bash
mkdir deepfold_work_structures
cp structure/*rank_001*.pdb deepfold_work_structures

# Iterar sobre los archivos en el directorio deepfold_work_structures/
for file in deepfold_work_structures/*; do
    if [[ -f "$file" ]]; then  # Verificar si el elemento es un archivo regular
        # Obtener el nombre del archivo sin la ruta
        filename=$(basename "$file")
        # Extraer los primeros dígitos después del guion bajo "_"
        new_filename=$(echo "$filename" | cut -d '_' -f 2) #esta linea saca el numero del nombre
        # Mover el archivo con el nuevo nombre
        mv "$file" "deepfold_work_structures/$new_filename.pdb" 
    fi
done
