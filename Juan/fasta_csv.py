import csv

def csv_to_fasta(input_file, output_file):
    with open(input_file, 'r') as csv_file, open(output_file, 'w') as fasta_file:
        csv_reader = csv.reader(csv_file)
        for row in csv_reader:
            sequence_id = row[0]
            sequence = row[1]
            fasta_file.write(f'>{sequence_id}\n{sequence}\n')

# Cambia 'input.csv' al nombre de tu archivo CSV y 'output.fasta' al nombre que desees para el archivo FASTA de salida
csv_to_fasta('data_prot.csv', 'sequence.fasta')
