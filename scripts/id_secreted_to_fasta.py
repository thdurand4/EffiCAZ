import re
import sys
from collections import defaultdict
from Bio import SeqIO
import click


@click.command(context_settings={'help_option_names': ('-h', '--help'), "max_content_width": 800})
@click.option('--protein_file', '-fasta', default=None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to fasta protein sorted with tmhmm')
@click.option('--id_secreted_protein', '-id', default=None,
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to id of secreted protein (output of wolfpsort parsed)')
@click.option('--fasta_output', '-o', default=None,
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to fasta protein secreted output')

def main(protein_file, id_secreted_protein, fasta_output):
    """This programme use ID of secreted protein with signal peptide and protein fasta file to generate fasta protein with this ID"""
    id_secreted = []
    fasta_prot = defaultdict(str)

    with open(id_secreted_protein, "r") as f1:
        for lignes in f1:
            without_backspace = lignes.rstrip("\n")
            id_secreted.append(without_backspace)


    for record in SeqIO.parse(protein_file, "fasta"):
        fasta_prot[record.id] = record

    fasta_secretion = []

    for cle in fasta_prot:
        for elem in id_secreted:
            if elem == cle:
                fasta_secretion.append(fasta_prot[elem])

    SeqIO.write(fasta_secretion,fasta_output,"fasta")

if __name__ == '__main__':
    main()