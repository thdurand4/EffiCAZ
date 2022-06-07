import re
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
import click
import os

@click.command(context_settings={'help_option_names': ('-h', '--help'), "max_content_width": 800})
@click.option('--protein_file', '-p', default=None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to fasta protein')
@click.option('--fasta_output', '-o', default=None,
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to create fasta sorted protein')
@click.option('--strain_name', '-name', default=None,
              type=click.STRING,
              required=True, show_default=True, help='Name of the strain')



def main(protein_file, fasta_output, strain_name):
    """This programme remove * character on protein seq and rename it"""
    # read fasta and save to dict
    sorted_prot = []
    for record in SeqIO.parse(protein_file, "fasta"):
        record.id = str(strain_name)+"_"+str(record.id)
        no_stop = re.sub("\*","",str(record.seq))
        record.seq = Seq(no_stop)
        sorted_prot.append(record)
        record.description =""
        record.name = ""

    SeqIO.write(sorted_prot,fasta_output,"fasta")

if __name__ == '__main__':
    main()
