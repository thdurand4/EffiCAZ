import re
import sys
from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import click


@click.command(context_settings={'help_option_names': ('-h', '--help'), "max_content_width": 800})
@click.option('--predgpi_file', '-p', default=None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to output of predgpi tool')
@click.option('--predgpi_output', '-o', default=None,
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to create ID of sorted protein with no membrane anchor')

def main(predgpi_file, predgpi_output):
    """ This program retrieve ID of secreted protein of predGPI with no anchor membrane"""
    predgpi_id = []
    with open(predgpi_file) as f1:
        for lignes in f1:
            col = lignes.split("\t")
            if col[2] != "GPI-anchor":
                #print(col[0],col[2])
                predgpi_id.append(col[0])

    output_id = open(predgpi_output, "w")
    for elem in predgpi_id:
        output_id.write(elem+"\n")
    output_id.close()

if __name__ == '__main__':
    main()



