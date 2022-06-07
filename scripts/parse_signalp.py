import re
import sys
from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import click


@click.command(context_settings={'help_option_names': ('-h', '--help'), "max_content_width": 800})
@click.option('--signalp_file_gff', '-s', default=None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to output of gff3 of signalp')
@click.option('--spsignalp_output', '-o', default=None,
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to create secreted ID of sorted protein')
@click.option('--treshold', '-th', default=None,
              type=click.FloatRange(min=0, max=1, min_open=True, max_open=True),
              required=True, show_default=True, help='treshold for protein with signal peptide (between 0 & 1)')

def main(signalp_file_gff, spsignalp_output, treshold):
    """ This program retrieve ID of secreted protein of signalp with cutoff between 0 & 1 """
    spsignalp_id = []
    with open(signalp_file_gff) as f1:
        for lignes in f1:
            if re.search("^[^#]",lignes):
                col = lignes.split("\t")
                if float(col[5]) >= treshold:
                    #print(col[0],col[5])
                    spsignalp_id.append(col[0])

    output_id = open(spsignalp_output, "w")
    for elem in spsignalp_id:
        output_id.write(elem+"\n")
    output_id.close()

if __name__ == '__main__':
    main()



