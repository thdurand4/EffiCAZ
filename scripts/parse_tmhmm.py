import re
import sys
from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import click

@click.command(context_settings={'help_option_names': ('-h', '--help'), "max_content_width": 800})
@click.option('--tmhmm_file', '-in', default=None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to input of tmhmm results')
@click.option('--parsetmhmm_output', '-o', default=None,
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to create secreted ID of sorted protein')
@click.option('--transmembranaire', '-tm', default=None,
              type=click.INT,
              required=True, show_default=True, help='Number of maxiumum of transmembranaire domain you want to keep')

def main(tmhmm_file, parsetmhmm_output, transmembranaire):
    """ This program parse the output of tool tmhmm with a cutoff and transmembranaire domain (max) """
    tm_parsing = []
    with open(tmhmm_file) as f1:
        for lignes in f1:
                ligne = lignes.rstrip("\n")
                col = ligne.split("\t")
                nb_tm = col[4].split("=")
                if  int(nb_tm[1]) <= transmembranaire:
                    tm_parsing.append(ligne)
    output_tmhmm = open(parsetmhmm_output, "w")
    for elem in tm_parsing:
        output_tmhmm.write(elem+"\n")
    output_tmhmm.close()

if __name__ == '__main__':
    main()

