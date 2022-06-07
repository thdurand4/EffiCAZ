import re
import sys
from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import click


@click.command(context_settings={'help_option_names': ('-h', '--help'), "max_content_width": 800})
@click.option('--wolfpsort_file', '-in', default=None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to the wolfpsort file')
@click.option('--spwolfpsort_output', '-o', default=None,
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to output of wolfpsort parsed')
@click.option('--treshold', '-th', default=None,
              type=click.INT,
              required=True, show_default=True, help='treshold for score wolfpsort')

def main(wolfpsort_file, spwolfpsort_output, treshold):
    """ This program retrieve ID of secreted protein of wolfpsort with score cutoff wolfpsort """
    wolfpsort_id = []
    with open(wolfpsort_file) as f1:
        for lignes in f1:
            if re.search("^[^#]",lignes):
                good_lignes = re.sub(",","",lignes)
                col = good_lignes.split(" ")
                if col[1] == "extr" and int(col[2])>= treshold:
                    print(col[0],col[1],col[2])
                    wolfpsort_id.append(col[0])

    output_id = open(spwolfpsort_output, "w")
    for elem in wolfpsort_id:
        output_id.write(elem+"\n")
    output_id.close()

if __name__ == '__main__':
    main()
