import re
import sys
from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import click


@click.command(context_settings={'help_option_names': ('-h', '--help'), "max_content_width": 800})
@click.option('--phobius_file', '-t', default=None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to output of phobius file')
@click.option('--spphobius_output', '-o', default=None,
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to create output of secreted ID of sorted protein')

def main(phobius_file, spphobius_output):
    """ This program retrieve ID of secreted protein of phobius with TM == 0 or TM == 1 """
    sptargetp_id = []
    with open(phobius_file) as f1:
        for lignes in f1:
            if re.search("^[^SEQUENCE]",lignes):
                col = lignes.split()
                if col[1] == "0" and col[2] == "Y":
                    #print(col[0],col[1],col[2])
                    sptargetp_id.append(col[0])
                if col[1] == "1" and col[2] == "Y":
                    #print(col[0], col[1], col[2])
                    sptargetp_id.append(col[0])

    output_id = open(spphobius_output, "w")
    for elem in sptargetp_id:
        output_id.write(elem+"\n")
    output_id.close()

if __name__ == '__main__':
    main()



