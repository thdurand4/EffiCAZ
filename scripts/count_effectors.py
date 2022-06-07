import re
import sys
from collections import defaultdict
from Bio import SeqIO
import pandas as pd
from collections import Counter
import click

@click.command(context_settings={'help_option_names': ('-h', '--help'), "max_content_width": 800})
@click.option('--gff', '-g', default=None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to input gff file')
@click.option('--output', '-o', default=None,
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to output file')
@click.option('--fasta_file', '-fasta', default=None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to fasta file of effector')


def main(gff, output, fasta_file):
    """This count the number of effector per contig"""


    id_effector = []
    gff_parse = []
    final_effectors = []
    dico_gff = defaultdict(list)

    for record in SeqIO.parse(fasta_file, "fasta"):
        id_effector.append(record.id)



    with open(gff,"r") as f1 :
        for lignes in f1:
                ligne = lignes.rstrip("\n")
                col = ligne.split("\t")
                id_1 = re.sub("ID=", "", col[8])
                id_2 = re.sub(";Name=\w+", "", id_1)
                if col[2] =="gene":
                    gff_parse.append(col[0]+" "+col[1]+" "+col[2]+" "+id_2)
                    dico_gff[col[0]].append(id_2)




    #print(len(id_effector))


    for cle in dico_gff:
        for elem in dico_gff[cle]:
            for id in id_effector:
                if elem == id :
                    final_effectors.append(cle)


    counts = pd.Series(final_effectors).value_counts()
    counts.to_csv(output, header = False, sep="\t")
    #print(counts)

    #print(df)


if __name__ == '__main__':
    main()

