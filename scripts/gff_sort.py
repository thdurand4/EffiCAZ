import pandas as pd
import click
import re
import os
import sys

@click.command(context_settings={'help_option_names': ('-h', '--help'), "max_content_width": 800})
@click.option('--gff', '-g', default=None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to input gff file')
@click.option('--output', '-o', default=None,
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to output file')
@click.option('--strain_name', '-name', default=None,
              type=click.STRING,
              required=True, show_default=True, help='Name of the strain')


def main(gff, output, strain_name):
    """This programme rename ID of the gff3 file"""
    gene_gff = []
    with open(gff, "r") as f1:
        for lignes in f1:
            col = lignes.split("\t")
            id_strain = re.sub("ID=","ID="+strain_name+"_",col[8])
            prot_gff = re.sub(";","T0;",id_strain)
            prot_gff = re.sub("T0T0", "T0", prot_gff)
            gene_gff.append(col[0]+"\t"+col[1]+"\t"+col[2]+"\t"+col[3]+"\t"+col[4]+"\t"+col[5]+"\t"+col[6]+"\t"+col[7]+"\t"+prot_gff)

    output_file = open(output,"w")
    for elem in gene_gff:
        output_file.write(elem)
    output_file.close()
if __name__ == '__main__':
    main()
