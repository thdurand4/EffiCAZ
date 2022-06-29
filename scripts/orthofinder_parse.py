import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import click


@click.command(context_settings={'help_option_names': ('-h', '--help'), "max_content_width": 800})
@click.option('--orthogroups_table', '-t', default=None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to input orthogroups table file')
@click.option('--strain_name', '-n', default=None,
              type=click.STRING, required=True)
@click.option('--output1', '-f', default=None,
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to output specific proteins sequences fasta file')
@click.option('--output2', '-c', default=None,
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to output specific orthogroups csv file')
@click.option('--path_seq', '-p', default=None,
              type=click.STRING, help='path to orthogroup sequences', required=True)

def main(orthogroups_table, strain_name, path_seq, output1, output2):
    gene_counts = pd.read_table(orthogroups_table, sep='\t')
    specific = gene_counts[gene_counts[strain_name] != 0]
    specific = specific[specific["Total"] == specific[strain_name]]

    file = ""
    list_record = []

    for name in specific["Orthogroup"]:
        file = str(path_seq) + str(name) + ".fa"
        for record in SeqIO.parse(file, "fasta"):
            list_record.append(record)

    SeqIO.write(list_record, str(output1), "fasta")
    specific.to_csv(str(output2))


if __name__ == '__main__':
    main()