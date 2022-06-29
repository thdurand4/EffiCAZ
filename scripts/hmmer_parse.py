import re
import sys
from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import click
from Bio import SearchIO


@click.command(context_settings={'help_option_names': ('-h', '--help'), "max_content_width": 800})
@click.option('--hmmer_file', '-f', default=None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to output of hmmer file')
@click.option('--parse_output', '-o', default=None,
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to create output of hmmer results parsed')
def main(hmmer_file, parse_output):
    """ This program parse the output .tbl of hmmer """
    hits = defaultdict(list)
    test = set()
    attribs = ['id','bitscore','evalue']
    query = defaultdict(list)
    with open(hmmer_file) as f1:
        for lignes in f1:
            if re.search('^[^#]',lignes):
                ligne = lignes.split()
                query["query_name"].append(ligne[2])
                query["accession"].append(ligne[3])
                test.add(ligne[0])

    with open(hmmer_file) as f2:
        for record in SearchIO.parse(f2, 'hmmer3-tab'):
            for hit in record.hits:
                for attrib in attribs:
                    hits[attrib].append(getattr(hit, attrib))
    table = pd.DataFrame.from_dict(hits)
    table2 = pd.DataFrame.from_dict(query)
    frames = [table,table2]
    results = pd.concat(frames,axis=1)
    print(results)
    sort_table = (results.sort_values(by='evalue'))
    sort_table = sort_table.drop_duplicates(subset=["id"],keep='first')
    #print(sort_table)
    #print(len(test))
    sort_table.to_csv(parse_output,index=False)





if __name__ == '__main__':
    main()
