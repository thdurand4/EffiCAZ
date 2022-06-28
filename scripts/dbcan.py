import pandas as pd
import click



@click.command(context_settings={'help_option_names': ('-h', '--help'), "max_content_width": 800})
@click.option('--dbcan', '-i', default=None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to input dbcan result file')
@click.option('--gff', '-g', default=None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to input gff file')
@click.option('--output', '-o', default=None,
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to output count csv file')

def main(dbcan, gff, output):

    dbcan = pd.read_table(dbcan)
    dbcan.columns = ["gene", "HMMER", "Hotpep", "Diamond", "#ofTools"]
    dbcan = dbcan[dbcan["#ofTools"] == 3]

    gff = pd.read_table(gff)
    gff.columns = ["contig", "source", "feature", "start", "stop", "score_1", "strand", "score_2", "id"]
    gff = gff[gff["feature"] == "gene"]
    gff["id"] = gff["id"].str.replace(r'Name=', "")
    gff["id"] = gff["id"].str.replace(r'Parent=', "")
    gff["id"] = gff["id"].str.replace(r'\;[0-9a-zA-Z_]*', "")
    gff["id"] = gff["id"].str.replace(r'ID=', "")
    print(gff.head())
    gff = gff[gff["id"].isin(dbcan["gene"])]
    #gff.reset_index(inplace=True)

    dbcan = pd.concat([dbcan, gff], axis=1)
    count = dbcan["contig"].value_counts()
    count.to_csv(str(output))

if __name__ == '__main__':
    main()