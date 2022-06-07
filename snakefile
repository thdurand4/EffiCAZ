configfile: "config/config.yaml"

fasta_prot_dir = config["DATA"]["PROTEIN"]
output_dir = config["DATA"]["OUTPUT"]
log_dir = f"{output_dir}LOGS/"
script_dir = config["DATA"]["SCRIPTS"]
gff_dir = config["DATA"]["GFF"]

PROTEIN, = glob_wildcards(fasta_prot_dir+"{samples}.fasta", followlinks=True)




def get_threads(rule, default):
    """
    give threads or 'cpus-per-task from cluster_config rule : threads to SGE and cpus-per-task to SLURM
    """
    if cluster_config:
        if rule in cluster_config and 'threads' in cluster_config[rule]:
            return int(cluster_config[rule]['threads'])
        elif rule in cluster_config and 'cpus-per-task' in cluster_config[rule]:
            return int(cluster_config[rule]['cpus-per-task'])
        elif '__default__' in cluster_config and 'cpus-per-task' in cluster_config['__default__']:
            return int(cluster_config['__default__']['cpus-per-task'])
        elif '__default__' in cluster_config and 'threads' in cluster_config['__default__']:
            return int(cluster_config['__default__']['threads'])
    if workflow.global_resources["_cores"]:
        return workflow.global_resources["_cores"]
    return default



rule finale:
    input:
        domain_prot = expand(f"{output_dir}3_HMMER_PFAM/{{samples}}_secreted.tbl", samples = PROTEIN),
        effector_contig = expand(f"{output_dir}5_FINAL_RESULT/EFFECTOR/{{samples}}/{{samples}}_effector_per_contig.txt", samples = PROTEIN)


rule rename_protein:
    threads: get_threads("rename_protein",1)
    input:
        protein = f"{fasta_prot_dir}{{samples}}.fasta"
    output:
        sorted_protein = f"{output_dir}1_PROTEIN_SORTED/{{samples}}.fasta"
    log :
        error =  f'{log_dir}protein_sorted/protein_sorted_{{samples}}.e',
        output = f'{log_dir}protein_sorted/protein_sorted_{{samples}}.o'
    message:
            f"""
             Running {{rule}}
                Input:
                    - Fasta : {{input.protein}}
                Output:
                    - Protein_sorted: {{output.sorted_protein}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}

            """
    shell:
        f"python {script_dir}rename_prot.py -p {{input.protein}} -o {{output.sorted_protein}} -name {{wildcards.samples}} 1>{{log.output}} 2>{{log.error}}"


rule phobius:
    threads: get_threads("phobius",5)
    input:
        protein = rules.rename_protein.output.sorted_protein
    output:
        output_phobius = f"{output_dir}2_SECRETED_PROTEIN/PHOBIUS/{{samples}}/{{samples}}_phobius.tsv"
    log :
        error =  f'{log_dir}phobius/phobius_{{samples}}.e',
        output = f'{log_dir}phobius/phobius_{{samples}}.o'
    message:
            f"""
             Running {{rule}}
                Input:
                    - Fasta : {{input.protein}}
                Output:
                    - Phobius_TSV: {{output.output_phobius}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}

            """
    envmodules:
            "phobius_local"
    shell:
        f"phobius.pl -short {{input.protein}} 1>{{output.output_phobius}} 2>{{log.error}}"

rule signalP:
    threads: get_threads("signalP",10)
    input:
        protein = rules.rename_protein.output.sorted_protein
    output:
        output_signalP = f"{output_dir}2_SECRETED_PROTEIN/SignalP/{{samples}}/output.gff3"
    params:
        output_dir_phobius =f"{output_dir}2_SECRETED_PROTEIN/SignalP/{{samples}}"
    log :
        error =  f'{log_dir}signalP/phobius_{{samples}}.e',
        output = f'{log_dir}signalP/phobius_{{samples}}.o'
    message:
            f"""
             Running {{rule}}
                Input:
                    - Fasta : {{input.protein}}
                Output:
                    - Phobius_TSV: {{output.output_signalP}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}

            """
    envmodules:
            "signalp/"
    shell:
        f"signalp6 --fastafile {{input.protein}} --output_dir {{params.output_dir_phobius}} --organism eukarya 1>{{log.output}} 2>{{log.error}}"
        f"\nrm -rf {{output_dir}}2_SECRETED_PROTEIN/SignalP/{{wildcards.samples}}/output_*.txt"

rule targetp:
    threads: get_threads("targetp",10)
    input:
        protein = rules.rename_protein.output.sorted_protein
    output:
        output_targetp = f"{output_dir}2_SECRETED_PROTEIN/TargetP/{{samples}}/{{samples}}_summary.targetp2"
    log :
        error =  f'{log_dir}targetp/targetp_{{samples}}.e',
        output = f'{log_dir}targetp/targetp_{{samples}}.o'
    message:
            f"""
             Running {{rule}}
                Input:
                    - Protein : {{input.protein}}
                Output:
                    - TargetP_summary: {{output.output_targetp}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}

            """
    envmodules:
            "targetp_local"
    shell:
        f"targetp -fasta {{input.protein}} -stdout 1>{{output.output_targetp}} 2>{{log.error}}"

rule predgpi:
    threads: get_threads("predgpi",10)
    input:
        protein = rules.rename_protein.output.sorted_protein
    output:
        output_predgpi = f"{output_dir}2_SECRETED_PROTEIN/PredGPI/{{samples}}/{{samples}}.predgpi"
    log:
        error =  f'{log_dir}predgpi/predgpi_{{samples}}.e',
        output = f'{log_dir}predgpi/predgpi_{{samples}}.o'
    message:
        f"""
                 Running {{rule}}
                    Input:
                        - Protein : {{input.protein}}
                    Output:
                        - Predgpi_summary: {{output.output_predgpi}}
                    Others
                        - Threads : {{threads}}
                        - LOG error: {{log.error}}
                        - LOG output: {{log.output}}

                """
    envmodules:
        "predgpi_local"
    shell:
        f"predgpi.py -f {{input.protein}} -o {{output.output_predgpi}} 1>{{log.output}} 2>{{log.error}}"

rule parse_phobius:
    threads: get_threads("parse_phobius",1)
    input:
        result_phobius = rules.phobius.output.output_phobius
    output:
        secreted_phobius = f"{output_dir}2_SECRETED_PROTEIN/ID_SECRETED/{{samples}}/PHOBIUS/{{samples}}_secretedID.phobius"
    log:
        error =  f'{log_dir}parse_phobius/parse_phobius_{{samples}}.e',
        output = f'{log_dir}parse_phobius/parse_phobius_{{samples}}.o'
    message:
        f"""
                 Running {{rule}}
                    Input:
                        - Results_phobius : {{input.result_phobius}}
                    Output:
                        - Parse_phobius: {{output.secreted_phobius}}
                    Others
                        - Threads : {{threads}}
                        - LOG error: {{log.error}}
                        - LOG output: {{log.output}}

                """
    shell:
        f"python {script_dir}parse_phobius.py -t {{input.result_phobius}} -o {{output.secreted_phobius}} 1>{{log.output}} 2>{{log.error}}"

rule parse_signalp:
    threads: get_threads("parse_signalp",1)
    input:
        result_signalp = rules.signalP.output.output_signalP
    output:
        secreted_signalp = f"{output_dir}2_SECRETED_PROTEIN/ID_SECRETED/{{samples}}/SIGNALP/{{samples}}_secretedID.signalp"
    params:
        threshold=config["TOOLS_PARAMS"]["PARSE_SIGNALP_TRESHOLD"]
    log:
        error =  f'{log_dir}parse_signalp/parse_signalp_{{samples}}.e',
        output = f'{log_dir}parse_signalp/parse_signalp_{{samples}}.o'
    message:
        f"""
                 Running {{rule}}
                    Input:
                        - Results_signalp : {{input.result_signalp}}
                    Output:
                        - Parse_signalp: {{output.secreted_signalp}}
                    Others
                        - Threads : {{threads}}
                        - LOG error: {{log.error}}
                        - LOG output: {{log.output}}

                """
    shell:
        f"python {script_dir}parse_signalp.py -s {{input.result_signalp}} -o {{output.secreted_signalp}} -th {{params.threshold}} 1>{{log.output}} 2>{{log.error}}"

rule parse_targetp:
    threads: get_threads("parse_targetp",1)
    input:
        result_targetp = rules.targetp.output.output_targetp
    output:
        secreted_targetp = f"{output_dir}2_SECRETED_PROTEIN/ID_SECRETED/{{samples}}/TARGETP/{{samples}}_secretedID.targetp"
    params:
        threshold=config["TOOLS_PARAMS"]["PARSE_TARGETP_TRESHOLD"]
    log:
        error =  f'{log_dir}parse_targetp/parse_targetp_{{samples}}.e',
        output = f'{log_dir}parse_targetp/parse_targetp_{{samples}}.o'
    message:
        f"""
                 Running {{rule}}
                    Input:
                        - Results_targetp : {{input.result_targetp}}
                    Output:
                        - Parse_targetp: {{output.secreted_targetp}}
                    Others
                        - Threads : {{threads}}
                        - LOG error: {{log.error}}
                        - LOG output: {{log.output}}

                """
    shell:
        f"python {script_dir}parse_targetp.py -t {{input.result_targetp}} -o {{output.secreted_targetp}} -th {{params.threshold}} 1>{{log.output}} 2>{{log.error}}"


rule parse_predgpi:
    threads: get_threads("parse_predgpi",1)
    input:
        result_predgpi = rules.predgpi.output.output_predgpi
    output:
        noanchor_predgpi = f"{output_dir}2_SECRETED_PROTEIN/ID_SECRETED/{{samples}}/PREDGPI/{{samples}}_noanchorID.predgpi"
    params:
        threshold=config["TOOLS_PARAMS"]["PARSE_TARGETP_TRESHOLD"]
    log:
        error =  f'{log_dir}parse_predgpi/parse_predgpi_{{samples}}.e',
        output = f'{log_dir}parse_predgpi/parse_predgpi_{{samples}}.o'
    message:
        f"""
                 Running {{rule}}
                    Input:
                        - Result_predgpi : {{input.result_predgpi}}
                    Output:
                        - Parse_predgpi: {{output.noanchor_predgpi}}
                    Others
                        - Threads : {{threads}}
                        - LOG error: {{log.error}}
                        - LOG output: {{log.output}}

                """
    shell:
        f"python {script_dir}parse_predgpi.py -p {{input.result_predgpi}} -o {{output.noanchor_predgpi}} 1>{{log.output}} 2>{{log.error}}"

rule intersect_tools:
    threads: get_threads("intersect_tools",1)
    input:
        result_parse_predgpi = rules.parse_predgpi.output.noanchor_predgpi,
        result_parse_targetp = rules.parse_targetp.output.secreted_targetp,
        result_parse_signalp = rules.parse_signalp.output.secreted_signalp,
        result_parse_phobius = rules.parse_phobius.output.secreted_phobius
    output:
        signalpeptide_id = f"{output_dir}2_SECRETED_PROTEIN/ID_SECRETED/{{samples}}/{{samples}}_intersect.signalpeptide"
    log:
        error=f'{log_dir}intersect_tools/intersect_tools_{{samples}}.e',
        output=f'{log_dir}intersect_tools/intersect_tools_{{samples}}.o'
    message:
        f"""
                     Running {{rule}}
                        Input:
                            - Parse_phobius : {{input.result_parse_phobius}}
                            - Parse_targetp : {{input.result_parse_targetp}}
                            - Parse_signalp : {{input.result_parse_signalp}}
                            - Parse_predgpi : {{input.result_parse_predgpi}}
                        Output:
                            - Intersect_signalpeptide: {{output.signalpeptide_id}}
                        Others
                            - Threads : {{threads}}
                            - LOG error: {{log.error}}
                            - LOG output: {{log.output}}

                    """
    shell:
        f"bash {script_dir}intersect.sh {{input.result_parse_predgpi}} {{input.result_parse_signalp}} {{input.result_parse_phobius}} {{input.result_parse_targetp}} 1>{{output.signalpeptide_id}} 2>{{log.error}}"


rule fasta_intersect:
    threads: get_threads("fasta_intersect",1)
    input:
        fasta_protein = rules.rename_protein.output.sorted_protein,
        intersect_id = rules.intersect_tools.output.signalpeptide_id
    output:
        intersect_fasta_prot = f"{output_dir}2_SECRETED_PROTEIN/ID_SECRETED/{{samples}}/{{samples}}_intersect.fasta"
    log:
        error=f'{log_dir}fasta_intersect/fasta_intersect_{{samples}}.e',
        output=f'{log_dir}fasta_intersect/fasta_intersect_{{samples}}.o'
    message:
        f"""
                         Running {{rule}}
                            Input:
                                - Fasta_protein : {{input.fasta_protein}}
                                - Intersect_id : {{input.intersect_id}}
                            Output:
                                - Fasta_intersect_protein : {{output.intersect_fasta_prot}}
                            Others
                                - Threads : {{threads}}
                                - LOG error: {{log.error}}
                                - LOG output: {{log.output}}

                        """
    shell:
        f"python {script_dir}fasta_intersect.py -p {{input.fasta_protein}} -s {{input.intersect_id}} -o {{output.intersect_fasta_prot}}  1>{{log.output}} 2>{{log.error}}"

rule tmhmm:
    threads: get_threads("tmhmm",5)
    input:
        fasta_intersect_prot = rules.fasta_intersect.output.intersect_fasta_prot
    output:
        tmhmm_output = f"{output_dir}2_SECRETED_PROTEIN/TMHMM/{{samples}}/{{samples}}.tmhmm"
    log:
        error=f'{log_dir}tmhmm/tmhmm_{{samples}}.e',
        output=f'{log_dir}tmhmm/tmhmm_{{samples}}.o'
    message:
        f"""
                             Running {{rule}}
                                Input:
                                    - Fasta_protein_intersect : {{input.fasta_intersect_prot}}
                                Output:
                                    - TMHMM output : {{output.tmhmm_output}}
                                Others
                                    - Threads : {{threads}}
                                    - LOG error: {{log.error}}
                                    - LOG output: {{log.output}}

                            """
    envmodules:
        "tmhmm/2.0c"
    shell:
        f"tmhmm -short {{input.fasta_intersect_prot}} 1>{{output.tmhmm_output}} 2>{{log.error}}"

rule parse_tmhmm:
    threads: get_threads("parse_tmhmm",1)
    input:
        tmhmm_outfile = rules.tmhmm.output.tmhmm_output
    output:
        tmhmm_parsed_file = f"{output_dir}2_SECRETED_PROTEIN/TMHMM/{{samples}}/{{samples}}_tmhmm_parsed.tsv"
    params:
        threshold=config["TOOLS_PARAMS"]["PARSE_TMHMM_TRESHOLD"]
    log:
        error=f'{log_dir}parse_tmhmm/parse_tmhmm_{{samples}}.e',
        output=f'{log_dir}parse_tmhmm/parse_tmhmm_{{samples}}.o'
    message:
        f"""
                                 Running {{rule}}
                                    Input:
                                        - TMHMM output : {{input.tmhmm_outfile}}
                                    Output:
                                        - TMHMM parsed : {{output.tmhmm_parsed_file}}
                                    Others
                                        - Threads : {{threads}}
                                        - LOG error: {{log.error}}
                                        - LOG output: {{log.output}}

                                """
    shell:
        f"python {script_dir}parse_tmhmm.py -in {{input.tmhmm_outfile}} -tm {{params.threshold}} -o {{output.tmhmm_parsed_file}}  1>{{log.output}} 2>{{log.error}}"

rule tmhmm_fasta:
    threads: get_threads("tmhmm_fasta",1)
    input:
        tmhmm_parsed = rules.parse_tmhmm.output.tmhmm_parsed_file,
        protein_intersected = rules.fasta_intersect.output.intersect_fasta_prot
    output:
        fasta_parsed = f"{output_dir}2_SECRETED_PROTEIN/TMHMM/{{samples}}/{{samples}}_tmhmm_parsed.fasta"
    log:
        error=f'{log_dir}tmhmm_fasta/tmhmm_fasta_{{samples}}.e',
        output=f'{log_dir}tmhmm_fasta/tmhmm_fasta{{samples}}.o'
    message:
        f"""
                                 Running {{rule}}
                                    Input:
                                        - TMHMM parsed : {{input.tmhmm_parsed}}
                                        - Fasta intersect : {{input.protein_intersected}}
                                    Output:
                                        - Fasta protein : {{output.fasta_parsed}}
                                    Others
                                        - Threads : {{threads}}
                                        - LOG error: {{log.error}}
                                        - LOG output: {{log.output}}

                                """
    shell:
        f"python {script_dir}tmhmm_to_fasta.py -p {{input.protein_intersected}} -tmhmm {{input.tmhmm_parsed}} -o {{output.fasta_parsed}}  1>{{log.output}} 2>{{log.error}}"

rule wolfpsort:
    threads: get_threads("wolfpsort",10)
    input:
        protein_tmhmm = rules.tmhmm_fasta.output.fasta_parsed
    output:
        result_wolfpsort = f"{output_dir}2_SECRETED_PROTEIN/WOLFPSORT/{{samples}}/{{samples}}_wolfpsort.txt"
    log:
        error=f'{log_dir}wolfpsort/wolfpsort_{{samples}}.e',
        output=f'{log_dir}wolfpsort/wolfpsort{{samples}}.o'
    message:
        f"""
                                    Running {{rule}}
                                       Input:
                                           - Fasta TMHMM : {{input.protein_tmhmm}}
                                       Output:
                                           - Result WOLFPSORT : {{output.result_wolfpsort}}
                                       Others
                                           - Threads : {{threads}}
                                           - LOG error: {{log.error}}
                                           - LOG output: {{log.output}}

                                   """
    envmodules:
        "wolfpsort/0.2"
    shell:
        f"runWolfPsortSummary fungi < {{input.protein_tmhmm}} 1>{{output.result_wolfpsort}} 2>{{log.error}}"
rule parse_wolfpsort:
    threads: get_threads("parse_wolfpsort",1)
    input:
        wolfpsort_output = rules.wolfpsort.output.result_wolfpsort
    output:
        id_secreted_prot = f"{output_dir}5_FINAL_RESULT/SECRETED_PROTEIN/{{samples}}/{{samples}}_secreted.id"
    params:
        threshold=config["TOOLS_PARAMS"]["PARSE_WOLFPOSORT_TRESHOLD"]
    log:
        error=f'{log_dir}parse_wolfpsort/parse_wolfpsort_{{samples}}.e',
        output=f'{log_dir}parse_wolfpsort/parse_wolfpsort{{samples}}.o'
    message:
        f"""
                                        Running {{rule}}
                                           Input:
                                               - WOLFPSORT OUTPUT : {{input.wolfpsort_output}}
                                           Output:
                                               - ID SECRETED PROTEIN : {{output.id_secreted_prot}}
                                           Others
                                               - Threads : {{threads}}
                                               - LOG error: {{log.error}}
                                               - LOG output: {{log.output}}

                                       """

    shell:
        f"python {script_dir}parse_wolfpsort.py -in {{input.wolfpsort_output}} -th {{params.threshold}} -o {{output.id_secreted_prot}}  1>{{log.output}} 2>{{log.error}}"

rule id_tofasta_secreted :
    threads: get_threads("id_tofasta_secreted",1)
    input:
        id_secreted = rules.parse_wolfpsort.output.id_secreted_prot,
        fasta_tmhmm = rules.tmhmm_fasta.output.fasta_parsed
    output:
        fasta_prot_secreted = f"{output_dir}5_FINAL_RESULT/SECRETED_PROTEIN/{{samples}}/{{samples}}_secreted.fasta"
    log:
        error=f'{log_dir}parse_wolfpsort/parse_wolfpsort_{{samples}}.e',
        output=f'{log_dir}parse_wolfpsort/parse_wolfpsort{{samples}}.o'
    message:
        f"""
                                            Running {{rule}}
                                               Input:
                                                   - FASTA PROTEIN : {{input.fasta_tmhmm}}
                                                   - ID SECRETED : {{input.id_secreted}}
                                               Output:
                                                   - FASTA SECRETED PROTEIN : {{output.fasta_prot_secreted}}
                                               Others
                                                   - Threads : {{threads}}
                                                   - LOG error: {{log.error}}
                                                   - LOG output: {{log.output}}

                                           """
    shell:
        f"python {script_dir}id_secreted_to_fasta.py -fasta {{input.fasta_tmhmm}} -id {{input.id_secreted}} -o {{output.fasta_prot_secreted}} 1>{{log.output}} 2>{{log.error}}"

rule hmmer_pfam :
    threads: get_threads("hmmer_pfam", 8)
    input:
        fasta_secreted = rules.id_tofasta_secreted.output.fasta_prot_secreted,
        bdd_pfam = config["DATA"]["BDD_PFAM"]
    output:
        protein_secreted_domain = f"{output_dir}3_HMMER_PFAM/{{samples}}_secreted.tbl"
    params:
        param_hmmer = config["TOOLS_PARAMS"]["HMMER"]
    log:
        error=f'{log_dir}hmmer_pfam/hmmer_pfam_{{samples}}.e',
        output=f'{log_dir}hmmer_pfam/hmmer_pfam_{{samples}}.o'
    message:
        f"""
                                                Running {{rule}}
                                                   Input:
                                                       - FASTA PROTEIN : {{input.fasta_secreted}}
                                                       - BDD PFAM : {{input.bdd_pfam}}
                                                   Output:
                                                       - DOMMAINES PROTEINES : {{output.protein_secreted_domain}}
                                                   Others
                                                       - Threads : {{threads}}
                                                       - LOG error: {{log.error}}
                                                       - LOG output: {{log.output}}

                                               """
    envmodules:
        "hmmer/3.2.1"
    shell:
        f"hmmsearch --tblout {{output.protein_secreted_domain}} {{params.param_hmmer}} {{input.bdd_pfam}} {{input.fasta_secreted}} 1>{{log.output}} 2>{{log.error}}"

rule effectorP :
    threads: get_threads("effectorP", 10)
    input:
        fasta_secreted = rules.id_tofasta_secreted.output.fasta_prot_secreted
    output:
        fasta_effectors = f"{output_dir}5_FINAL_RESULT/EFFECTOR/{{samples}}/{{samples}}_effector.fasta",
        effectorP_out = f"{output_dir}5_FINAL_RESULT/EFFECTOR/{{samples}}/{{samples}}_effectorP.out",
        no_effector_fasta = f"{output_dir}5_FINAL_RESULT/EFFECTOR/{{samples}}/{{samples}}_non_effector.fasta"
    log:
        error=f'{log_dir}effectorP/effectorP_{{samples}}.e',
        output=f'{log_dir}effectorP/effectorP_{{samples}}.o'
    message:
        f"""
                                                    Running {{rule}}
                                                       Input:
                                                           - FASTA PROTEIN : {{input.fasta_secreted}}
                                                       Output:
                                                           - EFFECTOR FASTA : {{output.fasta_effectors}}
                                                           - NON EFFECTOR FASTA : {{output.no_effector_fasta}}
                                                           - EFFECTORP_OUT : {{output.effectorP_out}}
                                                       Others
                                                           - Threads : {{threads}}
                                                           - LOG error: {{log.error}}
                                                           - LOG output: {{log.output}}

                                                   """
    envmodules:
        "effectorp_local"
    shell:
        f"EffectorP.py -o {{output.effectorP_out}} -E {{output.fasta_effectors}} -N {{output.no_effector_fasta}} -i {{input.fasta_secreted}} 1>{{log.output}} 2>{{log.error}}"

rule sort_gff:
    threads: get_threads("sort_gff",1)
    input:
        gff_file = f"{gff_dir}{{samples}}.gff3"
    output:
        gff_sorted = f"{output_dir}4_GFF_SORTED/{{samples}}/{{samples}}_sorted.gff3",
    log:
        error=f'{log_dir}sort_gff/sort_gff_{{samples}}.e',
        output=f'{log_dir}sort_gff/sort_gff_{{samples}}.o'
    message:
        f"""
                                                        Running {{rule}}
                                                           Input:
                                                               - GFF FILE : {{input.gff_file}}
                                                           Output:
                                                               - GFF SORTED : {{output.gff_sorted}}
                                                           Others
                                                               - Threads : {{threads}}
                                                               - LOG error: {{log.error}}
                                                               - LOG output: {{log.output}}

                                                       """
    shell:
        f"python {script_dir}gff_sort.py -g {{input.gff_file}} -o {{output.gff_sorted}} -name {{wildcards.samples}} 1>{{log.output}} 2>{{log.error}}"

rule count_effector:
    threads: get_threads("count_effector", 1)
    input:
        fasta_effectors = rules.effectorP.output.fasta_effectors,
        gff_protein = rules.sort_gff.output.gff_sorted
    output:
        effector_per_contig = f"{output_dir}5_FINAL_RESULT/EFFECTOR/{{samples}}/{{samples}}_effector_per_contig.txt"
    log:
        error=f'{log_dir}count_effector/count_effector_{{samples}}.e',
        output=f'{log_dir}count_effector/count_effector_{{samples}}.o'
    message:
        f"""
                                                        Running {{rule}}
                                                           Input:
                                                               - FASTA EFFECTOR : {{input.fasta_effectors}}
                                                               - GFF SORTED : {{input.gff_protein}}
                                                           Output:
                                                               - EFFECTOR PER CONTING : {{output.effector_per_contig}}
                                                           Others
                                                               - Threads : {{threads}}
                                                               - LOG error: {{log.error}}
                                                               - LOG output: {{log.output}}

                                                       """

    shell:
        """
        python {script_dir}count_effectors.py -g {input.gff_protein} -o {output.effector_per_contig} -fasta {input.fasta_effectors} 1>{log.output} 2>{log.error}
        sort -V {output.effector_per_contig} -o {output.effector_per_contig}
        """