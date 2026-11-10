# Shared assembly preparation and databases for every entrypoint.

rule extract_long_contigs:
    input:
        assembly = "results/{assembly}/assembly/{selected_assembly}/{assembly_name}.fa",
        config_stamp = "config/config.yml"
    output:
        "results/{assembly}/assembly_long_contigs/{selected_assembly}/{assembly_name}.fa"
    log:
        out = "logs/extract_long_contigs_{assembly}_{selected_assembly}_{assembly_name}.out",
        err = "logs/extract_long_contigs_{assembly}_{selected_assembly}_{assembly_name}.err"
    conda:
        "../envs/seqkit.yml"
    params:
        min_long_contig_length = config["min_long_contig_length"]
    shell:
        """
        (
            seqkit seq --min-len {params.min_long_contig_length} {input.assembly:q} \
            | seqkit sort --by-length --reverse \
            > {output:q}
        ) > {log.out:q} 2> {log.err:q}
        """


rule download_busco_database:
    output:
        directory("results/downloads/busco_downloads")
    log:
        out = "logs/download_busco_database.out",
        err = "logs/download_busco_database.err"
    conda:
        "../envs/busco.yml"
    params:
        lineage_dataset = config["busco_lineage_dataset"]
    shell:
        """
        (
            busco --download {params.lineage_dataset}
            mkdir -p {output:q}
            mv busco_downloads/* {output:q}
            rm -rf busco_downloads
        ) > {log.out:q} 2> {log.err:q}
        """
