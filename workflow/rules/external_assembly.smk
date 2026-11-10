# Preserve external FASTA identifiers so supplied annotations stay aligned.

rule stage_external_assembly:
    input:
        assembly = lambda wildcards: external_assembly,
        config_stamp = "config/config.yml"
    output:
        "results/external/assembly/{selected_assembly}/{assembly_name}.fa"
    log:
        out = "logs/stage_external_assembly_{selected_assembly}_{assembly_name}.out",
        err = "logs/stage_external_assembly_{selected_assembly}_{assembly_name}.err"
    conda:
        "../envs/seqkit.yml"
    shell:
        """
        (
            mkdir -p $(dirname {output:q})
            gzip -cdf {input.assembly:q} > {output:q}
        ) > {log.out:q} 2> {log.err:q}
        """
