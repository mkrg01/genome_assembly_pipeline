# These rules use a host Conda environment so existing published pipeline images
# can run this optional stage without a container rebuild.

rule purge_dups_map:
    input:
        assembly = "results/fcs/assembly/{selected_assembly}/{assembly_name}.fa",
        reads = "results/hifi_reads/merged/{assembly_name}_hifi_reads_curated.fastq.gz"
    output:
        paf = "results/purge_dups/coverage/{selected_assembly}/{assembly_name}/reads.paf.gz",
        stat = "results/purge_dups/coverage/{selected_assembly}/{assembly_name}/PB.stat",
        cov = "results/purge_dups/coverage/{selected_assembly}/{assembly_name}/PB.base.cov"
    wildcard_constraints:
        selected_assembly = "|".join(purge_dups_targets) or "(?!)"
    log:
        "logs/purge_dups_map_{selected_assembly}_{assembly_name}.log"
    conda:
        "../envs/purge_dups.yml"
    container:
        None
    threads:
        min(16, max(1, workflow.cores))
    shell:
        """
        (
            mkdir -p $(dirname {output.stat:q})
            genome_size=$(python3 workflow/scripts/purge_dups_support.py size {input.assembly:q})
            minimap2 -t {threads} -I "$genome_size" -x map-hifi \
                {input.assembly:q} {input.reads:q} | gzip -c > {output.paf:q}
            pbcstat -O $(dirname {output.stat:q}) {output.paf:q}
            test -s {output.stat:q}
            test -s {output.cov:q}
        ) > {log:q} 2>&1
        """


rule purge_dups_cutoffs:
    input:
        "results/purge_dups/coverage/{selected_assembly}/{assembly_name}/PB.stat"
    output:
        cutoffs = "results/purge_dups/cutoffs/{selected_assembly}/{assembly_name}.txt",
        metadata = "results/purge_dups/cutoffs/{selected_assembly}/{assembly_name}.json"
    params:
        manual = "" if purge_dups_cutoffs == "auto" else "--manual " + " ".join(map(str, purge_dups_cutoffs))
    log:
        "logs/purge_dups_cutoffs_{selected_assembly}_{assembly_name}.log"
    conda:
        "../envs/purge_dups.yml"
    container:
        None
    shell:
        """
        python3 workflow/scripts/purge_dups_support.py cutoffs \
            --stat {input:q} --output {output.cutoffs:q} \
            --metadata {output.metadata:q} {params.manual} > {log:q} 2>&1
        """


rule purge_dups_self_alignment:
    input:
        "results/fcs/assembly/{selected_assembly}/{assembly_name}.fa"
    output:
        split = "results/purge_dups/self_alignment/{selected_assembly}/{assembly_name}.split.fa",
        paf = "results/purge_dups/self_alignment/{selected_assembly}/{assembly_name}.paf.gz"
    log:
        "logs/purge_dups_self_alignment_{selected_assembly}_{assembly_name}.log"
    conda:
        "../envs/purge_dups.yml"
    container:
        None
    threads:
        min(16, max(1, workflow.cores))
    shell:
        """
        (
            genome_size=$(python3 workflow/scripts/purge_dups_support.py size {input:q})
            split_fa {input:q} > {output.split:q}
            minimap2 -t {threads} -I "$genome_size" -x asm5 -DP \
                {output.split:q} {output.split:q} | gzip -c > {output.paf:q}
        ) > {log:q} 2>&1
        """


rule purge_dups_candidate:
    input:
        assembly = "results/fcs/assembly/{selected_assembly}/{assembly_name}.fa",
        paf = "results/purge_dups/self_alignment/{selected_assembly}/{assembly_name}.paf.gz",
        cov = "results/purge_dups/coverage/{selected_assembly}/{assembly_name}/PB.base.cov",
        cutoffs = "results/purge_dups/cutoffs/{selected_assembly}/{assembly_name}.txt"
    output:
        assembly = "results/purge_dups/assembly/{selected_assembly}/{assembly_name}.fa",
        removed = "results/purge_dups/removed/{selected_assembly}/{assembly_name}.fa",
        raw_bed = "results/purge_dups/bed/{selected_assembly}/{assembly_name}.raw.bed",
        bed = "results/purge_dups/bed/{selected_assembly}/{assembly_name}.haplotypic.bed"
    log:
        "logs/purge_dups_candidate_{selected_assembly}_{assembly_name}.log"
    conda:
        "../envs/purge_dups.yml"
    container:
        None
    shell:
        """
        (
            purge_dups -2 -T {input.cutoffs:q} -c {input.cov:q} {input.paf:q} > {output.raw_bed:q}
            python3 workflow/scripts/purge_dups_support.py filter-bed \
                --assembly {input.assembly:q} --bed {output.raw_bed:q} --output {output.bed:q}
            if test -s {output.bed:q}; then
                get_seqs -e -c -l 0 -m 0 -g 0 -p {output.assembly:q} {output.bed:q} {input.assembly:q}
                mv {output.assembly:q}.purged.fa {output.assembly:q}
                mv {output.assembly:q}.hap.fa {output.removed:q}
            else
                cp {input.assembly:q} {output.assembly:q}
                : > {output.removed:q}
            fi
            python3 workflow/scripts/purge_dups_support.py size {output.assembly:q}
        ) > {log:q} 2>&1
        """
