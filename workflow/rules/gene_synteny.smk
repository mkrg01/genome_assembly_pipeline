# Independent self-comparisons of representative genes for each selected assembly.
# Host Conda also supports existing pipeline containers without rebuilding them.

rule prepare_gene_synteny:
    input:
        genome = "results/repeatmasker/{selected_assembly}/{assembly_name}.fa.masked",
        proteins = "results/longest_cds/{selected_assembly}/{assembly_name}_aa.fa",
        gff3 = "results/longest_cds/{selected_assembly}/{assembly_name}.gff3",
        script = os.path.join(workflow.basedir, "scripts/gene_synteny.py")
    output:
        proteins = "results/gene_synteny/{selected_assembly}/{assembly_name}/prepared/proteins.fa",
        genes = "results/gene_synteny/{selected_assembly}/{assembly_name}/prepared/genes.tsv",
        gff = "results/gene_synteny/{selected_assembly}/{assembly_name}/prepared/data.gff",
        metadata = "results/gene_synteny/{selected_assembly}/{assembly_name}/prepared/inputs.json"
    params:
        outdir = lambda w, output: os.path.dirname(output.proteins)
    log:
        "logs/prepare_gene_synteny_{selected_assembly}_{assembly_name}.log"
    conda:
        "../envs/gene_synteny.yml"
    container:
        None
    shell:
        """
        python3 {input.script:q} prepare --genome {input.genome:q} \
            --proteins {input.proteins:q} --gff3 {input.gff3:q} \
            --output-dir {params.outdir:q} > {log:q} 2>&1
        """


rule diamond_gene_synteny:
    input:
        proteins = rules.prepare_gene_synteny.output.proteins,
        script = os.path.join(workflow.basedir, "scripts/gene_synteny.py")
    output:
        hits = "results/gene_synteny/{selected_assembly}/{assembly_name}/diamond/hits.tsv",
        metadata = "results/gene_synteny/{selected_assembly}/{assembly_name}/diamond/search.json"
    threads:
        min(16, max(1, workflow.cores))
    resources:
        mem_mb = 16000
    log:
        "logs/diamond_gene_synteny_{selected_assembly}_{assembly_name}.log"
    conda:
        "../envs/gene_synteny.yml"
    container:
        None
    shell:
        """
        python3 {input.script:q} search --proteins {input.proteins:q} \
            --output {output.hits:q} \
            --threads {threads} > {log:q} 2>&1
        """


rule mcscanx_gene_synteny:
    input:
        genes = rules.prepare_gene_synteny.output.genes,
        gff = rules.prepare_gene_synteny.output.gff,
        hits = rules.diamond_gene_synteny.output.hits,
        search = rules.diamond_gene_synteny.output.metadata,
        script = os.path.join(workflow.basedir, "scripts/gene_synteny.py")
    output:
        homology = "results/gene_synteny/{selected_assembly}/{assembly_name}/blocks/data.blast",
        gff = "results/gene_synteny/{selected_assembly}/{assembly_name}/blocks/data.gff",
        collinearity = "results/gene_synteny/{selected_assembly}/{assembly_name}/blocks/data.collinearity",
        anchors = "results/gene_synteny/{selected_assembly}/{assembly_name}/blocks/anchors.tsv",
        bedpe = "results/gene_synteny/{selected_assembly}/{assembly_name}/blocks/anchors.bedpe",
        blocks = "results/gene_synteny/{selected_assembly}/{assembly_name}/blocks/blocks.tsv",
        summary = "results/gene_synteny/{selected_assembly}/{assembly_name}/blocks/summary.json"
    params:
        prepared = lambda w, input: os.path.dirname(input.genes),
        outdir = lambda w, output: os.path.dirname(output.collinearity)
    resources:
        mem_mb = 16000
    log:
        "logs/mcscanx_gene_synteny_{selected_assembly}_{assembly_name}.log"
    conda:
        "../envs/gene_synteny.yml"
    container:
        None
    shell:
        """
        python3 {input.script:q} blocks --prepared {params.prepared:q} \
            --raw {input.hits:q} --output-dir {params.outdir:q} > {log:q} 2>&1
        """


rule plot_gene_synteny:
    input:
        genes = rules.prepare_gene_synteny.output.genes,
        metadata = rules.prepare_gene_synteny.output.metadata,
        anchors = rules.mcscanx_gene_synteny.output.anchors,
        summary = rules.mcscanx_gene_synteny.output.summary,
        script = os.path.join(workflow.basedir, "scripts/plot_gene_synteny.py"),
        support = os.path.join(workflow.basedir, "scripts/gene_synteny.py"),
        style = os.path.join(workflow.basedir, "scripts/plot_style.py"),
        labels = os.path.join(workflow.basedir, "scripts/plot_self_dotplot.py")
    output:
        gene_order = "results/gene_synteny/{selected_assembly}/{assembly_name}/plots/gene_order.pdf",
        mb = "results/gene_synteny/{selected_assembly}/{assembly_name}/plots/mb.pdf",
        contigs = "results/gene_synteny/{selected_assembly}/{assembly_name}/plots/contigs.tsv",
        summary = "results/gene_synteny/{selected_assembly}/{assembly_name}/plots/plot.json"
    params:
        prepared = lambda w, input: os.path.dirname(input.genes),
        blocks = lambda w, input: os.path.dirname(input.anchors),
        outdir = lambda w, output: os.path.dirname(output.mb),
        organism = config["organism_name"],
        min_length = config["min_long_contig_length"]
    resources:
        mem_mb = 4000
    log:
        "logs/plot_gene_synteny_{selected_assembly}_{assembly_name}.log"
    conda:
        "../envs/gene_synteny.yml"
    container:
        None
    shell:
        """
        python3 {input.script:q} --prepared {params.prepared:q} \
            --blocks {params.blocks:q} --output-dir {params.outdir:q} \
            --organism-name {params.organism:q} \
            --min-long-contig-length {params.min_length} > {log:q} 2>&1
        """
