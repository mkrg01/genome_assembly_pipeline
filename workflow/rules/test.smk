rule test_conda:
    # input:
    #     "raw_data/test_conda.txt"
    output:
        "results/test_conda.txt"
    log:
        "logs/test_conda.log"
    conda:
        "../envs/test.yml"
    shell:
        "Rscript workflow/scripts/test.R > {output} 2> {log}"

rule test_container:
    # input:
    #     "raw_data/test_container.txt"
    params:
        cowsay = {config['cowsay']}
    output:
        "results/test_container.txt"
    log:
        "logs/test_container.log"
    container:
        "docker://aurelia01/cowsay:v1.1.1"
    shell:
        "fortune | cowsay -f {params.cowsay} > {output} 2> {log}"
