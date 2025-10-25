FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="9cc9119bbcb06e107174c61dbc4aae63e1597c8f41f9799c2c015657e1ef52c1"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/bam2fastq.yml
#   prefix: /conda-envs/a7401219cb36035d7c6438fc301a8525
#   name: bam2fastq
#   channels:
#     - conda-forge
#   dependencies:
#     - bioconda::pbtk=3.5.0
RUN mkdir -p /conda-envs/a7401219cb36035d7c6438fc301a8525
COPY workflow/envs/bam2fastq.yml /conda-envs/a7401219cb36035d7c6438fc301a8525/environment.yaml

# Conda environment:
#   source: workflow/envs/busco.yml
#   prefix: /conda-envs/6a29b1058e04ecf1abe493d36682b23f
#   name: busco
#   channels:
#     - conda-forge
#   dependencies:
#     - bioconda::busco=6.0.0
RUN mkdir -p /conda-envs/6a29b1058e04ecf1abe493d36682b23f
COPY workflow/envs/busco.yml /conda-envs/6a29b1058e04ecf1abe493d36682b23f/environment.yaml

# Conda environment:
#   source: workflow/envs/cdskit.yml
#   prefix: /conda-envs/e508d70ff7615bb9d0dce44e5a9aa7b8
#   name: cdskit
#   channels:
#     - conda-forge
#   dependencies:
#     - bioconda::cdskit=0.14.4
#     - bioconda::seqkit=2.10.0
RUN mkdir -p /conda-envs/e508d70ff7615bb9d0dce44e5a9aa7b8
COPY workflow/envs/cdskit.yml /conda-envs/e508d70ff7615bb9d0dce44e5a9aa7b8/environment.yaml

# Conda environment:
#   source: workflow/envs/fastp.yml
#   prefix: /conda-envs/6f277ab7b19cb3d303e76a43f5b601ac
#   name: fastp
#   channels:
#     - conda-forge
#   dependencies:
#     - bioconda::fastp=1.0.1
RUN mkdir -p /conda-envs/6f277ab7b19cb3d303e76a43f5b601ac
COPY workflow/envs/fastp.yml /conda-envs/6f277ab7b19cb3d303e76a43f5b601ac/environment.yaml

# Conda environment:
#   source: workflow/envs/fastplong.yml
#   prefix: /conda-envs/3b2cd99f0beba959d7ba47922f256189
#   name: fastplong
#   channels:
#     - conda-forge
#   dependencies:
#     - bioconda::fastplong=0.3.0
RUN mkdir -p /conda-envs/3b2cd99f0beba959d7ba47922f256189
COPY workflow/envs/fastplong.yml /conda-envs/3b2cd99f0beba959d7ba47922f256189/environment.yaml

# Conda environment:
#   source: workflow/envs/fcs.yml
#   prefix: /conda-envs/519414a323ee447d07e07ff75a64ca53
#   name: fcs
#   channels:
#     - conda-forge
#   dependencies:
#     - python=3.13.7
RUN mkdir -p /conda-envs/519414a323ee447d07e07ff75a64ca53
COPY workflow/envs/fcs.yml /conda-envs/519414a323ee447d07e07ff75a64ca53/environment.yaml

# Conda environment:
#   source: workflow/envs/genomescope2.yml
#   prefix: /conda-envs/4862331702f5c049997163f0eecce518
#   name: genomescope2
#   channels:
#     - conda-forge
#   dependencies:
#     - bioconda::genomescope2=2.0.1
RUN mkdir -p /conda-envs/4862331702f5c049997163f0eecce518
COPY workflow/envs/genomescope2.yml /conda-envs/4862331702f5c049997163f0eecce518/environment.yaml

# Conda environment:
#   source: workflow/envs/hifiasm.yml
#   prefix: /conda-envs/0b12060b07954552849bc096ff5185bc
#   name: hifiasm
#   channels:
#     - conda-forge
#   dependencies:
#     - bioconda::hifiasm=0.25.0
RUN mkdir -p /conda-envs/0b12060b07954552849bc096ff5185bc
COPY workflow/envs/hifiasm.yml /conda-envs/0b12060b07954552849bc096ff5185bc/environment.yaml

# Conda environment:
#   source: workflow/envs/inspector.yml
#   prefix: /conda-envs/d02e04a35d8a81f2d038a0274829cfe5
#   name: inspector
#   channels:
#     - conda-forge
#   dependencies:
#     - bioconda::inspector=1.3.1
RUN mkdir -p /conda-envs/d02e04a35d8a81f2d038a0274829cfe5
COPY workflow/envs/inspector.yml /conda-envs/d02e04a35d8a81f2d038a0274829cfe5/environment.yaml

# Conda environment:
#   source: workflow/envs/jellyfish.yml
#   prefix: /conda-envs/e543a641f50fbe22ca4cab969e00243c
#   name: jellyfish
#   channels:
#     - conda-forge
#   dependencies:
#     - bioconda::kmer-jellyfish=2.3.1
RUN mkdir -p /conda-envs/e543a641f50fbe22ca4cab969e00243c
COPY workflow/envs/jellyfish.yml /conda-envs/e543a641f50fbe22ca4cab969e00243c/environment.yaml

# Conda environment:
#   source: workflow/envs/merqury.yml
#   prefix: /conda-envs/773e31c806ebde20c2b24480afdb55f5
#   name: merqury
#   channels:
#     - conda-forge
#   dependencies:
#     - bioconda::merqury=1.3
RUN mkdir -p /conda-envs/773e31c806ebde20c2b24480afdb55f5
COPY workflow/envs/merqury.yml /conda-envs/773e31c806ebde20c2b24480afdb55f5/environment.yaml

# Conda environment:
#   source: workflow/envs/oatk.yml
#   prefix: /conda-envs/b36f6158ce0b17af5a10b61984e0eba7
#   name: oatk
#   channels:
#     - conda-forge
#   dependencies:
#     - bioconda::oatk=1.0
RUN mkdir -p /conda-envs/b36f6158ce0b17af5a10b61984e0eba7
COPY workflow/envs/oatk.yml /conda-envs/b36f6158ce0b17af5a10b61984e0eba7/environment.yaml

# Conda environment:
#   source: workflow/envs/omamer.yml
#   prefix: /conda-envs/18f814805861521fcafda0f4d20eab60
#   name: omamer
#   channels:
#     - conda-forge
#   dependencies:
#     - bioconda::omamer=2.1.0
RUN mkdir -p /conda-envs/18f814805861521fcafda0f4d20eab60
COPY workflow/envs/omamer.yml /conda-envs/18f814805861521fcafda0f4d20eab60/environment.yaml

# Conda environment:
#   source: workflow/envs/omark.yml
#   prefix: /conda-envs/1767be13a2e7b9afacdae6bfe22ff2d8
#   name: omark
#   channels:
#     - conda-forge
#   dependencies:
#     - bioconda::omark=0.3.1
RUN mkdir -p /conda-envs/1767be13a2e7b9afacdae6bfe22ff2d8
COPY workflow/envs/omark.yml /conda-envs/1767be13a2e7b9afacdae6bfe22ff2d8/environment.yaml

# Conda environment:
#   source: workflow/envs/seqkit.yml
#   prefix: /conda-envs/2eace2598424741072f07c109da9f230
#   name: seqkit
#   channels:
#     - conda-forge
#   dependencies:
#     - bioconda::seqkit=2.10.0
RUN mkdir -p /conda-envs/2eace2598424741072f07c109da9f230
COPY workflow/envs/seqkit.yml /conda-envs/2eace2598424741072f07c109da9f230/environment.yaml

# Conda environment:
#   source: workflow/envs/smudgeplot.yml
#   prefix: /conda-envs/fba3c46ee9bfcd31e9dc9a0e556e010e
#   name: smudgeplot
#   channels:
#     - conda-forge
#   dependencies:
#     - bioconda::fastk=1.1.0
#     - bioconda::smudgeplot=0.4.0
RUN mkdir -p /conda-envs/fba3c46ee9bfcd31e9dc9a0e556e010e
COPY workflow/envs/smudgeplot.yml /conda-envs/fba3c46ee9bfcd31e9dc9a0e556e010e/environment.yaml

# Conda environment:
#   source: workflow/envs/tidk.yml
#   prefix: /conda-envs/2e2d95d2893efc5ad203d140f0a56688
#   name: tidk
#   channels:
#     - conda-forge
#   dependencies:
#     - bioconda::tidk=0.2.65
RUN mkdir -p /conda-envs/2e2d95d2893efc5ad203d140f0a56688
COPY workflow/envs/tidk.yml /conda-envs/2e2d95d2893efc5ad203d140f0a56688/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/a7401219cb36035d7c6438fc301a8525 --file /conda-envs/a7401219cb36035d7c6438fc301a8525/environment.yaml && \
    conda env create --prefix /conda-envs/6a29b1058e04ecf1abe493d36682b23f --file /conda-envs/6a29b1058e04ecf1abe493d36682b23f/environment.yaml && \
    conda env create --prefix /conda-envs/e508d70ff7615bb9d0dce44e5a9aa7b8 --file /conda-envs/e508d70ff7615bb9d0dce44e5a9aa7b8/environment.yaml && \
    conda env create --prefix /conda-envs/6f277ab7b19cb3d303e76a43f5b601ac --file /conda-envs/6f277ab7b19cb3d303e76a43f5b601ac/environment.yaml && \
    conda env create --prefix /conda-envs/3b2cd99f0beba959d7ba47922f256189 --file /conda-envs/3b2cd99f0beba959d7ba47922f256189/environment.yaml && \
    conda env create --prefix /conda-envs/519414a323ee447d07e07ff75a64ca53 --file /conda-envs/519414a323ee447d07e07ff75a64ca53/environment.yaml && \
    conda env create --prefix /conda-envs/4862331702f5c049997163f0eecce518 --file /conda-envs/4862331702f5c049997163f0eecce518/environment.yaml && \
    conda env create --prefix /conda-envs/0b12060b07954552849bc096ff5185bc --file /conda-envs/0b12060b07954552849bc096ff5185bc/environment.yaml && \
    conda env create --prefix /conda-envs/d02e04a35d8a81f2d038a0274829cfe5 --file /conda-envs/d02e04a35d8a81f2d038a0274829cfe5/environment.yaml && \
    conda env create --prefix /conda-envs/e543a641f50fbe22ca4cab969e00243c --file /conda-envs/e543a641f50fbe22ca4cab969e00243c/environment.yaml && \
    conda env create --prefix /conda-envs/773e31c806ebde20c2b24480afdb55f5 --file /conda-envs/773e31c806ebde20c2b24480afdb55f5/environment.yaml && \
    conda env create --prefix /conda-envs/b36f6158ce0b17af5a10b61984e0eba7 --file /conda-envs/b36f6158ce0b17af5a10b61984e0eba7/environment.yaml && \
    conda env create --prefix /conda-envs/18f814805861521fcafda0f4d20eab60 --file /conda-envs/18f814805861521fcafda0f4d20eab60/environment.yaml && \
    conda env create --prefix /conda-envs/1767be13a2e7b9afacdae6bfe22ff2d8 --file /conda-envs/1767be13a2e7b9afacdae6bfe22ff2d8/environment.yaml && \
    conda env create --prefix /conda-envs/2eace2598424741072f07c109da9f230 --file /conda-envs/2eace2598424741072f07c109da9f230/environment.yaml && \
    conda env create --prefix /conda-envs/fba3c46ee9bfcd31e9dc9a0e556e010e --file /conda-envs/fba3c46ee9bfcd31e9dc9a0e556e010e/environment.yaml && \
    conda env create --prefix /conda-envs/2e2d95d2893efc5ad203d140f0a56688 --file /conda-envs/2e2d95d2893efc5ad203d140f0a56688/environment.yaml && \
    conda clean --all -y
