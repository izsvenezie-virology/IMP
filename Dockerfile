FROM mambaorg/micromamba:latest AS micromamba


FROM nextflow/nextflow:25.04.8

COPY --from=micromamba /usr/bin/micromamba /usr/bin/micromamba

RUN yum install -y git && \
    micromamba shell init -s bash -r /env
    
RUN curl -fsSL https://raw.githubusercontent.com/izsvenezie-virology/IMP/refs/heads/main/environment.yml -o /environment.yml && \
    micromamba create -p /env -f /environment.yml -y && \
    rm /environment.yml

RUN git clone https://github.com/izsvenezie-virology/IMP.git /imp

VOLUME /pipeline_input
VOLUME /pipeline_output
VOLUME /pipeline_resources
VOLUME /pipeline_cache

WORKDIR /pipeline_input

ENV NXF_LOG_FILE="/pipeline_cache/.nextflow.log"
ENV NXF_CACHE_DIR="/pipeline_cache/.nextflow"
ENV NXF_WORK="/pipeline_cache/work"

CMD ["nextflow", "/imp/src/main.nf", "-resume", "-profile", "docker", "-c", "/pipeline_input/imp.config"]
