FROM mambaorg/micromamba:latest AS micromamba


FROM nextflow/nextflow:25.10.2

COPY --from=micromamba /usr/bin/micromamba /usr/bin/micromamba

RUN yum install -y git && \
    micromamba shell init -s bash -r /env
    
RUN curl -fsSL https://raw.githubusercontent.com/izsvenezie-virology/IMP/refs/heads/main/environment.yml -o /environment.yml && \
    micromamba create -p /env -f /environment.yml -y && \
    rm /environment.yml

RUN git clone https://github.com/izsvenezie-virology/IMP.git /imp

VOLUME /pipeline_input
VOLUME /pipeline_output

WORKDIR /pipeline_input

CMD ["nextflow", "/imp/src/main.nf", "-resume", "-profile", "docker", "-c", "/pipeline_input/imp.config"]
