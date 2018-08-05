# Define base image
FROM ubuntu:14.04

WORKDIR /opt

#apt-get install dependencies
RUN apt-get update && apt-get install -y \
wget \
git \
zlib1g-dev \
software-properties-common \
python-software-properties \
python-dev \
python3-dev \
libmysqlclient-dev \
build-essential
RUN add-apt-repository ppa:fkrull/deadsnakes -y
RUN apt-get update && apt-get install python3.5 -y

#install pip
RUN wget https://bootstrap.pypa.io/get-pip.py
RUN python3.5 get-pip.py
RUN rm -rf get-pip.py

#install miniconda
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh
RUN wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O ~/miniconda.sh
RUN /bin/bash ~/miniconda.sh -b -p /opt/conda
RUN rm ~/miniconda.sh
ENV PATH /opt/conda/bin:$PATH

#install quagmir
RUN git clone https://github.com/duxan/quagmir
RUN conda env create -f quagmir/environment.yml

#restart following steps
RUN echo 13

WORKDIR /opt/quagmir
RUN git fetch
RUN git checkout adding_gff_13
WORKDIR /opt/

RUN rm -rf quagmir/motif-consensus.fa
RUN rm -rf quagmir/config.yaml
RUN rm -rf quagmir/data

COPY Dockerfile /opt/
MAINTAINER Nikola Tesic, Seven Bridges, <nikola.tesic@sbgenomics.com>