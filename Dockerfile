FROM continuumio/miniconda3

# Install Git
RUN apt-get update \
    && apt-get install -y --no-install-recommends git apt-transport-https gnupg2 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

RUN git clone https://github.com/Gu-Lab-RBL-NCI/QuagmiR/ /opt/quagmir && \
    cd /opt/quagmir && \
    conda env create -f environment.yml

# RUN conda install -y python=3.4 bioconda::dropbox=5.2.1=py35_0 \
#         bioconda::filechunkio=1.6=py35_0 \
#         bioconda::ftputil=3.2=py35_0 \
#         bioconda::pysftp=0.2.8=py35_0 \
#         bioconda::python-dateutil=2.3=py35_0 \
#         bioconda::snakemake=3.7.1=py35_0 \
#         bioconda::urllib3=1.12=py35_0 \
#         biopython=1.67=np111py35_0 \
#         docutils=0.12=py35_2 \
#         ecdsa=0.13=py35_0 \
#         mkl=11.3.3=0 \
#         numpy=1.11.0=py35_1 \
#         openssl=1.0.1 \
#         pandas=0.18.1=np111py35_0 \
#         paramiko=1.16.0=py35_0 \
#         pip=8.1.2=py35_0 \
#         pycrypto=2.6.1=py35_4 \
#         python=3.5.0=0 \
#         pytz=2016.4=py35_0 \
#         pyyaml=3.11=py35_4 \
#         r::zlib=1.2.8 \
#         readline=6.2=2 \
#         requests=2.10.0=py35_0 \
#         setuptools=23.0.0=py35_0 \
#         six=1.10.0=py35_0 \
#         sqlite=3.13.0=0 \
#         tk=8.5.18=0 \
#         wheel=0.29.0=py35_0 \
#         xz=5.0.5=1 \
#         yaml=0.1.6=0 \
#         pip \\
#     && conda clean -y --all

# RUN pip install pystan==2.17.0.0 \
#   dropbox==5.2.1 \
#   filechunkio==1.6 \
#   ftputil==3.2 \
#   pysftp==0.2.8 \
#   python-dateutil==2.3 \
#   snakemake==3.7.1 \
#   urllib3==1.12 \
#   git+https://github.com/infoscout/weighted-levenshtein.git