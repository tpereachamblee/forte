############################# BASE IMAGE ############################

FROM python:3.12.2

############################## METADATA #############################

LABEL author="Juan Blanco Heredia <blancoj@mskcc.org>"

LABEL version="7.4.0"

LABEL software="FusViz"

############################ INSTALLATION ###########################

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    time \
    procps \
    samtools \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    ghostscript \
    liblzma-dev \
    poppler-utils \
    build-essential \
    fonts-dejavu-core \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

ENV PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg

RUN pip install --upgrade pip setuptools wheel && \
    pip install --no-cache-dir \
        numpy==1.26.4 \
        pandas==2.2.3 \
        pysam==0.22.1 \
        argparse==1.1  \
        pyranges==0.1.4 \
        pycirclize==1.7.1 \
        pdf2image==1.17.0 \
        matplotlib==3.10.0

ENV MPLCONFIGDIR=/tmp/matplotlib

ENV MPLBACKEND=Agg

RUN mkdir -p /tmp/matplotlib && chmod 777 /tmp/matplotlib

WORKDIR /usr/local/bin/

COPY FusViz_v7.4.0.py /usr/local/bin/

RUN chmod 775 /usr/local/bin/FusViz_v7.4.0.py

RUN ln -s /usr/local/bin/FusViz_v7.4.0.py /usr/local/bin/FusViz

WORKDIR /data

CMD ["/bin/bash"]
