# Cellula image

################################################################

## Start from r-base:4.2.0
FROM rocker/r-base:4.2.0

## Change as u want
LABEL mantainer = "andrea.cossa@ieo.it"

## Compilers, utilities, git
RUN apt-get --quiet update --yes \
  && apt-get --quiet upgrade --yes \
  && apt-get --quiet install --yes \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    build-essential \
    libz-dev \
    wget \
    procps \
    git \
    cmake \
  && rm -rf /var/lib/apt/lists/*

## Mamba
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh \
    && bash Mambaforge-Linux-x86_64.sh -b -p /mambaforge
ENV PATH="/mambaforge/bin/:$PATH"

## Clone repo and install code
RUN git clone https://github.com/andrecossa5/Cellula.git 
RUN echo ciao
RUN cd Cellula \
    && git pull \
    && mamba env update -n base -f envs/new_Linux.yml \
    && mamba develop .

## Default command to initialize the container
CMD [ "bash" ]

#################################################################
