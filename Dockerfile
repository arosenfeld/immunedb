FROM ubuntu:18.04
# Get dependencies
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get install -y \
    python-setuptools \
    python3-venv \
    gcc \
    python3-dev \
    python3-setuptools \
    libdpkg-perl \
    mariadb-server \
    mariadb-client \
    make \
    wget \
    unzip \
    git \
    npm \
    r-base
WORKDIR /apps
# Get the frontend source, clearcut, and bowtie2
RUN git clone https://github.com/arosenfeld/immunedb-frontend
RUN wget http://bioinformatics.hungry.com/clearcut/clearcut-1.0.9.tar.gz && \
    tar xzf clearcut-1.0.9.tar.gz && mv clearcut-1.0.9 clearcut
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.3.4.1/bowtie2-2.3.4.1-linux-x86_64.zip && \
    unzip bowtie2-2.3.4.1-linux-x86_64.zip && \
    mv bowtie2-2.3.4.1-linux-x86_64 bowtie2
# Build the frontend, clearcut, bowtie2, and baseline
WORKDIR /apps/clearcut
RUN make
WORKDIR /apps/immunedb-frontend
RUN npm install
ENV BASENAME __BASENAME__
ENV API_ENDPOINT __ENDPOINT__
ENV NODE_ENV production
RUN npm run build
WORKDIR /apps/baseline
RUN wget http://immunedb.com/patched_baseline.tgz && \
    tar xzf patched_baseline.tgz
RUN Rscript -e 'install.packages(c("seqinr", "parallel"))'
# Copy ImmuneDB files and install
COPY requirements.txt setup.py README.md /apps/immunedb/
COPY lib/ /apps/immunedb/lib
COPY bin/ /apps/immunedb/bin
COPY immunedb/ /apps/immunedb/immunedb
WORKDIR /apps/immunedb
RUN python3 setup.py install
# Copy germlines and scripts
COPY docker/germlines/ /root/germlines
COPY docker/run.sh /root
COPY docker/proxy.py /apps/immunedb
COPY docker/mariadb/my.cnf /etc/mysql
COPY docker/setup_users.sql /tmp
COPY docker/example /example
# IgBLAST
ENV IGDATA /apps/igblast
WORKDIR /apps/igblast
COPY docker/run_igblast.sh /usr/sbin
COPY docker/make_db.sh /usr/sbin
RUN wget -q \
    ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/ncbi-igblast-\*-x64-linux.tar.gz && \
    tar xzf ncbi*.gz && \
    rm *.gz && \
    mv ncbi* src && \
    mv src/* . && \
    rm -r src
ENV PATH "${PATH}:/apps/bowtie2:/apps/clearcut:/apps/igblast/bin"
RUN make_db.sh database/human /root/germlines/igblast/human/*.fasta && \
    make_db.sh database/mouse /root/germlines/igblast/mouse/*.fasta && \
    mv database/human/*gapped* /root/germlines/igblast/human && \
    mv database/mouse/*gapped* /root/germlines/igblast/mouse
# Expose API and frontend ports
EXPOSE 5000 8080
# Setup MySQL volume
RUN mkdir -p /share
VOLUME /share
# Add the example data
WORKDIR /root
CMD bash -C 'run.sh';'bash'
