FROM ubuntu:20.04
MAINTAINER Amine Ghozlane "amine.ghozlane@pasteur.fr"
ARG CRAN_SOURCE
ENV TZ=Europe/Paris
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get update && apt-get install -y \
    wget \
    gdebi-core \
    libcurl4-openssl-dev\
    libcairo2-dev \
    libjpeg-dev \
    libtiff5-dev \
    libxt-dev \
    libxml2-dev \
    libxml2 \
    libreadline6-dev \
    git \
    libssl-dev \
    libssh2-1-dev \
    libnlopt-dev \
    python3-pip \
    python3-yaml \
    gcc \ 
    gfortran \ 
    g++ \
    make \
    openjdk-8-jdk \
    libmagick++-dev \
    tzdata \
    cmake
    
RUN git clone https://github.com/stevengj/nlopt.git && cd nlopt && mkdir build && cd build && cmake .. && make && make install

RUN pip3 install bioblend python-daemon==2.3.2

#Download and install shiny server
RUN wget --no-verbose ${CRAN_SOURCE}/src/base/R-3/R-3.6.2.tar.gz -P /opt/ && \
    tar -zxf /opt/R-3.6.2.tar.gz -C /opt && rm /opt/R-3.6.2.tar.gz && \
    cd /opt/R-3.6.2/ && ./configure --with-x=no && \
    make -j 4  && make install && cd / && rm -rf  /opt/R-3.6.2 && \  
    wget --no-verbose https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb


# Package R
COPY renv.lock /srv/shiny-server/
RUN R -e """install.packages('renv', repos='$CRAN_SOURCE');renv::restore(project='/srv/shiny-server/', prompt=F, repos='$CRAN_SOURCE')"""

# Configuration shiny
COPY docker_inst/shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY docker_inst/shiny-server.sh /usr/bin/shiny-server.sh

# Other software
COPY docker_inst/run_kronarshy.R /usr/bin/run_kronarshy.R
RUN git clone https://github.com/pierreLec/KronaRShy.git /srv/shiny-server/kronarshy && \
    git clone https://github.com/aghozlane/shaman_bioblend.git /usr/bin/shaman_bioblend && \
    chmod +x /usr/bin/shiny-server.sh &&  \
    mkdir -p /srv/shiny-server/www/masque/todo /srv/shiny-server/www/masque/doing /srv/shiny-server/www/masque/error /srv/shiny-server/www/masque/done && \
    chown -R shiny.shiny /srv/shiny-server/ && mkdir -p /var/log/shiny-server && chown shiny.shiny /var/log/shiny-server


# Copy of shaman
COPY . /srv/shiny-server/

EXPOSE 80

EXPOSE 5438

CMD ["/usr/bin/shiny-server.sh"] 
