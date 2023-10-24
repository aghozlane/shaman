FROM ubuntu:18.04
MAINTAINER Amine Ghozlane "amine.ghozlane@pasteur.fr"
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
    
RUN git clone https://github.com/stevengj/nlopt.git && cd nlopt && mkdir build && cd build && cmake .. && make && sudo make install

RUN pip3 install bioblend python-daemon==2.3.2

RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

#Download and install shiny server
RUN wget --no-verbose https://cran.r-project.org/src/base/R-3/R-3.6.1.tar.gz -P /opt/ && \
    tar -zxf /opt/R-3.6.1.tar.gz -C /opt && rm /opt/R-3.6.1.tar.gz && \
    cd /opt/R-3.6.1/ && ./configure --with-x=no && \
    make -j 4  && make install && cd / && rm -rf  /opt/R-3.6.1 && \  
    wget --no-verbose https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb



COPY docker_inst/shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY docker_inst/.Rprofile  /srv/shiny-server/
COPY . /srv/shiny-server/
COPY docker_inst/shiny-server.sh /usr/bin/shiny-server.sh
COPY docker_inst/run_kronarshy.R /usr/bin/run_kronarshy.R

RUN git clone https://github.com/pierreLec/KronaRShy.git /srv/shiny-server/kronarshy && \
    git clone https://github.com/aghozlane/shaman_bioblend.git /usr/bin/shaman_bioblend && \
    chown -R shiny.shiny  /srv/shiny-server/ && \
    cp /srv/shiny-server/.Rprofile /srv/shiny-server/kronarshy/.Rprofile && \
    chmod +x /usr/bin/shiny-server.sh
WORKDIR /srv/shiny-server/
RUN R -e """renv::restore(prompt=F)"""

EXPOSE 80

EXPOSE 5438

CMD ["/usr/bin/shiny-server.sh"] 
