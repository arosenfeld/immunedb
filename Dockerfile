FROM ubuntu:14.04
RUN apt-get install -y software-properties-common
RUN apt-key adv --recv-keys --keyserver hkp://keyserver.ubuntu.com:80 0xcbcb082a1bb943db
RUN add-apt-repository 'deb [arch=amd64,i386] http://mirror.jmu.edu/pub/mariadb/repo/10.1/ubuntu trusty main'
RUN apt-get update
RUN sudo apt-get install -y python-numpy python-scipy python-setuptools supervisor mariadb-server
COPY setup.py /app/
COPY sldb/ /app/sldb
COPY lib/ /app/lib
COPY bin/ /app/bin
WORKDIR /app
RUN python setup.py install
COPY docker/supervisord.conf /etc/supervisor
COPY docker/my.cnf /etc/mysql
CMD ["/usr/bin/supervisord", "-c", "/etc/supervisor/supervisord.conf"]
