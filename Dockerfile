FROM ubuntu:14.04
RUN apt-get install -y software-properties-common
RUN apt-key adv --recv-keys --keyserver hkp://keyserver.ubuntu.com:80 0xcbcb082a1bb943db
RUN apt-get update && apt-get install -y python-numpy python-scipy python-setuptools wget
COPY setup.py /app/
COPY sldb/ /app/sldb
COPY lib/ /app/lib
COPY bin/ /app/bin
WORKDIR /app
RUN python setup.py install
RUN wget https://raw.githubusercontent.com/vishnubob/wait-for-it/master/wait-for-it.sh
RUN chmod +x wait-for-it.sh
RUN mkdir /configs /data
COPY docker/configs/sldb.json /configs/sldb.json
WORKDIR /
CMD /app/./wait-for-it.sh -t 0 mariadb:3306 -- \
    sldb_admin create root sldb /configs --db-host mariadb --admin-pass insecure_password && \
    sldb_rest /configs/sldb.json
