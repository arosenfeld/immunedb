FROM ubuntu:14.04
RUN apt-get update
RUN apt-get install -y software-properties-common
RUN apt-key adv --recv-keys --keyserver hkp://keyserver.ubuntu.com:80 0xcbcb082a1bb943db
RUN apt-get update && apt-get install -y python-numpy python-setuptools wget python-dev
COPY setup.py /app/
COPY immunedb/ /app/immunedb
COPY lib/ /app/lib
COPY bin/ /app/bin
WORKDIR /app
RUN python setup.py install
RUN wget https://raw.githubusercontent.com/vishnubob/wait-for-it/master/wait-for-it.sh
RUN chmod +x wait-for-it.sh
RUN mkdir /root/configs /root/data
COPY docker/configs/immunedb.json /root/configs/immunedb.json
COPY docker/germlines/ /root/germlines
WORKDIR /root
CMD /app/./wait-for-it.sh -t 0 mariadb:3306 -- \
    immunedb_admin create immunedb /root/configs --db-host mariadb --admin-pass insecure_password && \
    immunedb_rest /root/configs/immunedb.json
