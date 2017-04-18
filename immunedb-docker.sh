#!/bin/bash
export API_ENDPOINT=${API_ENDPOINT:-http://localhost:5000}
export BASENAME=${BASENAME:-}
export DB_VOLUME=${DB_VOLUME:-~/immunedb/db}
export DATA_VOLUME=${DATA_VOLUME:-~/immunedb/data}
export API_PORT=${API_PORT:-5000}
export SERVE_PORT=${SERVE_PORT:-8080}
export PRODUCTION=${PRODUCTION:-0}

case $1 in
	'start' )
        docker-compose up;;
    'stop' )
        docker-compose down;;
    'shell' )
        ID=`docker ps | awk '{print $1" "$2}' | grep immunedb | grep -v frontend | awk '{print $1}'`
        instances=`echo $ID | wc -w`
        if [ $instances -eq 0 ]
        then
            echo '    There does not appear to be a running ImmuneDB instance'
        else
            docker exec -i -t $ID bash
        fi;;
    *)
        echo >&2 "usage: ./immunedb.sh COMMAND
Valid Commands:
    start      Starts an ImmuneDB docker instance
    stop       Stops an ImmuneDB docker instance
    shell      Launches shell into an ImmuneDB instance";
        exit 1;;
esac
