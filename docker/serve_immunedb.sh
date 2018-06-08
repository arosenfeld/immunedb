if [ ! "$#" -eq 1 ]
then
    echo "Database name must be provided"
    exit
fi

function teardown() {
    trap '' INT TERM
    echo 'Stopping REST API'
    kill -9 $RPID
    echo 'Stopping frontend'
    kill -9 $FPID
}

trap teardown SIGINT SIGTERM EXIT

echo "Running for database $1"
IP=`awk 'END{print $1}' /etc/hosts`
immunedb_rest $1 -p 5000 &
RPID=$!
cd /apps/immunedb-frontend
NODE_END=production API_ENDPOINT=http://$IP:5000 npm run serve &
FPID=$!

echo "*** Website is now served at http://$IP:8080 ***"
wait $FPID
