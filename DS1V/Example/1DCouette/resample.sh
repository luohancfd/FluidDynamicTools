#!/bin/bash
set +e

TSLEEP="5s"

FILE=$1
PIDKILL=$2
NLINE=$3
TORUN=$4
TOLOG=$5

echo $4
n=0
while [ $n -lt $NLINE ];do
  sleep $TSLEEP
  if [ -f $FILE ];then
    n=$(wc -l < $FILE)
  else
    n=0
  fi
done
echo "$(date) kill $2" >> ${TOLOG}
kill -9 $PIDKILL
nohup bash -c "${TORUN}" > /dev/null 2>&1 &
PID=$!
echo "$(date) run ${PID}" >> ${TOLOG}

