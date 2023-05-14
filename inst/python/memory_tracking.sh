#!/bin/bash

SCRIPT=$1
shift
ARGS="$@"

python "$SCRIPT" $ARGS &
PID=$!

while [ -e /proc/$PID ]; do
    RSS=$(ps -p $PID -o rss=)
    RSS_MB=$((RSS / 1024))
    echo "$(date +%s) $RSS_MB MB"
    echo "$(date +%s) $RSS_MB MB" >> memory_usage.log
    sleep 5
done
