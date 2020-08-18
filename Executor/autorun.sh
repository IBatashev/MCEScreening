#!/bin/bash

# A simple script to keep submitting screening jobs as long as queue is not too full
# Set time until this script is killed and maximum number of the jobs to be allowed in the queue
# run the script by: " nohup ./autorun.sh [timeout] [job_limit] & " to detach it from ssh session

timeout=${1:-4h}
job_limit=${2:-30}

{
    sleep "$timeout"
    echo "Autorun script killed  at $(date)."
    echo "Autorun script killed  at $(date)." >> auto_report
    echo "---" >> auto_report

    kill $$
} &
echo "Autorun script started at $(date). Timeout is $timeout. Job limit is $job_limit"
echo "Autorun script started at $(date). Timeout is $timeout. Job limit is $job_limit" >> auto_report

while true; do
currenttime=$(date)
jobs_in_queue=$(squeue | grep "ibatashe" | wc -l)
if (( $jobs_in_queue < $job_limit )) ; then
  ./batch_submit_30
  echo "30 jobs submitted      at $currenttime. Number of active jobs was $jobs_in_queue"
  echo "30 jobs submitted      at $currenttime. Number of active jobs was $jobs_in_queue" >> auto_report
else
  echo "Queue full             at $currenttime. Number of active jobs was $jobs_in_queue ...waiting"
  echo "Queue full             at $currenttime. Number of active jobs was $jobs_in_queue ...waiting" >> auto_report
fi

sleep 20m

done
