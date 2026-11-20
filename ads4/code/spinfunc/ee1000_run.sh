#!/bin/bash
# run an analysis command while pausing ONE queue process (keeps <= 2 busy processes on this agent's two cores); usage: ./ee1000_run.sh cmd...
P=$(pgrep -f "ee1000_(rows|direct).py" | tail -1)
[ -n "$P" ] && kill -STOP $P
trap '[ -n "$P" ] && kill -CONT $P' EXIT
"$@"
