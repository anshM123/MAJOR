#!/bin/bash
# fw_after.sh PATTERN SCRIPT: wait until no process matches PATTERN (pgrep -f), then run SCRIPT (resumable fw queues).
cd "$(dirname "$0")"
while pgrep -f "^python3 $1" > /dev/null; do sleep 20; done
exec "./$2"
