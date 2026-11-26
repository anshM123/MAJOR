#!/bin/bash
cd /tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/direct_ext
s=$(date +%s.%N); python3 direct_sign_tight_ext.py $1 $2 $3 $4 > runs/tight_$1_$2_$3.log 2>&1; rc=$?
echo "rc=$rc wall=$(echo "$(date +%s.%N)-$s"|bc)" >> runs/tight_$1_$2_$3.log
