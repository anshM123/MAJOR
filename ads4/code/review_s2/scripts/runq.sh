#!/bin/bash
# review rerun queue: runs the ANALYTIC_SLIVER2 pipeline from a private copy (ws/code/analytic2), 2 jobs at a time
R=/home/user/MAJOR/ads4/code/review_s2
export PYTHONDONTWRITEBYTECODE=1
cd $R/ws/code/analytic2
python3 rmass3.py 1 40 > $R/logs/rmass3_1_40.log 2>&1
python3 rmass3.py 40 898 > $R/logs/rmass3_40_898.log 2>&1
grep -v '^$' $R/scripts/jobs.txt | while read tag cmd; do echo "$tag|$cmd"; done | \
  xargs -P 2 -I{} bash -c 'l="{}"; tag=${l%%|*}; cmd=${l#*|}; cd '$R'/ws/code/analytic2; echo "start $tag $(date +%T)" >> '$R'/logs/queue.log; eval "$cmd" > '$R'/logs/$tag.log 2>&1; echo "done $tag rc=$? $(date +%T)" >> '$R'/logs/queue.log'
echo ALLDONE >> $R/logs/queue.log
