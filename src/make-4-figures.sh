#!/bin/bash
for t in II-A II-C I-E I-F ; do src/make-figure.sh data/figure-data/$t-systems-for-figure.tab $t-figure; done
