#!/bin/bash
# src/make-4-figures.sh <flags-for-Python-script>
# flags are, of course, optional
set -e
for t in II-A II-C I-E I-F ; do src/make-figure.sh data/figure-data/$t-systems-for-figure.tab $t-figure "$1"; done
