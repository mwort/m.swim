#!/usr/bin/env bash

for f in ../m.swim.*/m.swim.*.html; do
  ${f%.*}.py --html-description > tmp_description
  l=$(wc -l tmp_description)
  echo $l $f
  lt=$((${l% *} - 2))
  head -n $lt tmp_description > tmp_header
  tail -n 2 tmp_description > tmp_footer
  echo "<p>Last updated (docs/code): $(git log -1 --format=%cd $f) / $(git log -1 --format=%cd ${f%.*}.py) </p>" > tmp_date
  cat tmp_header $f tmp_date tmp_footer > ${f##*/}
  rm tmp_*
done
