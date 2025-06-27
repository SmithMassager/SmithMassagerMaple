#!/bin/bash

for dir in $(ls -d */)
do
  files=$(ls ${dir}/*.mpl)
  for file in ${files}
  do
    echo "running ${file}"
    maple -e0 -q ${file} | grep -i "error"
    if [ $? -eq 0 ]; then
      echo ${file}
    fi
  done
done
