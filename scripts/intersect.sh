#!/bin/sh

intersect() {
    for file in "$@"; do
        sort -u "$file"
    done | sort | uniq -cd | grep "^[^0-9]*$# " | sed -e "s/      [0-9]//" | sed -e "s/ //"
}


intersect $1 $2 $3 $4
