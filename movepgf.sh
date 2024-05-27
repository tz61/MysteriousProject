#!/bin/zsh

for i in *.pgf
do
    cp "$i" "pre/$i"
    cp "$i" "report/$i"
    rm "$i"
done