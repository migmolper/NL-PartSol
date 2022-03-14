#!/bin/sh

clear

clang-format -i Frictional-Monolithic.c

gcc -DDEBUG_MODE=1 Frictional-Monolithic.c -o Frictional-Monolithic  -llapack -lm

./Frictional-Monolithic