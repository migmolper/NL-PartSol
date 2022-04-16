#!/bin/sh

clear

clang-format -i Matsuoka_Nakai.c

gcc -DDEBUG_MODE=0 Matsuoka_Nakai.c -o Matsuoka_Nakai  -llapack -lm

./Matsuoka_Nakai