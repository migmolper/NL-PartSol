#!/bin/sh

clear

clang-format -i Matsuoka_Nakai.c

gcc -DDEBUG_MODE=1 Matsuoka_Nakai.c -o Matsuoka_Nakai  -llapack -lm

./Matsuoka_Nakai