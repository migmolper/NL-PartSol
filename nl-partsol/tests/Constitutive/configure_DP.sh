#!/bin/sh

clear

clang-format -i Drucker-Prager-Backward-Euler.c

gcc -DDEBUG_MODE=1 Drucker-Prager-Backward-Euler.c -o Drucker-Prager-Backward-Euler  -llapack -lm

./Drucker-Prager-Backward-Euler 