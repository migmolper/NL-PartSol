#!/bin/sh

clear

clang-format -i Von-Mises.c

gcc -DDEBUG_MODE=1 Von-Mises.c -o Von-Mises  -llapack -lm

./Von-Mises