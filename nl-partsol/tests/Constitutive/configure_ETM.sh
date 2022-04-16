#!/bin/sh

clear

clang-format -i Elastoplastic-Tangent-Matrix.c

gcc -DDEBUG_MODE=1 Elastoplastic-Tangent-Matrix.c -o Elastoplastic-Tangent-Matrix  -llapack -lm

./Elastoplastic-Tangent-Matrix 