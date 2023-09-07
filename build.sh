#!/bin/sh

set -xe

CFLAGS="-Wall -Wextra"
LIBS="-lm"

mkdir -p ./build

gcc $CFLAGS -o ./build/out.a stats.c $LIBS
