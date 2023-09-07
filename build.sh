#!/bin/sh

set -xe

CFLAGS="-Wall -Wextra -g"
LIBS="-lm"

mkdir -p ./build

gcc $CFLAGS -o ./build/out.a stats.c $LIBS
