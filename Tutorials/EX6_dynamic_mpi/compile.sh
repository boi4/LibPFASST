#!/usr/bin/env bash
set -euo pipefail
#set -x
#

dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
cd "$dir"/../../
echo "Making LibPFASST"
make DEBUG=TRUE VERBOSE=TRUE
cd "$dir"
#MPIDYNRES=../../../libmpidynres
#cp $MPIDYNRES/*.mod .
#cp $MPIDYNRES/build/tmp/*.f90.o
echo "Making Adv. Program"
make DEBUG=TRUE VERBOSE=TRUE
