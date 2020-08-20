#!/bin/bash

CMD=$(command -v pip3)
if [ -z "$CMD" ]; then
    echo "pip3 not found! Exiting"
    exit 10
fi

while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
  -V | --version )
    cat requirements.txt
    pip3 -V
    exit
    ;;
  --prefix )
    shift; PREFIX=$1
    ;;
  --tests )
    shift; TESTS=1
    ;;
  --pipoptions )
    shift; PIPOPT=$1
    ;;
  -v | --verbose )
    echo 'Installing:'
    cat requirements.txt
    verbose='-v'
    ;;
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi

if [ -z "$PREFIX" ]; then
  sudo pip3 $verbose install -r requirements.txt ${PIPOPT}
else
  pip3 $verbose install -r requirements.txt --user ${PIPOPT}
fi 