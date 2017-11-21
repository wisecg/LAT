#!/bin/bash
echo "start:"
date
echo "shell: "$SHELL

if [[ $- == *i* ]]; then
  echo "interactive shell"
else
  # we're here.
  echo "noninteractive shell"
fi

echo "rootsys: "$ROOTSYS
echo "projdir: "$PROJDIR
echo $MANUALSW
echo ""