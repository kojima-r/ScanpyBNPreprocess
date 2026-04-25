#!/bin/sh
# Backward-compatible wrapper: bootstrap BN with N=100 over 03data_bbknn/all2.txt.
sh "$(dirname "$0")/05run.sh" 03data_bbknn/all2.txt 100
