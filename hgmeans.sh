#!/bin/bash

# Population size (\Pi_min): {5, 10}
pi_min=10

# Max number of iterations (N2): {500, 5000}
n2=5000

# Boolean for clustering evaluation: {0, 1}
# IMPORTANT: If 1 is set, ground-truth label must be provided in /labels folder
evaluate=0

# Benchmark instances A1
# ./exec.sh "data/fisher.txt $pi_min $n2 $evaluate 2 3 4 5 6 7 8 9 10"
./exec.sh "data/ionosphere.txt $pi_min $n2 $evaluate 2 5 10 15 20 25 30 40 50"