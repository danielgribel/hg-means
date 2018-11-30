#!/bin/bash

# Population size (\Pi_min): {5, 10}
pi_min=10

# Max number of iterations (N2): {500, 5000}
n2=5000

# Boolean for clustering evaluation: {0, 1}. IMPORTANT: If 1 is set, ground-truth label must be provided in /labels folder
evaluate=0

# Benchmark instance
./exec.sh "data/fisher.txt $pi_min $n2 $evaluate 2 3 4 5 6 7 8 9 10"
