#!/bin/bash

# Population size (\Pi_min): {5, 10}
pi_min=10

# Max number of iterations (N2): {500, 5000}
n2=5000

# Boolean for clustering evaluation: {0, 1}
# IMPORTANT: If 1 is set, ground-truth labels must be provided in /labels folder with the same name of the dataset
evaluate=1

# Benchmark instances A1
# ./exec.sh "data/bavaria1.txt $pi_min $n2 $evaluate 2 3 4 5 6 7 8 9 10"
# ./exec.sh "data/bavaria2.txt $pi_min $n2 $evaluate 2 3 4 5 6 7 8 9 10"
# ./exec.sh "data/fisher.txt $pi_min $n2 $evaluate 2 3 4 5 6 7 8 9 10"
# ./exec.sh "data/german.txt $pi_min $n2 $evaluate 2 3 4 5 6 7 8 9 10"

# Benchmark instances A2
# ./exec.sh "data/breast.txt $pi_min $n2 $evaluate 2 5 10 15 20 25 30 40 50"
# ./exec.sh "data/congress.txt $pi_min $n2 $evaluate 2 5 10 15 20 25 30 40 50"
# ./exec.sh "data/heart.txt $pi_min $n2 $evaluate 2 5 10 15 20 25 30 40 50"
# ./exec.sh "data/ionosphere.txt $pi_min $n2 $evaluate 2 5 10 15 20 25 30 40 50"
# ./exec.sh "data/liver.txt $pi_min $n2 $evaluate 2 5 10 15 20 25 30 40 50"
# ./exec.sh "data/pima.txt $pi_min $n2 $evaluate 2 5 10 15 20 25 30 40 50"

# Mixture of Gaussians instances
./exec.sh "data/Gau-10-20-1.txt $pi_min $n2 $evaluate 10"
./exec.sh "data/Gau-10-50-1.txt $pi_min $n2 $evaluate 10"
./exec.sh "data/Gau-10-100-1.txt $pi_min $n2 $evaluate 10"
./exec.sh "data/Gau-10-200-1.txt $pi_min $n2 $evaluate 10"
./exec.sh "data/Gau-10-500-1.txt $pi_min $n2 $evaluate 10"
./exec.sh "data/Gau-20-20-1.txt $pi_min $n2 $evaluate 20"
./exec.sh "data/Gau-20-50-1.txt $pi_min $n2 $evaluate 20"
./exec.sh "data/Gau-20-100-1.txt $pi_min $n2 $evaluate 20"
./exec.sh "data/Gau-20-200-1.txt $pi_min $n2 $evaluate 20"
./exec.sh "data/Gau-20-500-1.txt $pi_min $n2 $evaluate 20"