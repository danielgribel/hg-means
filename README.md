# HG-means
Source code of HG-means clustering, an efficient hybrid genetic algorithm proposed for minimum sum-of-squares clustering (MSSC). This population-based metaheuristic uses K-means as a local search in combination with crossover, mutation, and diversification operators. As the proposed algorithm uses K-means, we included the fundamental source files of the fast K-means implementation of Greg Hamerly in this repository. Original files and complete source code of this K-means version can be found at: https://github.com/ghamerly/fast-kmeans

# Related Article
"HG-means: A scalable hybrid genetic algorithm for minimum sum-of-squares clustering"

# Run
Firstly, make sure that the datasets are placed in `/data` folder.
The program parameters as well as the input data should be defined in `hgmeans.sh` script,
which executes the algorithm.
Important: If the parameter `evaluate` -- for calculating external clustering measures -- is enabled in `hgmeans.sh` script, make sure that the ground-truth file is provided in `/labels` folder with the same name of the dataset.
To run the algorithm, try the following sequence of commands:
`> make`

`> ./hgmeans.sh`
After the execution of the algorithm, output files will be saved in `/out` folder.