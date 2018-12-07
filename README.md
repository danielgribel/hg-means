# HG-means

Source code of HG-means clustering, an efficient hybrid genetic algorithm proposed for minimum sum-of-squares clustering (MSSC). This population-based metaheuristic uses K-means as a local search in combination with crossover, mutation, and diversification operators. As HG-means algorithm uses K-means, we included the fundamental source files of the fast K-means implementation of Greg Hamerly (to whom we are grateful to make available the source code) in this repository.

Original files and complete source code of Greg Hamerly K-means can be found at: https://github.com/ghamerly/fast-kmeans

# Related Article

"HG-means: A scalable hybrid genetic algorithm for minimum sum-of-squares clustering". Technical Report PUC-Rio -- arXiv 1804.09813. Gribel, D.; Vidal, T.

https://arxiv.org/abs/1804.09813

# Run

Firstly, make sure that the datasets are placed in `/data` folder.

The program parameters as well as the input data should be defined in `hgmeans.sh` script,
which executes the algorithm.

Important: If the parameter `evaluate` -- for calculating external clustering measures -- is enabled in `hgmeans.sh` script, make sure that the ground-truth file is provided in `/labels` folder with the same name of the dataset.

To run the algorithm, try the following sequence of commands:

`> make`

`> ./hgmeans.sh`

After the execution of the algorithm, output files will be saved in `/out` folder.

# Data format

DATASET files: In the first line of a dataset file, the number of data points (n) and the dimensionality of the data (d) is set, separated by a single space. The remaining lines correspond to the coordinates of data points. Each line contains the values of the d features of a sample, where x_ij correspond to the j-th feature of the i-th sample of the data. Each feature value is separated by a single space, as depicted in the scheme below:

---------------------------
| n d                     |
---------------------------
| x_11 x_12 x_13 ... x_1d |
---------------------------
| x_21 x_22 x_23 ... x_2d |
---------------------------
| ...  ...  ...  ... ...  |
---------------------------
| x_n1 x_n2 x_n3 ... x_nd |
---------------------------

LABELS files: The content of a labels file exhibits the cluster of each sample of the dataset according to ground-truth, where y_i correspond to the label of the i-th sample:

-----
 y_1
-----
 y_2
-----
 ...
-----
 y_n
-----