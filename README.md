# HG-means

Source code of HG-means clustering, an efficient hybrid genetic algorithm proposed for the minimum sum-of-squares clustering (MSSC). This population-based metaheuristic uses K-means as a local search in combination with crossover, mutation, and diversification operators.

As HG-means algorithm uses K-means, we included the fundamental source files of the fast K-means implementation of Greg Hamerly (to whom we are grateful for making the source code available) in this repository, under the folder `/hamerly`. Original files and complete source code of Greg Hamerly K-means can be found at: https://github.com/ghamerly/fast-kmeans.

For the exact crossover, HG-means uses the implementation of Dlib (https://github.com/davisking/dlib) for solving an assignment problem. Dlib files are included in `/dlib-master` folder.

HG-means clustering is available as a C++ code, as well as a Python package.

## Related Article

*HG-means: A scalable hybrid genetic algorithm for minimum sum-of-squares clustering*. D. Gribel and T. Vidal, 2019. Pattern Recognition, https://doi.org/10.1016/j.patcog.2018.12.022

## Installation and Run

### C++

To run the algorithm in C++, go to `/hgmeans` folder and try the following sequence of commands:

`> make`

`> ./hgmeans 'dataset_path' pi_min n2 it [nb_clusters] 'w'`

### Example

`> make`

`> ./hgmeans 'data/iris.txt' 10 5000 1 2 5 10 'w'`

This script executes HG-means clustering for "iris" dataset, with 10 solutions in population, a maximum of 5000 iterations, 1 iteration (algorithm repetition), and 2, 5 and 10 clusters.

**Important:** You can provide a ground-truth file with the labels of clusters. In this case, make sure that a file with the same name of the dataset and '.label' extension is placed in the same folder of the dataset. If this file is provided, HG-means clustering will compute clustering performance metrics. See section **Data format** to check the expected data format for datasets and labels files.

### Parameters of the algorithm

`dataset_path`: The path of dataset.

`pi_min` (default = 10): Population size. Determines the size of the population in the genetic algorithm.

`n2` (default = 5000): Maximum number of iterations. Determines the total number of iterations the algorithm will take.

`it` (default = 1): The number of independent repetitions of the algorithm.

`[nb_clusters]`: The list with number of clusters. You can pass multiple values, separated by a single space.

`w`: A flag for saving the results in a file. Use 'w' if you wish to active this feature, or leave it blank. Important: the output file is saved in your current directory, within the folder `hgm_out`.

### Python
<!-- Firstly, you should have Cython installed. To install Cython, please refer to the official installation page:
https://cython.readthedocs.io/en/latest/src/quickstart/install.html -->

HG-means is also available as a Python package. To install HG-means, run the following installation command:

`> python -m pip install hgmeans`

**Important:** Check your user permissions for pip installation. You may need root user credentials.

For Windows users who do not have a C++ compiler, it may be required an installation of C++ Build tools, which can be downloaded here: https://go.microsoft.com/fwlink/?LinkId=691126

That is it! Now, open your Python interface, import the package and create an instance of HG-means. To execute it, just call the function `run()` with the corresponding parameters. See an example below:

`>>> import hgmeans`

`>>> demo = hgmeans.PyHGMeans()`

`>>> demo.run('data/iris.txt', 10, 5000, 1, [2,5,10], 'w')`

This script executes HG-means clustering for the "iris" dataset, with 10 solutions in population, a maximum of 5000 iterations, 1 repetition, and 2, 5 and 10 clusters. Here the number of clusters is passed in an array, so values are separated by commas.

## Data format

**Dataset files.** In the first line of a dataset file, the number of data points (n) and the dimensionality of the data (d) is set, separated by a single space. The remaining lines correspond to the coordinates of data points. Each line contains the values of the d features of a sample, where x_ij correspond to the j-th feature of the i-th sample of the data. Each feature value is separated by a single space, as depicted in the scheme below:

|  n   |   d  |      |     |      |
|------|------|------|-----|------|
| x_11 | x_12 | x_13 | ... | x_1d |
| x_21 | x_22 | x_23 | ... | x_2d |
| .... | .... | .... | ... | .... |
| x_n1 | x_n2 | x_n3 | ... | x_nd |

Some datasets are provided in `/data` folder in HG-means repository.

**Labels files.** The content of a labels file exhibits the cluster of each sample of the dataset according to ground-truth, where y_i correspond to the label of the i-th sample:

y_1

y_2

...

y_n

**Important**: Labels files must have the '.label' extension. Some labels are provided in `/data` folder in HG-means repository.
