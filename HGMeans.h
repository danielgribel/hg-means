#ifndef HGMeans_H
#define HGMeans_H

/* Authors: Daniel Gribel and Thibaut Vidal
 * Contact: dgribel@inf.puc-rio.br
 *
 */

#include "GeneticOperations.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>

#define SAVE_FILE false

class InputValidator {
    
    private:

        int nb_arg;

        string dataset_path;

        string dataset_name;

        string label_path;

        int pi_min;

        int max_it;

        vector<int> nb_clusters;

    public:

        InputValidator(std::string filename, int size_population, int it, std::vector<int> m) {

            nb_arg = 3 + m.size();

            // Dataset path
            dataset_path = filename;

            // Population size (Pi_min paramter)
            pi_min = size_population;

            // Maxixum number of iterations parameter
            max_it = it;

            // List with number of clusters
            nb_clusters = m;

            // Find the position of last slash in the path
            unsigned slash_pos = dataset_path.find_last_of("/\\");
            
            // Find the position of last dot in the path
            unsigned point_pos = dataset_path.find_last_of(".");
            
            // Builds string of partial path
            string partial_path = dataset_path.substr(0, slash_pos + 1);
            
            // Builds string of dataset name
            dataset_name = dataset_path.substr(slash_pos + 1, point_pos - slash_pos - 1);
            
            // Builds string of labels path
            label_path = partial_path + dataset_name + ".label";
        }

        string GetDatasetPath() { return dataset_path; };

        string GetDatasetName() { return dataset_name; };

        string GetLabelPath() { return label_path; };

        int GetPiMin() { return pi_min; };

        int GetMaxIt() { return max_it; };

        vector<int> GetNbClusters() { return nb_clusters; };

        bool Validate() {
            // At least 4 parameters should be provided:
            // Dataset, size of population, maximum iterations and number of clusters
            if(nb_arg < 4) {
                cerr << "Insufficient number of parameters provided. Please use the following input format:" << 
                endl << "./hgmeans DatasetPath Pi_min N2 [M]" << endl;
                return false;
            }
            
            ifstream input(dataset_path);
            if (!input) {
                cerr << "Unable to open data file: " << dataset_path << endl;
                return false;
            }

            // Checking bounds of variables
            if(pi_min < 1 || pi_min > 100000) {
                cerr << "Pi_min out of bounds" << endl;
                return false;
            }
            
            if(max_it < 1 || max_it > 100000) {
                cerr << "Max_it out of bounds" << endl;
                return false;
            }
                
            for(int i = 0; i < nb_clusters.size(); i++) {
                if(nb_clusters[i] < 1 || nb_clusters[i] > USHRT_MAX) { // Range for m must be [1, 65535]
                    cerr << "The number of clusters is out of the limit [1, " << USHRT_MAX << "]" << endl;
                    return false;
                }
            }
            return true;
        }
};

class HGMeans {
    public:
        HGMeans();
        ~HGMeans();
        void Go(char* filename, int size_population, int it, const std::vector<int>& m); 
};

#endif