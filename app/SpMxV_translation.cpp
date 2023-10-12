/******************************************************************************
 * SpMxV_translation.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Marcelo Fonseca Faraj <marcelofaraj@gmail.com>
 *****************************************************************************/

#include <argtable3.h>
#include <iostream>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h> 
#include <vector>

#include "balance_configuration.h"
#include "data_structure/graph_access.h"
#include "data_structure/matrix/normal_matrix.h"
#include "data_structure/matrix/online_distance_matrix.h"
#include "data_structure/matrix/online_precalc_matrix.h"
#include "data_structure/matrix/online_binary_matrix.h"
#include "data_structure/matrix/full_matrix.h"
#include "graph_io.h"
#include "macros_assertions.h"
#include "mapping/mapping_algorithms.h"
#include "parse_parameters.h"
#include "partition/graph_partitioner.h"
#include "partition/partition_config.h"
#include "partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"
#include "partition/uncoarsening/refinement/mixed_refinement.h"
#include "partition/uncoarsening/refinement/label_propagation_refinement/label_propagation_refinement.h"

int main(int argn, char **argv) {

        PartitionConfig partition_config;
        std::string graph_filename;

        bool is_graph_weighted = false;
        bool suppress_output   = false;
        bool recursive         = false;

        int ret_code = parse_parameters(argn, argv, 
                                        partition_config, 
                                        graph_filename, 
                                        is_graph_weighted, 
                                        suppress_output, recursive); 

        if(ret_code) {
                return 0;
        }

        std::streambuf* backup = std::cout.rdbuf();
        std::ofstream ofs;
        ofs.open("/dev/null");
        if(suppress_output) {
                std::cout.rdbuf(ofs.rdbuf()); 
        }

        partition_config.LogDump(stdout);
        graph_access G;     

	timer processing_t;
        timer t;
        graph_io::readMatrixToGraph(partition_config, G, graph_filename);
        std::cout << "io time: " << t.elapsed()  << std::endl;


        std::vector<PartitionID> input_partition;
        if(partition_config.input_partition != "") {
                std::cout <<  "reading input partition" << std::endl;
                graph_io::readPartition(G, partition_config.input_partition);

		// write the partition to the disc 
		std::stringstream filename;
		if(!partition_config.filename_output.compare("")) {
			filename << graph_filename << partition_config.k;
		} else {
			filename << partition_config.filename_output;
		}
		graph_io::writeSpMxVPartition(partition_config, G, filename.str());

        } else {
		std::stringstream filename;
		if(!partition_config.filename_output.compare("")) {
			filename << graph_filename << ".graph";
		} else {
			filename << partition_config.filename_output << ".graph";
		}
		graph_io::writeGraph(G, filename.str());
	}

	return 0;
}
