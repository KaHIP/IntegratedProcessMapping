/******************************************************************************
 * integratedmapping.cpp 
 * 
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

        timer t;
        graph_io::readGraphWeighted(G, graph_filename);
        std::cout << "io time: " << t.elapsed()  << std::endl;

        G.set_partition_count(partition_config.k); 

        balance_configuration bc;
        bc.configurate_balance( partition_config, G);

        srand(partition_config.seed);
        random_functions::setSeed(partition_config.seed);
        if (partition_config.use_delta_gains) {
                partition_config.has_gains = new std::vector<bool>(G.number_of_nodes(),true);
                partition_config.delta = new std::vector<std::pair<int,std::vector<DELTA*>>>(G.number_of_nodes(),std::make_pair(0,std::vector<DELTA*>()));
                partition_config.ref_layer = new int();
                (*partition_config.ref_layer)=1;
        }
        
        std::cout <<  "graph has " <<  G.number_of_nodes() <<  " nodes and " <<  G.number_of_edges() <<  " edges"  << std::endl;
        
        matrix* D=NULL;
        std::vector< NodeID > *perm_rank = NULL;
        if (partition_config.enable_mapping || partition_config.integrated_mapping) {
                perm_rank = new std::vector< NodeID >(partition_config.k);
                for( unsigned i = 0; i < perm_rank->size(); i++) {
                        (*perm_rank)[i] = i;
                }
                partition_config.perm_rank = perm_rank;
        }
        bool power_of_two = (partition_config.k & (partition_config.k-1)) == 0;
        if (partition_config.enable_mapping || partition_config.integrated_mapping){
                if (power_of_two && !partition_config.enable_mapping && !partition_config.multisection) {
                        partition_config.construction_algorithm = MAP_CONST_IDENTITY;
                }

                if (partition_config.use_bin_id) {
                        D = new online_precalc_matrix(partition_config.k, partition_config.k);
                        D->setPartitionConfig(partition_config);
                } else if (partition_config.use_compact_bin_id) {
                        D = new online_binary_matrix(partition_config.k, partition_config.k);
                        D->setPartitionConfig(partition_config);
                } else if (partition_config.full_matrix) {
                        D = new full_matrix(partition_config.k, partition_config.k);
                        D->setPartitionConfig(partition_config);
                } else  
                        if( partition_config.distance_construction_algorithm != DIST_CONST_HIERARCHY_ONLINE) {
                                D = new normal_matrix(partition_config.k, partition_config.k);
                        } else {
                                D = new online_distance_matrix(partition_config.k, partition_config.k);
                                D->setPartitionConfig(partition_config);
                        }
                partition_config.D = D;
        }

        // ***************************** perform partitioning ***************************************       
        t.restart();
        graph_partitioner partitioner;
        quality_metrics qm;

        std::cout <<  "performing integrated mapping!"  << std::endl;
        if(partition_config.time_limit == 0) {
                if (partition_config.global_msec) {
                        partitioner.perform_partitioning_krec_hierarchy(partition_config, G);
                } else {
                        partitioner.perform_partitioning(partition_config, G);
                }
        } else {
                PartitionID* map = new PartitionID[G.number_of_nodes()];
                EdgeWeight best_cut = std::numeric_limits<EdgeWeight>::max();
                while(t.elapsed() < partition_config.time_limit) {
                        partition_config.graph_allready_partitioned = false;
                        if (partition_config.global_msec) {
                                partitioner.perform_partitioning_krec_hierarchy(partition_config, G);
                        } else {
                                partitioner.perform_partitioning(partition_config, G);
                        }
                        EdgeWeight cut = qm.edge_cut(G);
                        if(cut < best_cut) {
                                best_cut = cut;
                                forall_nodes(G, node) {
                                        map[node] = G.getPartitionIndex(node);
                                } endfor
                        }
                }

                forall_nodes(G, node) {
                        G.setPartitionIndex(node, map[node]);
                } endfor
        }

        if( partition_config.kaffpa_perfectly_balance && !partition_config.integrated_mapping) {
                double epsilon                         = partition_config.imbalance/100.0;
                partition_config.upper_bound_partition = (1+epsilon)*ceil(partition_config.largest_graph_weight/(double)partition_config.k);

                complete_boundary boundary(&G);
                boundary.build();

                cycle_refinement cr;
                cr.perform_refinement(partition_config, G, boundary);
        }
        ofs.close();
        std::cout.rdbuf(backup);
        double runtime = t.elapsed();

        NodeWeight qap = 0;
        if(partition_config.enable_mapping || partition_config.integrated_mapping) {
                int qap_tmp = 0;
                //check if k is a power of 2 
                bool power_of_two = (partition_config.k & (partition_config.k-1)) == 0;
                std::vector< NodeID > *perm_rank_tmp = new std::vector< NodeID >(partition_config.k);
                graph_access C;
                complete_boundary boundary(&G);
                boundary.build();
                boundary.getUnderlyingQuotientGraph(C);

                forall_nodes(C, node) {
                        C.setNodeWeight(node, 1);
                } endfor


                if (partition_config.integrated_mapping) {
                        qap = qm.total_qap_noquotient(G, *D, *perm_rank );
                }


                t.restart();

                if (partition_config.enable_mapping) {
                        
                        mapping_algorithms ma;
                        if( partition_config.distance_construction_algorithm != DIST_CONST_HIERARCHY_ONLINE) {
                                ma.construct_a_mapping(partition_config, C, *D, *perm_rank_tmp);
                        } else {
                                ma.construct_a_mapping(partition_config, C, *D, *perm_rank_tmp);
                        }
                       
                }


                if (partition_config.integrated_mapping) {
                    

                        refinement* refine = new mixed_refinement();
                        complete_boundary* coarser_boundary = new complete_boundary(&G);
                        coarser_boundary->build();
                        refine->perform_refinement(partition_config, G, *coarser_boundary);
                } else {
                        qap = qap_tmp;
                        for( unsigned i = 0; i < perm_rank->size(); i++) {
                                (*perm_rank)[i] = (*perm_rank_tmp)[i];
                        }
                }

                runtime += t.elapsed();

                std::cout <<  "time spent for integrated mapping " << runtime << std::endl;

                if (partition_config.integrated_mapping) {
                        qap = qm.total_qap_noquotient(G, *D, *perm_rank );
                } else {
                        qap = qm.total_qap(C, *D, *perm_rank );
                }

                forall_nodes(G, node) {
                        G.setPartitionIndex(node, (*perm_rank)[G.getPartitionIndex(node)]);
                } endfor
        }
        // ******************************* done partitioning *****************************************       
        // output some information about the partition that we have computed 
        std::cout << "cut \t\t"         << qm.edge_cut(G)                 << std::endl;
        std::cout << "finalobjective  " << qm.edge_cut(G)                 << std::endl;
        std::cout << "bnd \t\t"         << qm.boundary_nodes(G)           << std::endl;
        std::cout << "balance \t"       << qm.balance(G)                  << std::endl;
        std::cout << "max_comm_vol \t"  << qm.max_communication_volume(G) << std::endl;
        if(partition_config.enable_mapping || partition_config.integrated_mapping) 
                std::cout <<  "quadratic assignment objective J(C,D,Pi') = " << qap << std::endl;

        // write the partition to the disc 
        std::stringstream filename;
        if(!partition_config.filename_output.compare("")) {
                filename << "tmppartition" << partition_config.k;
        } else {
                filename << partition_config.filename_output;
        }

        if (!partition_config.suppress_output) {
                graph_io::writePartition(G, filename.str());
        } else {
                std::cout << "No partition will be written as output." << std::endl;
        }

}
