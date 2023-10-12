/******************************************************************************
 * wcycle_partitioner.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/



#include "coarsening/coarsening_configurator.h"
#include "coarsening/contraction.h"
#include "coarsening/edge_rating/edge_ratings.h"
#include "coarsening/matching/gpa/gpa_matching.h"
#include "coarsening/matching/random_matching.h"
#include "coarsening/stop_rules/stop_rules.h"
#include "data_structure/graph_hierarchy.h"
#include "definitions.h"
#include "graph_partition_assertions.h"
#include "initial_partitioning/initial_partitioning.h"
#include "misc.h"
#include "random_functions.h"
#include "uncoarsening/refinement/mixed_refinement.h"
#include "uncoarsening/refinement/label_propagation_refinement/label_propagation_refinement.h"
#include "uncoarsening/refinement/refinement.h"
#include "wcycle_partitioner.h"
#include "mapping/mapping_algorithms.h"
#include "data_structure/matrix/online_distance_matrix.h"
#include "data_structure/matrix/online_precalc_matrix.h"
#include "data_structure/matrix/online_binary_matrix.h"
#include "partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.h"



int wcycle_partitioner::perform_partitioning(const PartitionConfig & config, graph_access & G) {
        PartitionConfig  cfg = config; 

        if(config.stop_rule == STOP_RULE_SIMPLE) {
                m_coarsening_stop_rule = new simple_stop_rule(cfg, G.number_of_nodes());
        } else {
                m_coarsening_stop_rule = new multiple_k_stop_rule(cfg, G.number_of_nodes());
        }
        int improvement = 0;
        improvement = (int) perform_partitioning_recursive(cfg, G, NULL);

        delete m_coarsening_stop_rule;

        return improvement;
}

int wcycle_partitioner::perform_partitioning_recursive( PartitionConfig & partition_config, 
                                                        graph_access & G, 
                                                        complete_boundary ** c_boundary) {

        //if graph not small enough
        //      perform matching two times
        //      perform coarsening two times
        //      call rekursive
        //else 
        //      initial partitioning
        //
        //refinement
        NodeID no_of_coarser_vertices = G.number_of_nodes();
        NodeID no_of_finer_vertices   = G.number_of_nodes();
        int improvement = 0;

        edge_ratings rating(partition_config);
        CoarseMapping* coarse_mapping =  new CoarseMapping();

        graph_access* finer                      = &G;
        matching* edge_matcher                   = NULL;
        contraction* contracter                  = new contraction();
        PartitionConfig copy_of_partition_config = partition_config;
        graph_access* coarser                    = new graph_access();

        Matching edge_matching;
        NodePermutationMap permutation;

        coarsening_configurator coarsening_config;
        coarsening_config.configure_coarsening(partition_config, &edge_matcher, m_level);
        
        rating.rate(*finer, m_level);

        edge_matcher->match(partition_config, *finer, edge_matching, *coarse_mapping, no_of_coarser_vertices, permutation);
        delete edge_matcher; 

        if(partition_config.graph_allready_partitioned) {
                contracter->contract_partitioned(partition_config, *finer, 
                                                 *coarser, edge_matching, 
                                                 *coarse_mapping, no_of_coarser_vertices, 
                                                 permutation);
        } else {
                contracter->contract(partition_config, *finer, 
                                     *coarser, edge_matching, 
                                     *coarse_mapping, no_of_coarser_vertices, 
                                     permutation);
        }

        coarser->set_partition_count(partition_config.k);
        complete_boundary* coarser_boundary =  NULL;
        refinement* refine = NULL;

        if(!partition_config.label_propagation_refinement || (partition_config.integrated_mapping && 
                                partition_config.initial_partitioning_type == INITIAL_PARTITIONING_RECPARTITION) ) {
                coarser_boundary = new complete_boundary(coarser);
                refine = new mixed_refinement();
        } else {
                refine = new label_propagation_refinement();
        }

        // Variables used for integrated_mapping
        std::vector< NodeID >* perm_rank = NULL;
        matrix *D = NULL;
        quality_metrics qm;
        int qap = 0;
        bool power_of_two = (partition_config.k & (partition_config.k-1)) == 0;
        if (partition_config.integrated_mapping) {
                perm_rank = partition_config.perm_rank;
                D         = partition_config.D;
        }
        //if (partition_config.integrated_mapping && partition_config.perm_rank == NULL) {
        //        perm_rank = new std::vector< NodeID >(partition_config.k);
        //        for( unsigned i = 0; i < perm_rank->size(); i++) {
        //                (*perm_rank)[i] = i;
        //        }
        //        partition_config.perm_rank = perm_rank;
        //}        
        // if (partition_config.integrated_mapping && D == NULL){
        //         if(!power_of_two !! partition_config.multisection) {
        //                 if( partition_config.distance_construction_algorithm != DIST_CONST_HIERARCHY_ONLINE) {
        //                         D = new normal_matrix(partition_config.k, partition_config.k);
        //                 } else {
        //                         D = new online_distance_matrix(partition_config.k, partition_config.k);
        //                         D->setPartitionConfig(partition_config);
        //                 }
        //         } else {
        //                 D = new online_distance_matrix(partition_config.k, partition_config.k);
        //                 D->setPartitionConfig(partition_config);
        //         }
        // }

        if(!m_coarsening_stop_rule->stop(no_of_finer_vertices, no_of_coarser_vertices)) {
                PartitionConfig cfg; cfg = partition_config;

                double factor = partition_config.balance_factor;
                cfg.upper_bound_partition = (factor +1.0)*partition_config.upper_bound_partition;

	        initial_partitioning init_part;

                // Initial partitioning followed by initial mapping (in coarsest level)
                if(cfg.integrated_mapping) {
                        init_part.perform_initial_partitioning(cfg, *coarser);

                        // qap = qm.total_qap_noquotient(*coarser, *D, *perm_rank );
                        // std::cout <<  "BEFORE - quadratic assignment objective J(C,D,Pi') = " << qap << std::endl;
                        //improvement += refine->perform_refinement_mapping(cfg, *coarser, *D, *perm_rank);

                        coarser_boundary->build();
                        improvement += refine->perform_refinement(cfg, *coarser, *coarser_boundary);

                        // qap = qm.total_qap_noquotient(*coarser, *D, *perm_rank );
                        // std::cout <<  "AFTER  - quadratic assignment objective J(C,D,Pi') = " << qap << std::endl;
                } else {
                        init_part.perform_initial_partitioning(cfg, *coarser);
                        if(!partition_config.label_propagation_refinement) 
                                coarser_boundary->build();
                        improvement += refine->perform_refinement(cfg, *coarser, *coarser_boundary);
                }

                m_deepest_level = m_level + 1;
        } else {
                m_level++;

                improvement += perform_partitioning_recursive( partition_config, *coarser, &coarser_boundary); 
                partition_config.graph_allready_partitioned = true;

                if(m_level % partition_config.level_split == 0 ) {

                        if(!partition_config.use_fullmultigrid 
                        || m_have_been_level_down.find(m_level) == m_have_been_level_down.end())  { 

                                if(!partition_config.label_propagation_refinement) {
                                        delete coarser_boundary;

                                        coarser_boundary                = new complete_boundary(coarser);
                                }
                                m_have_been_level_down[m_level] = true;

                                // configurate the algorithm to use the same amount
                                // of imbalance as was allowed on this level 
                                PartitionConfig cfg;
                                cfg = partition_config;
                                cfg.set_upperbound = false;

                                double cur_factor = partition_config.balance_factor/(m_deepest_level-m_level);
                                cfg.upper_bound_partition = ( (m_level != 0) * cur_factor+1.0)*partition_config.upper_bound_partition;

                                // Integrated mapping in the v corners of w-cycle and f-cycle
                                if(cfg.integrated_mapping && m_level == 0) {
                                        graph_access C;
                                        complete_boundary boundary(coarser);
                                        boundary.build();
                                        boundary.getUnderlyingQuotientGraph(C);
                                        int qap_tmp = 0;

                                        std::vector< NodeID >* perm_rank_tmp = new std::vector< NodeID >(partition_config.k); 

                                        forall_nodes(C, node) {
                                                C.setNodeWeight(node, 1);
                                        } endfor

                                        // if(!power_of_two || cfg.multisection) {
                                                mapping_algorithms ma;
                                                if( cfg.distance_construction_algorithm != DIST_CONST_HIERARCHY_ONLINE) {
                                                        ma.construct_a_mapping(cfg, C, *D, *perm_rank_tmp);
                                                        qap_tmp = qm.total_qap(C, *D, *perm_rank_tmp );
                                                } else {
                                                        ma.construct_a_mapping(cfg, C, *D, *perm_rank_tmp);
                                                        qap_tmp = qm.total_qap(C, *D, *perm_rank_tmp );
                                                }
                                        // } else {
                                        //         for( unsigned i = 0; i < perm_rank_tmp->size(); i++) {
                                        //                 (*perm_rank_tmp)[i] = i;
                                        //         }
                                        //         qap_tmp = qm.total_qap(C, *D, *perm_rank_tmp );
                                        // }

                                        qap = qm.total_qap_noquotient(*coarser, *D, *perm_rank );
                                        if (qap_tmp < qap) {
                                                qap = qap_tmp;
                                                for( unsigned i = 0; i < perm_rank->size(); i++) {
                                                        (*perm_rank)[i] = (*perm_rank_tmp)[i];
                                                }
                                                //std::cout <<  "quadratic assignment objective J(C,D,Pi') = " << qap << std::endl;
                                        }
                                        delete perm_rank_tmp;
                                }

                                // do the next arm of the F-cycle
                                improvement += perform_partitioning_recursive( cfg, *coarser, &coarser_boundary); 
                        }
                }

                m_level--;

        }

        if(partition_config.use_balance_singletons && !partition_config.label_propagation_refinement) {
                coarser_boundary->balance_singletons( partition_config, *coarser );
        }
        
        //project
        graph_access& fRef = *finer;
        graph_access& cRef = *coarser;
        forall_nodes(fRef, n) {
                NodeID coarser_node              = (*coarse_mapping)[n];
                PartitionID coarser_partition_id = cRef.getPartitionIndex(coarser_node);
                fRef.setPartitionIndex(n, coarser_partition_id);
        } endfor

        finer->set_partition_count(coarser->get_partition_count());
        complete_boundary* current_boundary = NULL;
        if(!partition_config.label_propagation_refinement) {
                current_boundary = new complete_boundary(finer);
                current_boundary->build_from_coarser(coarser_boundary, no_of_coarser_vertices, coarse_mapping ); 
        }

        PartitionConfig cfg; cfg = partition_config;
        double cur_factor = partition_config.balance_factor/(m_deepest_level-m_level);

        //only set the upperbound if it is the first time 
        //we go down the F-cycle
        if( partition_config.set_upperbound ) {
                cfg.upper_bound_partition = ( (m_level != 0) * cur_factor+1.0)*partition_config.upper_bound_partition;
        } else {
                cfg.upper_bound_partition = partition_config.upper_bound_partition;
        }

        
        if(cfg.integrated_mapping) {
                // qap = qm.total_qap_noquotient(*finer, *D, *perm_rank );
                // std::cout <<  "BEFORE - quadratic assignment objective J(C,D,Pi') = " << qap << std::endl;
                
                improvement += refine->perform_refinement(cfg, *finer, *current_boundary);
                //improvement += refine->perform_refinement_mapping(cfg, *finer, *D, *perm_rank);
                // qap = qm.total_qap_noquotient(*finer, *D, *perm_rank );
                // std::cout <<  "AFTER  - quadratic assignment objective J(C,D,Pi') = " << qap << std::endl;
        } else {
                improvement += refine->perform_refinement(cfg, *finer, *current_boundary);
        }

        if(c_boundary != NULL) {
                delete *c_boundary;
                *c_boundary = current_boundary;
        } else {
		if( current_boundary != NULL ) delete current_boundary;
	}

        delete contracter;
        delete coarse_mapping;
        delete coarser_boundary;
        delete coarser;
        delete refine;

        return improvement;
}
