/******************************************************************************
 * label_propagation_refinement.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/


#include "label_propagation_refinement.h"
#include "partition/coarsening/clustering/node_ordering.h"
#include "tools/random_functions.h"
#include <map>
#include <set>



label_propagation_refinement::label_propagation_refinement() {
                
}

label_propagation_refinement::~label_propagation_refinement() {
                
}

EdgeWeight label_propagation_refinement::perform_refinement(PartitionConfig & partition_config, 
                                                            graph_access & G, 
                                                            complete_boundary & boundary) {
        NodeWeight block_upperbound = partition_config.upper_bound_partition;

        // in this case the _matching paramter is not used 
        // coarse_mappng stores cluster id and the mapping (it is identical)
        std::vector<PartitionID> hash_map(partition_config.k,0);
        std::vector<NodeID> permutation(G.number_of_nodes());
        std::vector<NodeWeight> cluster_sizes(partition_config.k, 0);

        node_ordering n_ordering;
        n_ordering.order_nodes(partition_config, G, permutation);

        std::queue< NodeID > * Q             = new std::queue< NodeID >();
        std::queue< NodeID > * next_Q        = new std::queue< NodeID >();
        std::vector<bool> * Q_contained      = new std::vector<bool>(G.number_of_nodes(), false);
        std::vector<bool> * next_Q_contained = new std::vector<bool> (G.number_of_nodes(), false);
        forall_nodes(G, node) {
                cluster_sizes[G.getPartitionIndex(node)] += G.getNodeWeight(node);
                Q->push(permutation[node]);
        } endfor

        for( int j = 0; j < partition_config.label_iterations_refinement; j++) {
                unsigned int change_counter = 0;
                while( !Q->empty() ) {
                        NodeID node = Q->front();
                        Q->pop();
                        (*Q_contained)[node] = false;

                        //now move the node to the cluster that is most common in the neighborhood
                        forall_out_edges(G, e, node) {
                                NodeID target = G.getEdgeTarget(e);
                                hash_map[G.getPartitionIndex(target)]+=G.getEdgeWeight(e);
                                //std::cout <<  "curblock " <<  G.getPartitionIndex(target)  << std::endl;
                        } endfor

                        //second sweep for finding max and resetting array
                        PartitionID max_block = G.getPartitionIndex(node);
                        PartitionID my_block  = G.getPartitionIndex(node);

                        PartitionID max_value = 0;
                        forall_out_edges(G, e, node) {
                                NodeID target             = G.getEdgeTarget(e);
                                PartitionID cur_block     = G.getPartitionIndex(target);
                                PartitionID cur_value     = hash_map[cur_block];
                                if((cur_value > max_value  || (cur_value == max_value && random_functions::nextBool())) 
                                && (cluster_sizes[cur_block] + G.getNodeWeight(node) < block_upperbound || (cur_block == my_block && cluster_sizes[my_block] <= partition_config.upper_bound_partition)))
                                {
                                        max_value = cur_value;
                                        max_block = cur_block;
                                }

                                hash_map[cur_block] = 0;
                        } endfor

                        cluster_sizes[G.getPartitionIndex(node)]  -= G.getNodeWeight(node);
                        cluster_sizes[max_block]         += G.getNodeWeight(node);
                        bool changed_label                = G.getPartitionIndex(node) != max_block; 
                        change_counter                   += changed_label;
                        G.setPartitionIndex(node, max_block);
                        //std::cout <<  "maxblock " <<  max_block  << std::endl;

                        if(changed_label) {
                                forall_out_edges(G, e, node) {
                                        NodeID target = G.getEdgeTarget(e);
                                        if(!(*next_Q_contained)[target]) {
                                                next_Q->push(target);
                                                (*next_Q_contained)[target] = true;
                                        } 
                                } endfor
                        }
                } 

                std::swap( Q, next_Q);
                std::swap( Q_contained, next_Q_contained);

        }
        

        delete Q;
        delete next_Q;
        delete Q_contained;
        delete next_Q_contained;
        
        return 0;

}



EdgeWeight label_propagation_refinement::perform_refinement_mapping(PartitionConfig & config, graph_access & G, matrix & D, 
                                                std::vector< NodeID > & perm_rank, complete_boundary & boundary){
        
        if (commons == NULL) {
                commons = new kway_graph_refinement_commons(config);
        }
        NodeWeight block_upperbound = config.upper_bound_partition;

        // in this case the _matching paramter is not used 
        // coarse_mappng stores cluster id and the mapping (it is identical)
        std::vector<NodeID> permutation(G.number_of_nodes());
        std::vector<NodeWeight> cluster_sizes(config.k, 0);

        node_ordering n_ordering;
        n_ordering.order_nodes(config, G, permutation);

        std::queue< NodeID > * Q             = new std::queue< NodeID >();
        std::queue< NodeID > * next_Q        = new std::queue< NodeID >();
        std::vector<bool> * Q_contained      = new std::vector<bool>(G.number_of_nodes(), false);
        std::vector<bool> * next_Q_contained = new std::vector<bool> (G.number_of_nodes(), false);
        forall_nodes(G, node) {
                Q->push(permutation[node]);
        } endfor

        for( int j = 0; j < config.label_iterations_refinement_map; j++) {
                unsigned int change_counter = 0;
                if ( Q->empty() ) {
                        forall_nodes(G, node) {
                                Q->push(permutation[node]);
                        } endfor
                }
                while( !Q->empty() ) {
                        NodeID my_node = Q->front();
                        Q->pop();
                        (*Q_contained)[my_node] = false;

                        if (!is_boundary(my_node,G)) {
                                continue;
                        }

                        PartitionID my_block   = G.getPartitionIndex(my_node);
                        PartitionID min_block  = my_block;
                        int min_value = 0;
                        NodeWeight my_weight = G.getNodeWeight(my_node);
                        EdgeWeight min_edge_degree = 0;

                        if (config.use_delta_gains) {
                                min_value = commons->best_delta_gain(config, boundary, G, my_node, min_block, min_edge_degree, false, true);
                                if (min_block == INVALID_PARTITION) {
                                        min_block  = my_block;
                                }
                        } else {
                                std::vector<PartitionID> hash_map(config.k,0);  // stores sum(edges) incident from my_node to each partition
                                std::set<PartitionID> disputing_partitions;
                                disputing_partitions.insert(my_block);
                                NodeID perm_rank_node  = perm_rank[my_block];
                                std::map<PartitionID, int> node_objective_function;
                                node_objective_function[my_block] = 0;
                                //now move my_node to the cluster that is most common in the neighborhood
                                forall_out_edges(G, e, my_node) {
                                        NodeID target           = G.getEdgeTarget(e);
                                        PartitionID cur_block   = G.getPartitionIndex(target);
                                        hash_map[cur_block]+=G.getEdgeWeight(e);
                                        if (cur_block != my_block) {
                                                disputing_partitions.insert(cur_block);
                                                NodeID perm_rank_target   = perm_rank[cur_block];
                                                // if (cluster_sizes[my_block] > block_upperbound) {
                                                if (boundary.getBlockWeight(my_block) > block_upperbound) {
                                                // We keep my_block as infinite but still among the possibilities in case no other move is feasible
                                                        node_objective_function[my_block] = std::numeric_limits<int>::max();
                                                }
                                                else {
                                                        node_objective_function[my_block] += G.getEdgeWeight(e) * D.get_xy(perm_rank_node, perm_rank_target);
                                                }
                                        }
                                } endfor

                                //second sweep for finding min and resetting array
                                min_value = node_objective_function[my_block];
                                forall_out_edges(G, e, my_node) {
                                        NodeID target             = G.getEdgeTarget(e);
                                        PartitionID cur_block     = G.getPartitionIndex(target);
                                        NodeID perm_rank_target   = perm_rank[cur_block];
                                        // Avoid useless calculations
                                        if ( node_objective_function.find(cur_block) != node_objective_function.end() ||
                                                ((boundary.getBlockWeight(cur_block) + my_weight > block_upperbound) 
                                        &&   (cur_block != my_block || boundary.getBlockWeight(my_block) > config.upper_bound_partition)) ) 
                                        {
                                                continue;
                                        }
                                        // Calculate partial mapping objective function in case 'my_node' goes to 'cur_block'
                                        node_objective_function[cur_block] = 0;
                                        std::set<PartitionID>::iterator it;
                                        for (it = disputing_partitions.begin(); it != disputing_partitions.end(); ++it)
                                        {
                                                NodeID cur_block_inner = *it; 
                                                if (cur_block_inner == cur_block)       
                                                        continue;
                                                PartitionID cur_value_inner          = hash_map[cur_block_inner];
                                                NodeID perm_rank_target_inner        = perm_rank[cur_block_inner];
                                                node_objective_function[cur_block]  += cur_value_inner * D.get_xy(perm_rank_target, perm_rank_target_inner);
                                        }
                                        if( node_objective_function[cur_block] < min_value  
                                                || (node_objective_function[cur_block] == min_value && random_functions::nextBool()))
                                        {
                                                min_value = node_objective_function[cur_block];
                                                min_block = cur_block;
                                        }
                                } endfor

                                // cluster_sizes[my_block]     -= my_weight;
                                // cluster_sizes[min_block]    += my_weight;
                        }
                        
                        bool changed_label           = my_block != min_block; 
                        change_counter              += changed_label;


                        if(changed_label) {
                                forall_out_edges(G, e, my_node) {
                                        NodeID target = G.getEdgeTarget(e);
                                        if(!(*next_Q_contained)[target]) {
                                                next_Q->push(target);
                                                (*next_Q_contained)[target] = true;
                                        } 
                                } endfor

                                G.setPartitionIndex(my_node, min_block);        

                                boundary_pair pair;
                                pair.k = config.k;
                                pair.lhs = my_block;
                                pair.rhs = min_block;

                                boundary.postMovedBoundaryNodeUpdates(my_node, &pair, true, true);

                                boundary.setBlockNoNodes(my_block, boundary.getBlockNoNodes(my_block)-1);
                                boundary.setBlockNoNodes(min_block,   boundary.getBlockNoNodes(min_block)+1);
                                boundary.setBlockWeight( my_block, boundary.getBlockWeight(my_block)-my_weight);
                                boundary.setBlockWeight( min_block,   boundary.getBlockWeight(min_block)+my_weight);

                                if (config.use_delta_gains) {
                                        commons->move_delta_gain (config, boundary, G, my_node, my_block, min_block);
                                        // commons->move_delta_degrees (config, boundary, G, my_node, my_block, min_block);       
                                }
                        }
                } 

                std::swap( Q, next_Q);
                std::swap( Q_contained, next_Q_contained);

        }

        delete Q;
        delete next_Q;
        delete Q_contained;
        delete next_Q_contained;
   
        return 0;
}




bool label_propagation_refinement::is_boundary(NodeID node, graph_access & G) {

        PartitionID my_part   = G.getPartitionIndex(node);

        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                if (G.getPartitionIndex(target) == my_part) {
                        return true;
                }
        } endfor
        return false;
}  