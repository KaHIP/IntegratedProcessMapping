/******************************************************************************
 * kway_graph_refinement_commons.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef KWAY_GRAPH_REFINEMENT_COMMONS_PVGY97EW
#define KWAY_GRAPH_REFINEMENT_COMMONS_PVGY97EW

#include <vector>

#include "data_structure/priority_queues/priority_queue_interface.h"
#include "definitions.h"
#include "random_functions.h"
#include "uncoarsening/refinement/refinement.h"
#include "uncoarsening/refinement/quotient_graph_refinement/2way_fm_refinement/vertex_moved_hashtable.h"
#include <unordered_map>
#include <set>
#include <unordered_set>
#include "quality_metrics.h"

class kway_graph_refinement_commons  {
        public:

                kway_graph_refinement_commons( PartitionConfig & config );
                virtual ~kway_graph_refinement_commons();

                void init( PartitionConfig & config );

                bool incident_to_more_than_two_partitions(graph_access & G, NodeID & node);

                EdgeWeight compute_gain(graph_access & G, 
                                        NodeID & node, 
                                        PartitionID & max_gainer, 
                                        EdgeWeight & ext_degree);

                EdgeWeight compute_gain_map(const PartitionConfig & config,
                                            complete_boundary & boundary,
                                            graph_access & G, 
                                            NodeID & node, 
                                            PartitionID & max_gainer, 
                                            EdgeWeight & max_degree,
                                            bool ignore_balance = true);

                EdgeWeight best_delta_gain(PartitionConfig & config,
                                           complete_boundary & boundary,
                                           graph_access & G, 
                                           NodeID & node, 
                                           PartitionID & max_gainer, 
                                           EdgeWeight & max_degree,
                                           bool ignore_balance = true,
                                           bool label_prop = false);

                EdgeWeight init_delta_gain(PartitionConfig & config,
                                           complete_boundary & boundary,
                                           graph_access & G, 
                                           NodeID & node, 
                                           PartitionID & max_gainer, 
                                           EdgeWeight & max_degree,
                                           bool ignore_balance = true,
                                           bool label_prop = false);  

                void move_delta_gain(PartitionConfig & config,
                                     complete_boundary & boundary,
                                     graph_access & G,
                                     NodeID node,
                                     const PartitionID par_old,
                                     const PartitionID par_new);  

                int par_new_delta_gain(PartitionConfig & config,
                                       NodeID node,
                                       NodeID neighb,
                                       graph_access & G);

                inline int init_delta_degrees(PartitionConfig & config,
                                              graph_access & G, 
                                              NodeID & node,
                                              const PartitionID lhs,
                                              const PartitionID rhs,
                                              EdgeWeight & int_obj_func,
                                              EdgeWeight & ext_obj_func,
                                              bool & has_ext_edge);

                inline void move_delta_degrees (PartitionConfig & config,
                                                complete_boundary & boundary,
                                                graph_access & G,
                                                NodeID node,
                                                const PartitionID par_old,
                                                const PartitionID par_new);
                
                inline int init_gain_from_degree (PartitionConfig & config,
                                                  complete_boundary & boundary,
                                                  graph_access & G, 
                                                  NodeID & node, 
                                                  PartitionID & max_gainer, 
                                                  EdgeWeight & max_degree,
                                                  bool ignore_balance = true,
                                                  bool label_prop = false);

                inline void ignore_delta_gain(PartitionConfig & config,
                                              graph_access & G, 
                                              NodeID & node, 
                                              int increment);

                bool int_ext_degree(graph_access & G, 
                                    const NodeID & node,
                                    const PartitionID lhs,
                                    const PartitionID rhs,
                                    EdgeWeight & int_degree,
                                    EdgeWeight & ext_degree);
                

                inline unsigned getUnderlyingK();

                inline void erase_block_delta(std::vector<DELTA*> &delta_node, size_t position);

        private:

                //for efficient computation of internal and external degrees
                struct round_struct {
                        unsigned round;
                        EdgeWeight local_degree;
                };

                std::vector<round_struct>                    m_local_degrees;
                unsigned                                     m_round;
};

inline unsigned kway_graph_refinement_commons::getUnderlyingK() {
        return m_local_degrees.size();
}

inline void kway_graph_refinement_commons::init(PartitionConfig & config) {
        m_local_degrees.resize(config.k);
        for( PartitionID i = 0; i < config.k; i++) {
                m_local_degrees[i].round        = 0;
                m_local_degrees[i].local_degree = 0;
        }

        m_round = 0;//needed for the computation of internal and external degrees
}

inline bool kway_graph_refinement_commons::incident_to_more_than_two_partitions(graph_access & G, NodeID & node) {
        bool ret_value = false;
        PartitionID own_partition = G.getPartitionIndex(node);
        PartitionID second_partition = INVALID_PARTITION;

        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                PartitionID target_partition = G.getPartitionIndex(target);
                if(target_partition != own_partition) {
                        if(second_partition == INVALID_PARTITION) {
                                second_partition = target_partition;
                        } else if(target_partition != second_partition) {
                                ret_value = true;
                                break;
                        }
                }

        } endfor

        return ret_value;
}

inline bool kway_graph_refinement_commons::int_ext_degree( graph_access & G, 
                                                           const NodeID & node,
                                                           const PartitionID lhs,
                                                           const PartitionID rhs,
                                                           EdgeWeight & int_degree,
                                                           EdgeWeight & ext_degree) {


        ASSERT_EQ(lhs, G.getPartitionIndex(node));

        int_degree               = 0;
        ext_degree               = 0;
        bool update_is_difficult = false;

        forall_out_edges(G, e, node) {
                NodeID target                 = G.getEdgeTarget(e);
                PartitionID targets_partition = G.getPartitionIndex(target);

                if(targets_partition == lhs) {
                        int_degree += G.getEdgeWeight(e); 
                } else if(targets_partition == rhs) {
                        ext_degree += G.getEdgeWeight(e);
                }

                if(targets_partition != lhs && targets_partition != rhs) {
                        update_is_difficult = true;
                } 
        } endfor

        return update_is_difficult;
}

inline Gain kway_graph_refinement_commons::compute_gain(graph_access & G, 
                                                        NodeID & node, 
                                                        PartitionID & max_gainer, 
                                                        EdgeWeight & ext_degree) {
        //for all incident partitions compute gain
        //return max gain and max_gainer partition
        PartitionID source_partition = G.getPartitionIndex(node);
        EdgeWeight max_degree        = 0;
        max_gainer                   = INVALID_PARTITION;

        m_round++;//can become zero again
        forall_out_edges(G, e, node) {
                NodeID target                = G.getEdgeTarget(e);
                PartitionID target_partition = G.getPartitionIndex(target);

                if(m_local_degrees[target_partition].round == m_round) {
                        m_local_degrees[target_partition].local_degree += G.getEdgeWeight(e);
                } else {
                        m_local_degrees[target_partition].local_degree = G.getEdgeWeight(e);
                        m_local_degrees[target_partition].round = m_round;
                }


                if(m_local_degrees[target_partition].local_degree >= max_degree && target_partition != source_partition) {
                        if(m_local_degrees[target_partition].local_degree > max_degree) {
                                max_degree = m_local_degrees[target_partition].local_degree;
                                max_gainer = target_partition;
                        } else {
                                //break ties randomly
                                bool accept = random_functions::nextBool();
                                if(accept) {
                                        max_degree = m_local_degrees[target_partition].local_degree;
                                        max_gainer = target_partition;
                                }
                        }
                }
        } endfor

        if(max_gainer != INVALID_PARTITION) {
                ext_degree = max_degree;
        } else {
                ext_degree = 0;
        }

        if(m_local_degrees[source_partition].round != m_round) {
                m_local_degrees[source_partition].local_degree = 0;
        } 

        return max_degree-m_local_degrees[source_partition].local_degree;
}



 
inline Gain kway_graph_refinement_commons::compute_gain_map (const PartitionConfig & config,
                                                             complete_boundary & boundary,
                                                             graph_access & G, 
                                                             NodeID & node, 
                                                             PartitionID & max_gainer, 
                                                             EdgeWeight & max_degree,
                                                             bool ignore_balance = true) {

        //for all incident partitions compute gain
        //return max gain and max_gainer partition
        PartitionID source_partition = G.getPartitionIndex(node);
        EdgeWeight min_value         = 0;
        max_gainer                = INVALID_PARTITION;

        std::vector< NodeID >* perm_rank = config.perm_rank;
        matrix* D = config.D;
        std::unordered_map<PartitionID, int> choice_objective_function;
        int & objective_source = choice_objective_function[source_partition];
        objective_source = 0;
        NodeID perm_rank_source                     = (*perm_rank)[source_partition];
        std::set<PartitionID> disputing_partitions;
        disputing_partitions.insert(source_partition);


        m_round++;//can become zero again

        m_local_degrees[source_partition].local_degree = 0;
        m_local_degrees[source_partition].round = m_round;

        // m_local_degrees stores sum(edges) incident from my_node to each partition
        forall_out_edges(G, e, node) {
                NodeID target                = G.getEdgeTarget(e);
                PartitionID target_partition = G.getPartitionIndex(target);
                NodeID perm_rank_target   = (*perm_rank)[target_partition];

                if(m_local_degrees[target_partition].round == m_round) {
                        m_local_degrees[target_partition].local_degree += G.getEdgeWeight(e);
                } else {
                        m_local_degrees[target_partition].local_degree = G.getEdgeWeight(e);
                        m_local_degrees[target_partition].round = m_round;
                }

                if (target_partition != source_partition) {
                        disputing_partitions.insert(target_partition);
                        objective_source += G.getEdgeWeight(e) * D->get_xy(perm_rank_source, perm_rank_target);
                }
        } endfor

        // Second sweep for finding min
        NodeWeight this_nodes_weight = G.getNodeWeight(node);
        // int min_value_tmp = objective_source;
        int min_value_tmp = std::numeric_limits<int>::max();

        // if (boundary.getBlockWeight(source_partition) > config.upper_bound_partition) {
        //         // Keep my_block as infinite but still among the possibilities in case no other move is feasible
        //         min_value_tmp = std::numeric_limits<int>::max();
        // }
        forall_out_edges(G, e, node) {
                NodeID target                = G.getEdgeTarget(e);
                PartitionID target_partition = G.getPartitionIndex(target);
                NodeID perm_rank_target   = (*perm_rank)[target_partition];

                if ( choice_objective_function.find(target_partition) != choice_objective_function.end() ||
                        // boundary.getBlockWeight(target_partition) + this_nodes_weight > config.upper_bound_partition ||
                        target_partition == source_partition) {
                        continue;
                }

                // Calculate partial mapping objective function in case 'source_partition' goes to 'target_partition'
                choice_objective_function[target_partition] = 0;
                std::set<PartitionID>::iterator it;
                for (it = disputing_partitions.begin(); it != disputing_partitions.end(); ++it)
                {
                        NodeID target_partition_inner = *it; 
                        if (target_partition_inner == target_partition)       
                                continue;
                        PartitionID target_value_inner       = m_local_degrees[target_partition_inner].local_degree;
                        NodeID perm_rank_target_inner        = (*perm_rank)[target_partition_inner];
                        choice_objective_function[target_partition]  += target_value_inner * 
                                                                        D->get_xy(perm_rank_target, perm_rank_target_inner);
                }
                if( (ignore_balance || boundary.getBlockWeight(target_partition) + this_nodes_weight <= config.upper_bound_partition) 
                  && ( choice_objective_function[target_partition] < min_value_tmp  
                  || (choice_objective_function[target_partition] == min_value_tmp && random_functions::nextBool()) ))// pseudorandom
                {
                        min_value_tmp = choice_objective_function[target_partition];
                        min_value     = min_value_tmp;
                        max_gainer = target_partition;
                }
        } endfor

        if(max_gainer != INVALID_PARTITION) {
                max_degree = m_local_degrees[max_gainer].local_degree;
        } else {
                max_degree = 0;
                min_value = 0;
                objective_source = 0;
        }

        if(m_local_degrees[source_partition].round != m_round) {
                m_local_degrees[source_partition].local_degree = 0;
        } 

        return (Gain)objective_source - min_value;
}




inline Gain kway_graph_refinement_commons::best_delta_gain (PartitionConfig & config,
                                                            complete_boundary & boundary,
                                                            graph_access & G, 
                                                            NodeID & node, 
                                                            PartitionID & max_gainer, 
                                                            EdgeWeight & max_degree,
                                                            bool ignore_balance = true,
                                                            bool label_prop = false) {
        if (config.skip_delta_gains) {
                return compute_gain_map (config, boundary, G, node, max_gainer, max_degree,ignore_balance);
        }


        if ((*config.delta)[node].first != (*config.ref_layer)) {
                return init_delta_gain (config, boundary, G, node, max_gainer, max_degree, ignore_balance, label_prop);
        } else if (!(*config.has_gains)[node]) {
                return init_gain_from_degree (config, boundary, G, node, max_gainer, max_degree, ignore_balance, label_prop);
        }

        //for all incident partitions compute gain
        //return max gain and max_gainer partition
        PartitionID source_partition = G.getPartitionIndex(node);
        max_gainer                   = INVALID_PARTITION;
        NodeWeight this_nodes_weight = G.getNodeWeight(node);
        int max_gain                 = (!label_prop || boundary.getBlockWeight(source_partition) > config.upper_bound_partition) 
                                        * std::numeric_limits<int>::min();
        int max_gain_final           = 0;
        max_degree                   = 0;

        std::vector<DELTA*> &delta_node = (*config.delta)[node].second;
        std::vector<DELTA*>::iterator it;
        for (it = delta_node.begin(); it != delta_node.end(); ++it) {
                DELTA* curr_delta = (*it);
                PartitionID p = curr_delta->block;
                if ( (!ignore_balance && boundary.getBlockWeight(p) + this_nodes_weight > config.upper_bound_partition )
                    || p == source_partition) {
                        continue;
                }

                int gain_p = curr_delta->gain;
                if( gain_p > max_gain || (gain_p == max_gain && random_functions::nextBool()) ) {
                        max_gain   = gain_p;
                        max_gain_final = max_gain;
                        max_degree = curr_delta->degree;
                        max_gainer = p;
                }
        }

        return (Gain)max_gain_final;
}







// for all incident partitions in 'node', compute delta gain from scratch
inline EdgeWeight kway_graph_refinement_commons::init_delta_gain (PartitionConfig & config,
                                                                  complete_boundary & boundary,
                                                                  graph_access & G, 
                                                                  NodeID & node, 
                                                                  PartitionID & max_gainer, 
                                                                  EdgeWeight & max_degree,
                                                                  bool ignore_balance = true,
                                                                  bool label_prop=false) {
        
        PartitionID source_partition = G.getPartitionIndex(node);
        std::vector< NodeID >* perm_rank = config.perm_rank;
        matrix* D = config.D;
        NodeID perm_rank_source = (*perm_rank)[source_partition];

        std::vector<DELTA*> &delta_node = (*config.delta)[node].second;
        // to make it consistent
        for (DELTA* obj : delta_node) {
                delete obj;
        }
        delta_node.clear();
        
        (*config.delta)[node].first = (*config.ref_layer);
        
        // int max_gain = (label_prop)? 0 : std::numeric_limits<int>::min();
        int max_gain = (!label_prop || boundary.getBlockWeight(source_partition) > config.upper_bound_partition) 
                       * std::numeric_limits<int>::min();
        int max_gain_final = 0;

        max_degree = 0;
        max_gainer = INVALID_PARTITION;

        
        DELTA* delta_node_source = new DELTA {
                source_partition,       // block
                0,                      // gain
                0                       // degree
        };
        delta_node.push_back(delta_node_source);
        std::unordered_map<PartitionID,DELTA*> stored;
        stored[source_partition] = delta_node_source;

        int curr_OF = 0;
        NodeWeight this_nodes_weight = G.getNodeWeight(node);

        // config.delta[node] stores sum(weight(edge from 'node')) incident in each partition  
        forall_out_edges(G, e, node) {
                NodeID target                = G.getEdgeTarget(e);
                PartitionID target_partition = G.getPartitionIndex(target);
                NodeID perm_rank_target   = (*perm_rank)[target_partition];

                if (stored.find(target_partition) == stored.end()) {
                        DELTA * delta_node_target = new DELTA {
                                target_partition,       // block
                                0,                      // gain
                                G.getEdgeWeight(e)      // degree
                        };      
                        stored[target_partition] = delta_node_target;
                        delta_node.push_back(delta_node_target);
                } else {
                        stored[target_partition]->degree += G.getEdgeWeight(e);
                }

                if (target_partition != source_partition) {
                        curr_OF  += G.getEdgeWeight(e) * D->get_xy(perm_rank_source, perm_rank_target);
                }
        } endfor
        (*config.has_gains)[node] = true;

        // std::set<PartitionID> checked;
        
        // forall_out_edges(G, e, node) {
        //         NodeID target                = G.getEdgeTarget(e);
        //         PartitionID target_partition = G.getPartitionIndex(target);
        //         if (target_partition == source_partition || 
        //                 // boundary.getBlockWeight(target_partition) + this_nodes_weight > config.upper_bound_partition ||
        //                 checked.find(target_partition) != checked.end()) {
        //                 continue;
        //         }
        //         checked.insert(target_partition);
        std::vector<DELTA*>::iterator it;
        for (it = delta_node.begin(); it != delta_node.end(); ++it) {
                DELTA* delta_node_target = (*it);
                PartitionID target_partition = delta_node_target->block;
                NodeID perm_rank_target   = (*perm_rank)[target_partition];

                if (target_partition == source_partition) {
                        continue;
                }

                int target_gain = curr_OF;
                std::vector<DELTA*>::iterator it_inner;
                for (it_inner = delta_node.begin(); it_inner != delta_node.end(); ++it_inner)
                {
                        DELTA* delta_node_target_inner = (*it_inner);
                        PartitionID target_partition_inner = delta_node_target_inner->block; 
                        if (target_partition_inner == target_partition)       
                                continue;
                        PartitionID target_value_inner       = delta_node_target_inner->degree;
                        NodeID perm_rank_target_inner        = (*perm_rank)[target_partition_inner];
                        target_gain  -= target_value_inner * D->get_xy(perm_rank_target, perm_rank_target_inner);
                }
                delta_node_target->gain = target_gain;
                if( (ignore_balance||boundary.getBlockWeight(target_partition)+this_nodes_weight<=config.upper_bound_partition)
                 && ( target_gain > max_gain || (target_gain == max_gain && random_functions::nextBool())) ) {
                        max_gain   = target_gain;
                        max_gain_final = max_gain;
                        max_degree = delta_node_target->degree;
                        max_gainer = target_partition;
                }
        }
        // } endfor

        return (Gain)max_gain_final;
}






// Update (*config.delta) in terms of gain and degree (DON'T update them before calling this method)
inline void kway_graph_refinement_commons::move_delta_gain (PartitionConfig & config,
                                                            complete_boundary & boundary,
                                                            graph_access & G,
                                                            NodeID node,
                                                            const PartitionID par_old,
                                                            const PartitionID par_new) {

        ////////////////////// Update gains and degrees of 'node' //////////////////////
        // DEGREES in 'node''S NEIGHBORING BLOCKS DON'T CHANGE AFTER MOVEMENT. 
        // Check if 'par_new' queue already had 'node'
        std::vector<DELTA*> &delta_node = (*config.delta)[node].second;
        int size = delta_node.size();
        int old_part_pos = -1;
        int new_part_pos = -1;
        for (int i=0; i<size; i++) {
                if (delta_node[i]->block == par_old) {
                        old_part_pos = i;
                }
                if (delta_node[i]->block == par_new) {
                        new_part_pos = i;
                }
        }

        if (new_part_pos == -1) {
                DELTA* delta_node_new = new DELTA {
                        par_new,        // block
                        0,              // gain
                        0               // degree
                };
                delta_node.push_back(delta_node_new);
                new_part_pos = size++;
                // Below, notice that: gain_old_new == -gain_new_old
                delta_node_new->gain = (-1) * par_new_delta_gain(config, node, par_old, G); 
        } 
        int gain_par_new = delta_node[new_part_pos]->gain;    // cache gain value

        // CORRECTION OF GAINS IN 'node'
        std::vector<DELTA*>::iterator it;
        for (it = delta_node.begin(); it != delta_node.end(); ++it) {
                DELTA* visited_delta = *it;
                PartitionID p = visited_delta->block;
                visited_delta->gain -= gain_par_new;
        }

        // remove 'node' from 'par_old' queue if its degree is zero. 
        if (delta_node[old_part_pos]->degree == 0) {
                erase_block_delta(delta_node, old_part_pos);
        }




        std::vector< NodeID >* perm_rank = config.perm_rank;
        const NodeID perm_rank_par_old   = (*perm_rank)[par_old];
        const NodeID perm_rank_par_new   = (*perm_rank)[par_new];
        NodeID perm_rank_neighb;

        ////////////////////// Update gains and degrees of the neighboring nodes //////////////////////
        forall_out_edges(G, edge, node) {
                int w_edge = G.getEdgeWeight(edge);
                NodeID neighb = G.getEdgeTarget(edge);
                PartitionID par_neighb = G.getPartitionIndex(neighb);

                if ((*config.delta)[neighb].first != (*config.ref_layer) ) {   
                        // init_delta_gain (config, G, neighb);
                        continue;
                }

                std::vector<DELTA*> &delta_neighb = (*config.delta)[neighb].second;

                DELTA * delta_neighb_old = NULL;
                DELTA * delta_neighb_new = NULL;

                int size = delta_neighb.size();

                int old_part_pos = -1;
                int new_part_pos = -1;
                for (int i=0; i<size; i++) {
                        if (delta_neighb[i]->block == par_old) {
                                old_part_pos = i;
                                delta_neighb_old = delta_neighb[i];
                        }
                        if (delta_neighb[i]->block == par_new) {
                                new_part_pos = i;
                                delta_neighb_new = delta_neighb[i];
                        }
                }


                if ((*config.has_gains)[neighb]) {
                        int D_i_j_old = 0;
                        int D_i_j_new = 0;

                        matrix* D = config.D;
                        if (par_neighb == par_old) {
                                perm_rank_neighb = perm_rank_par_old;
                                D_i_j_new = D->get_xy(perm_rank_neighb, perm_rank_par_new) * (par_neighb != par_new);
                        } else if (par_neighb == par_new) {
                                perm_rank_neighb = perm_rank_par_new;
                                D_i_j_old = D->get_xy(perm_rank_neighb, perm_rank_par_old) * (par_neighb != par_old);
                        } else {
                                perm_rank_neighb = (*perm_rank)[par_neighb];
                                D_i_j_old = D->get_xy(perm_rank_neighb, perm_rank_par_old) * (par_neighb != par_old);
                                D_i_j_new = D->get_xy(perm_rank_neighb, perm_rank_par_new) * (par_neighb != par_new);
                        }

                        // CORRECTION OF GAINS IN 'neighb'
                        int w_j = w_edge;
                        int Dw_i_j_j_old = D_i_j_old * w_j;
                        int Dw_i_j_j_new = D_i_j_new * w_j;

                        (*config.has_gains)[neighb] = true;

                        if (par_old != par_neighb) {
                                delta_neighb_old->gain -= Dw_i_j_j_old;
                        }
                        if (new_part_pos != -1 && par_new != par_neighb) {
                                delta_neighb_new->gain += (Dw_i_j_j_new);
                        }

                        // Contribution gain in case 'neighb' goes from 'par_neighb' to 'target_part_k'
                        std::vector<DELTA*>::iterator it_inn;
                        for (it_inn = delta_neighb.begin(); it_inn != delta_neighb.end(); ++it_inn)
                        {
                                DELTA* delta_neighb_k = *it_inn;
                                PartitionID target_part_k = delta_neighb_k->block; 
                                if (target_part_k == par_neighb)       
                                        continue;

                                int D_j_k_old = 0;
                                int D_j_k_new = 0;
                                NodeID perm_rank_target_k;
                                if (target_part_k == par_old) {
                                        perm_rank_target_k = perm_rank_par_old;
                                        D_j_k_new = D->get_xy(perm_rank_par_new, perm_rank_target_k);
                                } else if (target_part_k == par_new) {
                                        perm_rank_target_k = perm_rank_par_new;
                                        D_j_k_old = D->get_xy(perm_rank_par_old, perm_rank_target_k);
                                } else {
                                        perm_rank_target_k = (*perm_rank)[target_part_k];
                                        D_j_k_old = D->get_xy(perm_rank_par_old, perm_rank_target_k);
                                        D_j_k_new = D->get_xy(perm_rank_par_new, perm_rank_target_k);
                                }

                                int tmp_eq = (par_old == par_neighb);
                                delta_neighb_k->gain  += tmp_eq * (D_j_k_old*w_j)
                                                        -(!tmp_eq)*(target_part_k != par_old) * (Dw_i_j_j_old - D_j_k_old*w_j);

                                tmp_eq = (par_new == par_neighb);
                                delta_neighb_k->gain  += - tmp_eq *(D_j_k_new*w_j)
                                                        +(!tmp_eq)*(target_part_k != par_new)*(Dw_i_j_j_new - D_j_k_new*w_j);
                        }
                }

                // CORRECTION OF DEGREES IN 'neighb'
                delta_neighb_old->degree -= w_edge; 
                // remove neighb from par_old queue if resulting degree is zero. 
                if (delta_neighb_old->degree == 0 && par_old != par_neighb) {
                        erase_block_delta(delta_neighb, old_part_pos);
                }
                // if neighb is not in par_new queue, add it
                if (new_part_pos == -1) {
                        DELTA* delta_neighb_new = new DELTA {
                                par_new,        // block
                                0,              // gain
                                w_edge          // degree
                        };
                        delta_neighb.push_back(delta_neighb_new);
                        if ((*config.has_gains)[neighb]) {
                                delta_neighb_new->gain = par_new_delta_gain(config, neighb, par_new, G);
                        }
                } else {
                        delta_neighb_new->degree += w_edge;
                }
        } endfor
}



// Calculate (from scratch) the gain of 'node' if it moves to 'par_new'.
// PS: This method assumes that the degrees of all neighboring partitions 
// are correct, including 'par_new' and also 'G.getPartitionIndex(node)'
inline int kway_graph_refinement_commons::par_new_delta_gain ( PartitionConfig & config,
                                                               NodeID node,
                                                               PartitionID par_new,
                                                               graph_access & G) {
        
        PartitionID source_partition = G.getPartitionIndex(node);
        std::vector< NodeID >* perm_rank = config.perm_rank;
        matrix* D = config.D;
        NodeID perm_rank_source = (*perm_rank)[source_partition];

        std::vector<DELTA*> &delta_node = (*config.delta)[node].second;
        DELTA* delta_node_j = NULL;
        DELTA* delta_node_source = NULL;

        PartitionID target_part_j = par_new;

        if (target_part_j == source_partition) {
                std::cout << "UNEXPECTED BEHAVIOR 2 - notify Marcelo\n";
                return 0;
        }
        int new_gain = 0;
        
        NodeID perm_rank_target_j   = (*perm_rank)[target_part_j];

        // contribute to gain in case 'node' goes from 'source_partition' to 'target_part_k'
        std::vector<DELTA*>::iterator it_inn;
        for (it_inn = delta_node.begin(); it_inn != delta_node.end(); ++it_inn)
        {
                DELTA* visited_delta = (*it_inn);
                PartitionID target_part_k = visited_delta->block; 
                if (target_part_k == target_part_j || target_part_k == source_partition) {
                        if (target_part_k == target_part_j) {
                                delta_node_j = visited_delta;
                        }
                        if (target_part_k == source_partition) {
                                delta_node_source = visited_delta;
                        }
                        continue;
                }
                NodeID perm_rank_target_k        = (*perm_rank)[target_part_k];

                int D_i_k = D->get_xy(perm_rank_source, perm_rank_target_k);
                int w_k = visited_delta->degree;
                int D_j_k = D_i_k;
                if (perm_rank_source != perm_rank_target_j) {
                        D_j_k = D->get_xy(perm_rank_target_j, perm_rank_target_k);
                }
                new_gain  += w_k * (D_i_k - D_j_k);
        }

        int w_i = delta_node_source->degree;
        int D_i_j = D->get_xy(perm_rank_source, perm_rank_target_j);
        int w_j = delta_node_j->degree;
        int Dw_i_j_j = D_i_j * w_j;
        
        new_gain += Dw_i_j_j - D_i_j*w_i;

        return new_gain;
}






// for all incident partitions in 'node', compute degrees from scratch
inline int kway_graph_refinement_commons::init_delta_degrees (PartitionConfig & config,
                                                              graph_access & G, 
                                                              NodeID & node,
                                                              const PartitionID lhs,
                                                              const PartitionID rhs,
                                                              EdgeWeight & int_obj_func,
                                                              EdgeWeight & ext_obj_func,
                                                              bool & has_ext_edge) {
        
        std::vector< NodeID >* perm_rank = config.perm_rank;
        matrix* D = config.D;
        NodeID perm_rank_lhs = (*perm_rank)[lhs];
        NodeID perm_rank_rhs = (*perm_rank)[rhs];
        int D_lhs_rhs = D->get_xy(perm_rank_lhs,perm_rank_rhs);
        
        (*config.delta)[node].first = (*config.ref_layer);
        (*config.has_gains)[node] = false;

        std::vector<DELTA*> &delta_node = (*config.delta)[node].second;
        // to make it consistent
        for (DELTA* obj : delta_node) {
                delete obj;
        }
        delta_node.clear();

        std::unordered_map<PartitionID,DELTA*> stored;

        DELTA* delta_node_lhs = new DELTA {
                lhs,    // block
                0,      // gain
                0       // degree
        };
        delta_node.push_back(delta_node_lhs);
        stored[lhs] = delta_node_lhs;

        DELTA* delta_node_rhs = new DELTA {
                rhs,    // block
                0,      // gain
                0       // degree
        };
        delta_node.push_back(delta_node_rhs);
        stored[rhs] = delta_node_rhs;


        int & delta_node_lhs_degree = delta_node_lhs->degree;
        int & delta_node_rhs_degree = delta_node_rhs->degree;
        delta_node_lhs_degree = 0;
        delta_node_rhs_degree = 0;

        int curr_OF_lhs = 0;
        int curr_OF_rhs = 0;
        bool update_is_difficult = false;
        has_ext_edge = false;

        // config.delta[node] stores sum(weight(edge from 'node')) incident in each partition  
        forall_out_edges(G, e, node) {
                NodeID target                = G.getEdgeTarget(e);
                PartitionID target_partition = G.getPartitionIndex(target);
                EdgeWeight e_weight          = G.getEdgeWeight(e);

                if(target_partition == rhs) {
                        has_ext_edge = true;
                        curr_OF_lhs += e_weight * D_lhs_rhs;
                        delta_node_rhs_degree += e_weight;
                } else if (target_partition == lhs) {
                        delta_node_lhs_degree += e_weight;
                        curr_OF_rhs += e_weight * D_lhs_rhs;
                } else {
                        NodeID perm_rank_target = (*perm_rank)[target_partition];
                        update_is_difficult = true;
                        curr_OF_lhs += e_weight * D->get_xy(perm_rank_lhs, perm_rank_target);
                        curr_OF_rhs += e_weight * D->get_xy(perm_rank_rhs, perm_rank_target);
                        if (stored.find(target_partition) == stored.end()) {
                                DELTA * delta_node_target = new DELTA {
                                        target_partition,       // block
                                        0,                      // gain
                                        e_weight                // degree
                                };      
                                stored[target_partition] = delta_node_target;
                                delta_node.push_back(delta_node_target);
                        } else {
                                stored[target_partition]->degree += e_weight;
                        }
                }
        } endfor
        int_obj_func = curr_OF_lhs;
        ext_obj_func = curr_OF_rhs;
        return update_is_difficult;
}




inline void kway_graph_refinement_commons::move_delta_degrees (PartitionConfig & config,
                                                               complete_boundary & boundary,
                                                               graph_access & G,
                                                               NodeID node,
                                                               const PartitionID par_old,
                                                               const PartitionID par_new) {
        (*config.has_gains)[node] = false;

        ////////////////////// Update degrees of the neighboring nodes //////////////////////
        forall_out_edges(G, edge, node) {
                NodeID neighb = G.getEdgeTarget(edge);

                if ((*config.delta)[neighb].first != (*config.ref_layer)) {
                        // init_delta_gain (config, G, neighb);
                        continue;
                }

                PartitionID par_neighb = G.getPartitionIndex(neighb);
                std::vector<DELTA*> &delta_neighb = (*config.delta)[neighb].second;
                int w_edge = G.getEdgeWeight(edge);

                DELTA * delta_neighb_old = NULL;
                DELTA * delta_neighb_new = NULL;
                int size = delta_neighb.size();
                int old_part_pos = -1;
                int new_part_pos = -1;
                for (int i=0; i<size; i++) {
                        if (delta_neighb[i]->block == par_old) {
                                old_part_pos = i;
                                delta_neighb_old = delta_neighb[i];
                        }
                        if (delta_neighb[i]->block == par_new) {
                                new_part_pos = i;
                                delta_neighb_new = delta_neighb[i];
                        }
                }

                int & delta_neighb_old_degree = delta_neighb_old->degree;
                (*config.has_gains)[neighb] = false;

                // CORRECTION OF DEGREES IN 'neighb'
                delta_neighb_old_degree -= w_edge; 
                // remove neighb from par_old queue if resulting degree is zero. 
                if (delta_neighb_old_degree == 0) {
                        erase_block_delta(delta_neighb, old_part_pos);
                }
                // if neighb is not in par_new queue, add it
                if (new_part_pos == -1) {
                        delta_neighb_new = new DELTA {
                                par_new,        // block
                                0,              // gain
                                w_edge          // degree
                        };
                        delta_neighb.push_back(delta_neighb_new);
                        } else {
                        delta_neighb_new->degree += w_edge;
                }
        } endfor
        
       
}






// Compute delta gain for a node with degrees already computed correctly
inline int kway_graph_refinement_commons::init_gain_from_degree (PartitionConfig & config,
                                                                 complete_boundary & boundary,
                                                                 graph_access & G, 
                                                                 NodeID & node, 
                                                                 PartitionID & max_gainer, 
                                                                 EdgeWeight & max_degree,
                                                                 bool ignore_balance = true,
                                                                 bool label_prop = false) {
        
        PartitionID source_partition = G.getPartitionIndex(node);
        std::vector< NodeID >* perm_rank = config.perm_rank;
        matrix* D = config.D;
        NodeID perm_rank_source = (*perm_rank)[source_partition];

        if ((*config.delta)[node].first != (*config.ref_layer)) {
                std::cout << "ERROR 20 (Marcelo): It should not enter this method";
        }
        (*config.delta)[node].first = (*config.ref_layer);
        
        int max_gain = (!label_prop || boundary.getBlockWeight(source_partition) > config.upper_bound_partition) 
                       * std::numeric_limits<int>::min();
        int max_gain_final = 0;
        
        max_degree   = 0;
        max_gainer   = INVALID_PARTITION;

        std::vector<DELTA*> &delta_node = (*config.delta)[node].second;
        NodeWeight this_nodes_weight = G.getNodeWeight(node);
        DELTA * delta_node_source = NULL;

        int curr_OF = 0;
        std::vector<DELTA*>::iterator it;
        for (it = delta_node.begin(); it != delta_node.end(); ++it) {
                DELTA* delta_node_target = (*it);
                PartitionID target_partition = delta_node_target->block;
                NodeID perm_rank_target   = (*perm_rank)[target_partition];

                if (target_partition == source_partition) {
                        delta_node_source = delta_node_target;
                        continue;
                }

                curr_OF += delta_node_target->degree * D->get_xy(perm_rank_source, perm_rank_target);
                int target_gain = 0;
                std::vector<DELTA*>::iterator it_inner;
                for (it_inner = delta_node.begin(); it_inner != delta_node.end(); ++it_inner)
                {
                        DELTA* delta_node_target_inner = (*it_inner);
                        PartitionID target_partition_inner = delta_node_target_inner->block; 
                        if (target_partition_inner == target_partition)       
                                continue;
                        NodeID target_value_inner = delta_node_target_inner->degree;
                        NodeID perm_rank_target_inner = (*perm_rank)[target_partition_inner];
                        target_gain -= target_value_inner * D->get_xy(perm_rank_target, perm_rank_target_inner);
                }
                delta_node_target->gain = target_gain;
        }

        if (delta_node_source == NULL) {
                delta_node_source = new DELTA {
                        source_partition,       // block
                        -curr_OF,               // gain
                        0                       // degree
                };
                delta_node.push_back(delta_node_source);
        } else {
                delta_node_source->gain = -curr_OF;
        }


        for (it = delta_node.begin(); it != delta_node.end(); ++it) {
                DELTA* delta_node_target = (*it);
                PartitionID target_partition = delta_node_target->block;
                delta_node_target->gain += curr_OF;
                Gain target_gain = delta_node_target->gain;
                if( target_partition!=source_partition && (
                   (ignore_balance||boundary.getBlockWeight(target_partition)+this_nodes_weight<=config.upper_bound_partition) 
                    && ( target_gain > max_gain || (target_gain == max_gain && random_functions::nextBool()))) ) {
                        max_gain       = target_gain;
                        max_gain_final = max_gain; 
                        max_degree     = delta_node_target->degree;
                        max_gainer     = target_partition;
                }
        }

        (*config.has_gains)[node] = true;

        return (Gain)max_gain_final;
}


// Erases in O(1) time a delta block from some node's vector of neighboring delta blocks
inline void kway_graph_refinement_commons::erase_block_delta(std::vector<DELTA*> &delta_node, size_t position) {
        size_t len = delta_node.size();
        if (len == 0) {
                std::cout << "Unexpected Error Marcelo 77\n"; 
                return;
        }
        delete delta_node[position];
        if (position != len-1) {
                delta_node[position] = delta_node[len-1];
        } 
        delta_node.pop_back();
}




inline void kway_graph_refinement_commons::ignore_delta_gain(PartitionConfig & config,
                                                             graph_access & G, 
                                                             NodeID & node, 
                                                             int increment) {
        std::vector<std::pair<int,std::vector<DELTA*>>> &gen_delta = (*config.delta);
        gen_delta[node].first += increment;
        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                gen_delta[target].first += increment;
        } endfor
}






#endif /* end of include guard: KWAY_GRAPH_REFINEMENT_COMMONS_PVGY97EW */

