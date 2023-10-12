/******************************************************************************
 * two_way_fm.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef TWOWAY_FM_YLYN82Y1
#define TWOWAY_FM_YLYN82Y1

#include <vector>

#include "data_structure/graph_access.h"
#include "data_structure/priority_queues/priority_queue_interface.h"
#include "definitions.h"
#include "partition_config.h"
#include "uncoarsening/refinement/quotient_graph_refinement/complete_boundary.h"
#include "uncoarsening/refinement/quotient_graph_refinement/partial_boundary.h"
#include "uncoarsening/refinement/quotient_graph_refinement/two_way_refinement.h"
#include "vertex_moved_hashtable.h"
#include "data_structure/matrix/matrix.h"
#include "uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_commons.h"


class two_way_fm : public two_way_refinement {
        public:
                two_way_fm( );
                virtual ~two_way_fm();
                EdgeWeight perform_refinement(PartitionConfig & config, 
                                graph_access & G, 
                                complete_boundary & boundary, 
                                std::vector<NodeID> & lhs_pq_start_nodes, 
                                std::vector<NodeID> & rhs_pq_start_nodes,
                                boundary_pair * refinement_pair,        
                                NodeWeight & lhs_part_weight,
                                NodeWeight & rhs_part_weight,
                                EdgeWeight & cut,
                                bool & something_changed);

                inline bool int_ext_degree(graph_access & G, 
                                const NodeID & node,
                                const PartitionID lhs,
                                const PartitionID rhs,
                                EdgeWeight & int_degree,
                                EdgeWeight & ext_degree);

                inline bool int_ext_degree_map( PartitionConfig & config,
                                                graph_access & G, 
                                                const NodeID & node,
                                                const PartitionID lhs,
                                                const PartitionID rhs,
                                                EdgeWeight & int_obj_func,
                                                EdgeWeight & ext_obj_func,
                                                bool & has_ext_edge);
                
                inline bool int_ext_degree_delta( PartitionConfig & config,
                                                  graph_access & G, 
                                                  NodeID & node,
                                                  const PartitionID lhs,
                                                  const PartitionID rhs,
                                                  EdgeWeight & int_obj_func,
                                                  EdgeWeight & ext_obj_func,
                                                  bool & has_ext_edge);


        private:
                void init_queue_with_boundary(PartitionConfig & config,
                                              graph_access & G,
                                              std::vector<NodeID> &bnd_nodes,
                                              refinement_pq * queue,                     
                                              PartitionID partition_of_boundary, 
                                              PartitionID other); 


                void move_node(PartitionConfig & config, 
                               graph_access & G,
                               NodeID & node,
                               vertex_moved_hashtable & moved_idx,
                               refinement_pq * from_queue,
                               refinement_pq * to_queue,
                               PartitionID from,
                               PartitionID to,
                               boundary_pair * pair,
                               NodeWeight * from_part_weight,
                               NodeWeight * to_part_weight,
                               complete_boundary & boundary);

                void move_node_back(PartitionConfig & config, 
                                    graph_access & G,
                                    NodeID & node,
                                    vertex_moved_hashtable & moved_idx,
                                    refinement_pq * from_queue,
                                    refinement_pq * to_queue,
                                    PartitionID from, 
                                    PartitionID to,
                                    boundary_pair * pair,
                                    NodeWeight * from_part_weight,
                                    NodeWeight * to_part_weight,
                                    complete_boundary & boundary); 


                kway_graph_refinement_commons * commons=NULL;

                ///////////////////////////////////////////////////////////////////////////
                //Assertions
                ///////////////////////////////////////////////////////////////////////////
#ifndef NDEBUG
                //assert that every node in the lhs boundary has external degree > 0
                bool assert_only_boundary_nodes(graph_access & G, 
                                                PartialBoundary & lhs_boundary, 
                                                PartitionID lhs, 
                                                PartitionID rhs);

                //assert that every node with ext degree > 0 is lhs boundary 
                bool assert_every_boundary_nodes(graph_access & G, 
                                                 PartialBoundary & lhs_boundary, 
                                                 PartitionID lhs, 
                                                 PartitionID rhs);

                //check all of the possible compinations of the two assertions above
                bool assert_directed_boundary_condition(graph_access & G, 
                                                        complete_boundary & boundary, 
                                                        PartitionID lhs, 
                                                        PartitionID rhs);
#endif

};

inline bool two_way_fm::int_ext_degree( graph_access & G, 
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
                NodeID target = G.getEdgeTarget(e);
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



inline bool two_way_fm::int_ext_degree_map( PartitionConfig & config,
                                            graph_access & G, 
                                            const NodeID & node,
                                            const PartitionID lhs,
                                            const PartitionID rhs,
                                            EdgeWeight & int_obj_func,
                                            EdgeWeight & ext_obj_func,
                                            bool & has_ext_edge) {


        ASSERT_EQ(lhs, G.getPartitionIndex(node));

        matrix* D = config.D;
        std::vector< NodeID >* perm_rank = config.perm_rank;
        NodeID perm_rank_lhs  = (*perm_rank)[lhs];
        NodeID perm_rank_rhs  = (*perm_rank)[rhs];
        int D_lhs_rhs = D->get_xy(perm_rank_lhs,perm_rank_rhs);

        int_obj_func               = 0;
        ext_obj_func               = 0;
        bool update_is_difficult = false;

        has_ext_edge = false;

        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                PartitionID targets_partition = G.getPartitionIndex(target);
                NodeID perm_rank_target  = (*perm_rank)[targets_partition];
                EdgeWeight e_weight = G.getEdgeWeight(e);

                if(targets_partition == rhs) {
                        has_ext_edge = true;
                        int_obj_func += e_weight * D_lhs_rhs;
                } else if (targets_partition == lhs) {
                        ext_obj_func += e_weight * D_lhs_rhs;
                } else {
                        NodeID perm_rank_target = (*config.perm_rank)[targets_partition];
                        update_is_difficult = true;
                        int_obj_func += e_weight * D->get_xy(perm_rank_lhs, perm_rank_target);
                        ext_obj_func += e_weight * D->get_xy(perm_rank_rhs, perm_rank_target);
                }
        } endfor

        return update_is_difficult;
}



inline bool two_way_fm::int_ext_degree_delta( PartitionConfig & config,
                                              graph_access & G, 
                                              NodeID & node,
                                              const PartitionID lhs,
                                              const PartitionID rhs,
                                              EdgeWeight & int_obj_func,
                                              EdgeWeight & ext_obj_func,
                                              bool & has_ext_edge) {

        ASSERT_EQ(lhs, G.getPartitionIndex(node));
        if (commons == NULL) {
                commons = new kway_graph_refinement_commons(config);
        }

        if ((*config.delta)[node].first != (*config.ref_layer)) {
                return commons->init_delta_degrees(config, G, node, lhs, rhs, int_obj_func, ext_obj_func, has_ext_edge);
        }

        has_ext_edge = false;
        bool update_is_difficult = false;
        matrix* D = config.D;
        std::vector< NodeID >* perm_rank = config.perm_rank;
        NodeID perm_rank_lhs  = (*perm_rank)[lhs];
        NodeID perm_rank_rhs  = (*perm_rank)[rhs];
        int D_lhs_rhs = D->get_xy(perm_rank_lhs,perm_rank_rhs);
        int_obj_func = 0;
        ext_obj_func = 0;

        std::vector<DELTA*> &delta_node = (*config.delta)[node].second;
        std::vector<DELTA*>::iterator it, end;
        end = delta_node.end();
        for (it = delta_node.begin(); it != end; ++it)
        {
                DELTA* visited_delta = *it;
                PartitionID target_part = visited_delta->block; 
                EdgeWeight e_weight = visited_delta->degree;

                if (e_weight <= 0) 
                        continue;
                
                if(target_part == rhs) {
                        has_ext_edge = true;
                        int_obj_func += e_weight * D_lhs_rhs;
                } else if (target_part == lhs) {
                        ext_obj_func += e_weight * D_lhs_rhs;
                } else {
                        NodeID perm_rank_target = (*config.perm_rank)[target_part];
                        update_is_difficult = true;
                        int_obj_func += e_weight * D->get_xy(perm_rank_lhs, perm_rank_target);
                        ext_obj_func += e_weight * D->get_xy(perm_rank_rhs, perm_rank_target);
                }
        }
        return update_is_difficult;
}


#endif /* end of include guard: TWO_WAY_FM_YLYN82Y1 */
