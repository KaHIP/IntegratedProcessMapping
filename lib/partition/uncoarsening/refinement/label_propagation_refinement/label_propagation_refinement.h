/******************************************************************************
 * label_propagation_refinement.h  
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/


#ifndef LABEL_PROPAGATION_REFINEMENT_R4XW141Y
#define LABEL_PROPAGATION_REFINEMENT_R4XW141Y

#include "definitions.h"
#include "../refinement.h"
#include "data_structure/matrix/matrix.h"
#include "partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_commons.h"

class label_propagation_refinement : public refinement {
public:
        label_propagation_refinement();
        virtual ~label_propagation_refinement();

        virtual EdgeWeight perform_refinement(PartitionConfig & config, 
                                              graph_access & G, 
                                              complete_boundary & boundary); 

        EdgeWeight perform_refinement_mapping(PartitionConfig & config, graph_access & G, matrix & D, 
                                              std::vector< NodeID > & perm_rank, complete_boundary & boundary); 
        
        bool is_boundary(NodeID node, graph_access & G);
private:
        kway_graph_refinement_commons * commons;
};


#endif /* end of include guard: LABEL_PROPAGATION_REFINEMENT_R4XW141Y */
