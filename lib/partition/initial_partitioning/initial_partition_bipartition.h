/******************************************************************************
 * initial_partition_bipartition.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef INITIAL_PARTITION_BIPARTITION_HMA7329W
#define INITIAL_PARTITION_BIPARTITION_HMA7329W

#include "initial_partitioner.h"

class initial_partition_bipartition : public initial_partitioner {
public:
        initial_partition_bipartition();
        virtual ~initial_partition_bipartition();

        void initial_partition( PartitionConfig & config, const unsigned int seed,  
                                graph_access & G, int* partition_map, int ismultisec=0); 

        void initial_partition( PartitionConfig & config, const unsigned int seed,  
                                graph_access & G, 
                                int* xadj,
                                int* adjncy, 
                                int* vwgt, 
                                int* adjwgt,
                                int* partition_map, 
                                int ismultisec=0); 

};


#endif /* end of include guard: INITIAL_PARTITION_BIPARTITION_HMA7329W */
