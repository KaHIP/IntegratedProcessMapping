/******************************************************************************
 * mixed_refinement.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "cycle_improvements/cycle_refinement.h"
#include "kway_graph_refinement/kway_graph_refinement.h"
#include "kway_graph_refinement/multitry_kway_fm.h"
#include "mixed_refinement.h"
#include "quotient_graph_refinement/quotient_graph_refinement.h"
#include "label_propagation_refinement/label_propagation_refinement.h"

mixed_refinement::mixed_refinement() {

}

mixed_refinement::~mixed_refinement() {

}

EdgeWeight mixed_refinement::perform_refinement(PartitionConfig & config, graph_access & G, complete_boundary & boundary) {
        refinement* refine                              = new quotient_graph_refinement();
        refinement* kway                                = new kway_graph_refinement();
        multitry_kway_fm* multitry_kway                 = new multitry_kway_fm();
        cycle_refinement* cycle_refine                  = new cycle_refinement();
        label_propagation_refinement* label_refine      = new label_propagation_refinement();

        EdgeWeight overall_improvement = 0; 
        //call refinement
        if(!config.integrated_mapping || (config.bipartition_gp_local_search && 
                                        config.initial_partitioning_type != INITIAL_PARTITIONING_RECPARTITION) ) {
                if(config.no_change_convergence) {
                        bool sth_changed = true;
                        while(sth_changed) {
                                EdgeWeight improvement = 0;
                                if(config.corner_refinement_enabled) {
                                        improvement += kway->perform_refinement(config, G, boundary);
                                }

                                if(!config.quotient_graph_refinement_disabled) {
                                        improvement += refine->perform_refinement(config, G, boundary);
                                }

                                overall_improvement += improvement;
                                sth_changed = improvement != 0;
                        }

                } else {
                        if(config.corner_refinement_enabled) {
                                overall_improvement += kway->perform_refinement(config, G, boundary);
                        } 

                        if(!config.quotient_graph_refinement_disabled) {
                                overall_improvement += refine->perform_refinement(config, G, boundary);
                        }

                        if(config.kaffpa_perfectly_balanced_refinement) {
                                overall_improvement += cycle_refine->perform_refinement(config, G, boundary);
                        }
                }
        } 


        


        if(config.qap_0quotient_ref && config.initial_partitioning_type == INITIAL_PARTITIONING_RECPARTITION) {
//                 quality_metrics qm;
//                 int qap;
//                 qap = qm.total_qap_noquotient(G, *config.D, *config.perm_rank );
//                 std::cout <<  "BEFORE-qap_0quotient_ref - J(C,D,Pi') = " << qap << std::endl;
// std::cout << "BEF - Balance \t" << qm.balance(G) << std::endl;
                overall_improvement += refine->perform_refinement(config, G, boundary);

//                 // std::cout << "AFT - Balance \t" << qm.balance(G) << std::endl;
//                 // qap = qm.total_qap_noquotient(G, *config.D, *config.perm_rank );
//                 // std::cout <<  "AFTER-qap_0quotient_ref  - J(C,D,Pi') = " << qap << std::endl;
        }




        
        if ( config.qap_blabel_propagation_refinement && config.initial_partitioning_type == INITIAL_PARTITIONING_RECPARTITION ) {
                overall_improvement += label_refine->perform_refinement_mapping(config, G, *config.D, *config.perm_rank, boundary);
        }




        if(config.qap_bquotient_ref && !config.qap_0quotient_ref && config.initial_partitioning_type == INITIAL_PARTITIONING_RECPARTITION) {
                overall_improvement += refine->perform_refinement(config, G, boundary);
        }





        if ( config.qap_bkway_fm && config.initial_partitioning_type == INITIAL_PARTITIONING_RECPARTITION ) {
                // quality_metrics qm;
                // int qap;
                // qap = qm.total_qap_noquotient(G, *config.D, *config.perm_rank );
                // std::cout <<  "BEFORE-qap_bkway_fm - J(C,D,Pi') = " << qap << std::endl;
                // std::cout << "BEF - Balance \t" << qm.balance(G) << std::endl;

                overall_improvement += kway->perform_refinement(config, G, boundary);

                // std::cout << "AFT - Balance \t" << qm.balance(G) << std::endl;
                // qap = qm.total_qap_noquotient(G, *config.D, *config.perm_rank );
                // std::cout <<  "AFTER-qap_bkway_fm  - J(C,D,Pi') = " << qap << std::endl;
        }

        if ( config.qap_bmultitry_kway_fm && config.initial_partitioning_type == INITIAL_PARTITIONING_RECPARTITION ) {
                overall_improvement += multitry_kway->perform_refinement(config, G, boundary, config.local_multitry_rounds, 
                                                                                true, config.local_multitry_fm_alpha); 
        }


        
        
        
        if(config.qap_quotient_ref && config.initial_partitioning_type == INITIAL_PARTITIONING_RECPARTITION) {
                overall_improvement += refine->perform_refinement(config, G, boundary);
        }




        
        if ( config.qap_label_propagation_refinement && config.initial_partitioning_type == INITIAL_PARTITIONING_RECPARTITION ) {
                overall_improvement += label_refine->perform_refinement_mapping(config, G, *config.D, *config.perm_rank, boundary);
                // quality_metrics qm;
                // int qap;
                // qap = qm.total_qap_noquotient(G, *config.D, *config.perm_rank );
                // std::cout <<  "AFTER-qap_label_propagation_refinement - J(C,D,Pi') = " << qap << std::endl;
                // std::cout << "BEF - Balance \t" << qm.balance(G) << std::endl;
        }


        if ( config.qap_multitry_kway_fm && config.initial_partitioning_type == INITIAL_PARTITIONING_RECPARTITION ) {
//                 quality_metrics qm;
//                 int qap;
//                 qap = qm.total_qap_noquotient(G, *config.D, *config.perm_rank );
//                 std::cout <<  "BEFORE-qap_multitry_kway_fm - J(C,D,Pi') = " << qap << std::endl;
// std::cout << "BEF - Balance \t" << qm.balance(G) << std::endl;
                
                overall_improvement += multitry_kway->perform_refinement(config, G, boundary, config.local_multitry_rounds, 
                                                                                true, config.local_multitry_fm_alpha); 

//  std::cout << "AFT - Balance \t" << qm.balance(G) << std::endl;
//                 qap = qm.total_qap_noquotient(G, *config.D, *config.perm_rank );
//                 std::cout <<  "AFTER-qap_multitry_kway_fm  - J(C,D,Pi') = " << qap << std::endl;                                                                               
        }


        if ( config.qap_kway_fm && config.initial_partitioning_type == INITIAL_PARTITIONING_RECPARTITION ) {
//                 quality_metrics qm;
//                 int qap;
//                 qap = qm.total_qap_noquotient(G, *config.D, *config.perm_rank );
//                 std::cout <<  "BEFORE-qap_kway_fm - J(C,D,Pi') = " << qap << std::endl;
// std::cout << "BEF - Balance \t" << qm.balance(G) << std::endl;
        
                overall_improvement += kway->perform_refinement(config, G, boundary);

// std::cout << "AFT - Balance \t" << qm.balance(G) << std::endl;
//                 qap = qm.total_qap_noquotient(G, *config.D, *config.perm_rank );
//                 std::cout <<  "AFTER-qap_kway_fm  - J(C,D,Pi') = " << qap << std::endl;
        }


        if ( config.qap_alabel_propagation_refinement && config.initial_partitioning_type == INITIAL_PARTITIONING_RECPARTITION ) {
                overall_improvement += label_refine->perform_refinement_mapping(config, G, *config.D, *config.perm_rank, boundary);
        }






        if (config.initial_partitioning_type == INITIAL_PARTITIONING_RECPARTITION && config.use_delta_gains) {
                // config.delta->clear();
                (*config.ref_layer)++;
        }


        delete refine;
        delete kway;
        delete multitry_kway;
        delete cycle_refine;

        return overall_improvement;
}

