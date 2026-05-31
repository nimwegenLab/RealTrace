#include <Eigen/Cholesky>
#include <chrono>

#include "correlation_tree.h"

/* --------------------------------------------------------------------------
* Sample trajectory
* -------------------------------------------------------------------------- */

void calc_first_joint_distribution(const std::vector<std::vector<double>> &params_vecs, 
                                    MOMAdata &cell, 
                                    int n){

    /* Calculate P(z_n+1, z_n | D_n) */

    Eigen::MatrixXd D(2,2);
    Gaussian combined_joint;
    Affine_gaussian conditional;
    int m=1;

    /* Same as joint calculation for m=1 */
    std::vector<double> params_vec = params_vecs[cell.segment[n+m]];

    /* ------------ Posterior ------------ */ 
    /* include x and g -> P(z_n+1, z_n | D_n1+1) */        
    if (cell.noise_model == "scaled"){
        D <<  params_vec[7], 0, 0,  abs(params_vec[8]*(cell.mean(1)+cell.fp_auto));
    }
    else {
        D <<  params_vec[7], 0, 0,  params_vec[8];
    }

    cell.joint = include_measurement(cell.joint, D, cell.log_length(n+m), cell.fp(n+m)); 

    combined_joint = incorporate_backward_prob(seperate_gaussian(cell.joint ), 
                                                cell.mean_backward[n+m], 
                                                cell.cov_backward[n+m], 
                                                params_vec);

    cell.consecutive_joints[n+m] = combined_joint.flip_xy(); // flip P(z_n+1, z_n)
    cell.conditionals[n+m] = seperate_gaussian(cell.consecutive_joints[n+m]).conditional;

}



/* ======================================================== */
/* Looping time points and cells */
/* ======================================================== */

void first_joint_distributions_recr(  const std::vector<std::vector<double>> &params_vecs, 
                                MOMAdata *cell, 
                                int n, 
                                bool is_joint_at_division){
    /*  
    * Recursive implementation that applies the function calc_joint_distributions to every cell in the genealogy
    */
    if (cell == nullptr)
        return;
    if (!is_joint_at_division){
        calc_first_joint_distribution(params_vecs, *cell, n);
        return;
    }
    first_joint_distributions_recr(params_vecs, cell->daughter1, -1, false);
    first_joint_distributions_recr(params_vecs, cell->daughter2, -1, false);
}


void sc_first_joint_distributions(const std::vector<std::vector<double>> &params_vecs, 
                            MOMAdata &cell){
    /* 
    * Calculates first joints for all time points in this cell. 
    */

    for (size_t n=0; n<cell.time.size(); ++n){

        // P(z_n+1, z_n | D_n) (using Theta_n)
        if (n<cell.time.size()-1){
            cell.joint = consecutive_joint(params_vecs[cell.segment[n]], cell, n);
            first_joint_distributions_recr(params_vecs, &cell, n, false);
        }
        else{
            if (cell.daughter1 != nullptr){
                cell.joint = consecutive_joint_cell_division(params_vecs[cell.segment[n]], cell, n);
                cell.daughter1->joint = cell.joint;
            }
            if (cell.daughter2 != nullptr){
                cell.daughter2->joint = cell.joint;
            }
            first_joint_distributions_recr(params_vecs, &cell, n, true);
        }
    }
}

void first_joint_distributions(const std::vector<std::vector<double>> &params_vecs, 
                                std::vector<MOMAdata> &cells){
    /* Go over all cells to call sc_first_joint_distributions for each cell */
    for (size_t i=0; i<cells.size(); ++i){
        sc_first_joint_distributions(params_vecs, cells[i]);
    }
}


void init_cell_variables(std::vector<MOMAdata> &cells){
    /* initialize consecutive_joints vecor of each cell to have the same length as time points*/
    Gaussian dummy_g;
    Affine_gaussian dummy_c;
    Eigen::VectorXd dummy_m = Eigen::VectorXd::Zero(4);

    for (size_t i=0; i<cells.size(); ++i){
        for (size_t t=0; t<cells[i].time.size(); ++t){
            cells[i].consecutive_joints.push_back(dummy_g);
            cells[i].conditionals.push_back(dummy_c);
            cells[i].sample.push_back(dummy_m);
        }
    }
}



/* --------------------------------------------------------------------------
* Sampling
* -------------------------------------------------------------------------- */

Eigen::VectorXd sample_gaussian(Eigen::VectorXd mean, Eigen::MatrixXd cov){
    /* Sample a 4 dim vector from 4 dim gaussian with mean and cov (non-deterministic, ie system time as seed) */

    // Cholesky decomposition
    Eigen::MatrixXd L = cov.llt().matrixL();

    // Generate samples
    unsigned seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    
    // std::mt19937 
    std::default_random_engine generator(seed);
    std::normal_distribution<double> distribution(0.0, 1.0);
    
    Eigen::VectorXd z(4);
    z << distribution(generator), distribution(generator), distribution(generator), distribution(generator);

    return mean + L * z; 
}

/* =====  ===== */

void sc_sample_trajectory(MOMAdata &cell){
    for (size_t t=0; t<cell.time.size(); ++t){
            if(cell.is_root() && t==0){
            // sample from one time point posterior 
            cell.sample[t] = sample_gaussian(cell.mean_prediction[t], cell.cov_prediction[t]);
        }
        else{
            /* get last sample, ie for zn-1*/
            Eigen::VectorXd last_sample;
            if(t==0){
                last_sample = cell.parent->sample[cell.parent->sample.size() - 1];
            }
            else{
                last_sample = cell.sample[t-1];
            }
            

            Eigen::VectorXd dummy = Eigen::VectorXd::Zero(4);

            Gaussian distr = cell.conditionals[t].evaluate(last_sample);
            cell.sample[t] = sample_gaussian(distr.m, distr.C);
            // cell.sample[t] = dummy;
            //sample from conditional
        }
    }
}

void sample_trajectory_recr(MOMAdata *cell){
    /*  
    * Recursive function
    */
    if (cell == nullptr)
        return;
    sc_sample_trajectory(*cell);
    sample_trajectory_recr(cell->daughter1);
    sample_trajectory_recr(cell->daughter2);
}


void sample_trajectory(std::vector<MOMAdata> &cells){
    /* applies to each cell going down the tree starting from all root cells */
    std::vector<MOMAdata *> p_roots = get_roots(cells);

    for(size_t i=0; i<p_roots.size(); ++i){
        sample_trajectory_recr(p_roots[i]);
    }
}


/* Output samples */

std::string outfile_name_samples(std::map<std::string, std::string> arguments, std::vector<Parameter_set>& params_list){
    /* Filename for a covariance file */
    std::string outfile = out_dir(arguments);
    outfile += file_base(arguments["infile"]);
    for (size_t i=0; i<params_list.size(); ++i){
        outfile += outfile_param_code(params_list[i]);
    }
    return outfile  + "_samples" + ".csv";
}

void setup_outfile_samples( std::string outfile, std::vector<Parameter_set> &params_list, int n_samples){    
    // same formating as prediction files
    for(size_t i=0; i<params_list.size(); ++i){
    if (i==0)
        params_list[i].to_csv(outfile);
    else
        params_list[i].to_csv(outfile, std::ios_base::app);
    }     

    std::ofstream file(outfile, std::ios_base::app);
    file << "\ncell_id,parent_id,time,sample,x,g,l,q\n";

    file.close();
}

void write_sample_to_file(std::string outfile, std::vector<MOMAdata> &cells, int sample_idx){
    std::ofstream file(outfile, std::ios_base::app);

    for(size_t i=0; i<cells.size();++i){
        for (size_t j=0; j<cells[i].time.size();++j ){
            file << cells[i].cell_id << "," << cells[i].parent_id << "," 
                 << cells[i].time[j] << "," << sample_idx << ",";
                 
            output_vector(file, cells[i].sample[j]);

            file << "\n"; 
        }
    }


    file.close();
}