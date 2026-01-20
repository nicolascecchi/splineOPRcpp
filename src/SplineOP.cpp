#include <Eigen/Dense>

#include "SplineOP.h"
#include "SpeedEstimator.h"
// Continue with initSpeeds that has a problem with the shape and type


// Define the construction
SplineOP::SplineOP(Eigen::MatrixXd data
                ,int nstates
                ,std::vector<int> sp
                ,double data_var
                ,int seed):
     nobs{static_cast<int>(data.cols())}
    ,ndims{static_cast<int>(data.rows())}
    ,sp{sp}
    ,nstates{nstates}   
    ,nspeeds{sp.size()}
    ,speeds(nobs, Eigen::MatrixXd::Zero(data.rows(), nstates)) // initialize with dim nobs, its elements are Eigen::MatrixXd
    ,costs{nstates, data.cols()}//, std::numeric_limits<double>::infinity())
    ,pruning_costs{nstates, data.cols()}//, std::numeric_limits<double>::infinity())
    ,times_for_states(nstates, {0})
    ,initSpeeds{data.rows(),nspeeds}
    ,states() // default initialization, overwritten in the body of the constructor
    ,argmin_i{nstates, data.cols()}
    ,argmin_s{nstates, data.cols()}  
    ,qc{data}
    ,changepoints(1,static_cast<int>(data.cols())) //??? place holder, will be superseeded afterwards
    {
        costs.setConstant(std::numeric_limits<double>::infinity());
        argmin_i.setConstant(-1);
        argmin_s.setConstant(-1);

        states = generate_states(nstates, data, data_var, seed);
        //const std::vector<int> sp{20,30,40,50,60}; //         const std::vector<int> sp{20,30,40,50,60};
        initSpeeds = EstimateSpeeds(data, sp);
    }


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Eigen::MatrixXd SplineOP::generate_matrix_of_noise(std::mt19937& gen
                                                    , double std_dev 
                                                    , int rows 
                                                    , int cols) 
{
    Eigen::MatrixXd matrix_of_noise(rows, cols);
    // Use N(0, std_dev) distribution
    std::normal_distribution<double> normal_dist(0.0, std_dev);
    
    // NOTE: This inner loop is unavoidable with std::normal_distribution.
    for (int j = 0; j < cols; ++j) 
    {
        for (int i = 0; i < rows; ++i) 
        {
            matrix_of_noise(i, j) = normal_dist(gen);
        }
    }
    return matrix_of_noise;
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


/**
 * @brief Generates a 3D structure of states from observed data.
 * * @param nstates Number of generated states per time step (rows).
 * @param data The observed data (ndims x nobs).
 * @param data_var Variance of the Gaussian noise.
 * @param seed Seed for random generation.
 * @return std::vector<Eigen::MatrixXd> The 3D state structure: 
 * vector index = time (t), Matrix = (nstates x ndims).
 */
std::vector<Eigen::MatrixXd> SplineOP::generate_states(int nstates
                                                    , Eigen::MatrixXd data // Input is now MatrixXd
                                                    , double data_var 
                                                    , int seed)
{
    //int nobs = data.cols(); // Observations/Time (N) already known
    //int ndims = data.rows(); // Coordinates/Dimensions (K)
    int noise_states = nstates - 1; // Number of states requiring noise
    
    // 3D structure: vector index = time (t), Matrix = (ndims x nstates)
    std::vector<Eigen::MatrixXd> states_3d;
    states_3d.reserve(nobs);

    std::mt19937 gen(seed);
    double std_dev = std::sqrt(data_var);

    // Loop over time (observations)
    for (int t = 0; t < nobs; t++)
    {
        // 1. Create the 2D matrix for the current time slice (ndims x nstates)
        Eigen::MatrixXd current_time_states(ndims, nstates);

        // 2. State 0: The observed data point (The Mean)
        // Assign the observed data column (ndims x 1 vector) to the first state column.
        current_time_states.col(0) = data.col(t); 
        Eigen::MatrixXd matrix_of_noise = SplineOP::generate_matrix_of_noise(gen,std_dev,ndims,noise_states);
        
        // 3. Generate remaining states (1 to nstates-1)
        current_time_states.rightCols(noise_states) = matrix_of_noise.colwise() + data.col(t);
        states_3d.push_back(std::move(current_time_states));
    }
    return states_3d;
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Penalized version
void SplineOP::predict(double beta)
{
    const double INF = std::numeric_limits<double>::infinity();
    // reinitialize changepoints and costs for new fit (small overhead for first time)
    changepoints = std::vector(1,static_cast<int>(nobs-1)); 
    costs.setConstant(INF);
    argmin_i.setConstant(-1);
    argmin_s.setConstant(-1);

    // Preallocate temporaries (all have ndims)
    Eigen::VectorXd p_t; // position at time t
    Eigen::VectorXd p_s; // position at time s
    Eigen::VectorXd v_s; // speeds at time s
    Eigen::VectorXd p_s_best(ndims); // best initial position
    Eigen::VectorXd v_s_best(ndims); // best initial speed
    
    // Variables for optimization
    double interval_cost; // Cost of interval being considered
    double candidate; // Total cost of the interval being considered
    double current_MIN = INF; // Current best cost
    int best_i = -1; // Current best position index
    int best_s = -1; // Current best previos time index
    int best_speed_idx = -1; // Best initial speed idx, only meaningful
                             // if best_s = 0.

    // Loop over end times
    for (int t = 1; t < nobs; t++)
    {
        // Get states for the ending time being evaluated
        // Loop over indexes of the end states at time t
        for (int j = 0; j < nstates; j++)
        { // current last state
            p_t = states[t].col(j); // Fix final position
            
            // Set optimization variables to defaults
            // We need to reset for each ending state j
            current_MIN = INF;
            best_speed_idx = -1;
            best_i = -1;
            best_s = -1;
            
            // Treat the case of s = 0 explicitly
            // Avoid evaluating the if(s==0) t times, when it is only true once.
            int s = 0;
            // Loop over indexes of initial state
            for (int i = 0; i < nstates; i++)
            { 
                p_s = states[s].col(i); // Fix start position in space
                // Loop over initial speeds
                for (int spdidx = 0; spdidx < nspeeds; spdidx++)
                { 
                    // Set current initial speed
                    v_s = initSpeeds.col(spdidx);

                    // Quadratic cost for interval [s, t)
                    interval_cost = qc.interval_cost(s, t, p_s, p_t, v_s);
                    candidate     = interval_cost; // Candidate cost (DP recurrence)
                    
                    if (candidate < current_MIN)
                    {
                        current_MIN = candidate;
                        best_speed_idx = spdidx;
                        best_i = i;
                        best_s = s;
                    }
                }                        
            }
            // Find the best solution (state j,time t)
            for (int s = 1; s < t; s++)
            { // previous times
                //Rcpp::checkUserInterrupt(); // allow user interruption

                // Loop over previous-time state indexes
                for (int i = 0; i < nstates; i++)
                { // previous state
                    p_s = states[s].col(i); // Get starting state position
                    v_s = speeds[s].col(i); // Get starting speed

                    // Quadratic cost for interval [s, t)
                    interval_cost = qc.interval_cost(s, t, p_s, p_t, v_s);
                    candidate     = costs(i, s) + interval_cost + beta; // Candidate cost (DP recurrence)

                    if (candidate < current_MIN)
                    {
                        current_MIN    = candidate;
                        best_i         = i;
                        best_s         = s;
                        best_speed_idx = -1;
                    }
                }
            }
            // Compute v_t ONCE using the best predecessor
            p_s_best = states[best_s].col(best_i);
            if (best_s == 0)
            {
                v_s_best = initSpeeds.col(best_speed_idx);
            }
            else
            {
                v_s_best = speeds[best_s].col(best_i);
            }

            const double inv_dt = 1.0 / static_cast<double>(t - best_s);
            // v_t = 2 * (p_t - p_s) / (t - s) - v_s
            speeds[t].col(j) = 2.0 * inv_dt * (p_t - p_s_best) - v_s_best;
                    
            // Store DP values
            costs(j, t) = current_MIN;
            argmin_i(j, t) = best_i;
            argmin_s(j, t) = best_s;
        }   
    }
    SplineOP::backtrack_changes();
}

void SplineOP::pruningv1(double beta, double margin)
{
    Eigen::VectorXd v_s;
    Eigen::VectorXd v_t;
    Eigen::VectorXd p_s;
    Eigen::VectorXd p_t;
    Eigen::VectorXd best_speed;
    double interval_cost;
    double candidate;
    double current_MIN = std::numeric_limits<double>::infinity();
    int best_i = -1;
    int best_s = -1;
    // reinitialize changepoints and costs for new fit (small overhead for first time)
    changepoints = std::vector(1,static_cast<int>(nobs-1)); 
    costs.setConstant(std::numeric_limits<double>::infinity());
    pruning_costs.setConstant(std::numeric_limits<double>::infinity());
    std::vector<std::vector<int>> initial_value(nstates, std::vector<int>{0});
    times_for_states.assign(initial_value.begin(), initial_value.end());

    // Big loop that goes through time
    for (int t = 1; t < nobs; t++)
    { // Loop over the ending states at current time t 
        for (int j = 0; j < nstates; j++)
        {
            p_t = states[t].col(j).eval(); // Fix final position
            current_MIN = std::numeric_limits<double>::infinity(); // set up variables used for optimization
            best_speed;
            best_i = -1;
            best_s = -1;
            // Double loop over all accessible previous time-state pairs
            // finds the best accessible solution for Q(j,t)
            // Loop 1: Over all the states
            for (int i=0; i<nstates; i++)
            {   // Loop 2: Over time positions that are still alive for that state index
                for (int s : times_for_states[i])
                {
                    p_s = states[s].col(i).eval(); // Get starting state position
                    if (s == 0)
                    {
                        for (int spdidx = 0; spdidx < nspeeds; spdidx++){ // init speed loop
                            v_s = initSpeeds.col(spdidx).eval();
                            v_t = 2*(p_t - p_s)/(t - s) - v_s; // simple slope rule
                            // Quadratic cost for interval [s, t)
                            interval_cost = qc.interval_cost(s, t, p_s, p_t, v_s);
                            // Candidate cost (DP recurrence)
                            candidate = interval_cost;
                            if (candidate < current_MIN){
                                current_MIN = candidate;
                                best_speed = v_t;
                                best_i = i;
                                best_s = s;
                            }
                        }                        
                    } // end of cost without changepoint before
                    else
                    {
                        // compute speed
                        v_s = speeds[s].col(i).eval();
                        v_t = 2*(p_t - p_s)/(t - s) - v_s; 
                        // Quadratic cost for interval [s, t)
                        // THIS INTERVAL COST IS BREAKING THE CODE 
                        interval_cost = qc.interval_cost(s, t, p_s, p_t, v_s);
                        // Candidate cost (DP recurrence)
                        candidate = costs(i, s) + interval_cost + beta;         
                        if (candidate < current_MIN)
                        {
                            current_MIN = candidate;
                            best_speed = v_t;
                            best_i = i;
                            best_s = s;
                        }
                    } // end of cost with at least 1 changepoint before
                } // end of times for a given position index
            } // end of position index
            
            // Fill best found solutions
            // std::cout<< "j: " << j << " t: " << t <<std::endl; 
            // std::cout<< "curr min: " << current_MIN << "best sp: " << best_speed <<std::endl; 
            costs(j, t) = current_MIN;
            speeds[t].col(j) = best_speed;
            argmin_i(j, t) = best_i;
            argmin_s(j, t) = best_s;
        } // end of state j for time t 
        prunev1(t, margin); 
    } // end of t
    SplineOP::backtrack_changes();
}


// Penalized version
void SplineOP::pruningv2(double beta)
{
    Eigen::VectorXd v_s;
    Eigen::VectorXd v_t;
    Eigen::VectorXd p_s;
    Eigen::VectorXd p_t;
    Eigen::VectorXd best_speed;
    double interval_cost;
    double candidate;
    double current_MIN = std::numeric_limits<double>::infinity();
    int best_i = -1;
    int best_s = -1;
    // reinitialize changepoints and costs for new fit (small overhead for first time)
    changepoints = std::vector(1,static_cast<int>(nobs-1)); 
    costs.setConstant(std::numeric_limits<double>::infinity());
    pruning_costs.setConstant(std::numeric_limits<double>::infinity());
    std::vector<std::vector<int>> initial_value(nstates, std::vector<int>{0});
    times_for_states.assign(initial_value.begin(), initial_value.end());

    // Big loop that goes through time
    for (int t = 1; t < nobs; t++)
    { // Loop over the ending states at current time t 
        for (int j = 0; j < nstates; j++)
        {
            p_t = states[t].col(j).eval(); // Fix final position
            current_MIN = std::numeric_limits<double>::infinity(); // set up variables used for optimization
            best_speed;
            best_i = -1;
            best_s = -1;
            // Double loop over all accessible previous time-state pairs
            // finds the best accessible solution for Q(j,t)
            // Loop 1: Over all the states
            for (int i=0; i<nstates; i++)
            {   // Loop 2: Over time positions that are still alive for that state index
                for (int s : times_for_states[i])
                {
                    p_s = states[s].col(i).eval(); // Get starting state position
                    if (s == 0)
                    {
                        for (int spdidx = 0; spdidx < nspeeds; spdidx++){ // init speed loop
                            v_s = initSpeeds.col(spdidx).eval();
                            v_t = 2*(p_t - p_s)/(t - s) - v_s; // simple slope rule
                            // Quadratic cost for interval [s, t)
                            interval_cost = qc.interval_cost(s, t, p_s, p_t, v_s);
                            // Candidate cost (DP recurrence)
                            candidate = interval_cost;
                            if (candidate < current_MIN){
                                current_MIN = candidate;
                                best_speed = v_t;
                                best_i = i;
                                best_s = s;
                            }
                            if (candidate < pruning_costs(i,s))
                            {
                            pruning_costs(i,s) = candidate;
                            }
                        }                        
                    } // end of cost without changepoint before
                    else
                    {
                        // compute speed
                        v_s = speeds[s].col(i).eval();
                        v_t = 2*(p_t - p_s)/(t - s) - v_s; 
                        // Quadratic cost for interval [s, t)
                        // THIS INTERVAL COST IS BREAKING THE CODE 
                        interval_cost = qc.interval_cost(s, t, p_s, p_t, v_s);
                        // Candidate cost (DP recurrence)
                        candidate = costs(i, s) + interval_cost + beta;         
                        if (candidate < current_MIN)
                        {
                            current_MIN = candidate;
                            best_speed = v_t;
                            best_i = i;
                            best_s = s;
                        }
                        if ((candidate-beta) < pruning_costs(i,s))
                        {
                            pruning_costs(i,s) = candidate;
                        }
                    } // end of cost with at least 1 changepoint before
                } // end of times for a given position index
            } // end of position index
            
            // Fill best found solutions
            // std::cout<< "j: " << j << " t: " << t <<std::endl; 
            // std::cout<< "curr min: " << current_MIN << "best sp: " << best_speed <<std::endl; 
            costs(j, t) = current_MIN;
            speeds[t].col(j) = best_speed;
            argmin_i(j, t) = best_i;
            argmin_s(j, t) = best_s;
        } // end of state j for time t 
        prunev2(t); 
    } // end of t
    SplineOP::backtrack_changes();
}

void SplineOP::prunev1(int t, double margin)
{
    Eigen::VectorXd v_s;
    Eigen::VectorXd v_t;
    Eigen::VectorXd p_s;
    Eigen::VectorXd p_t;
    double interval_cost;
    double candidate;
    double best_cost_jt;
    int counter;

    // PRUNING
        //std::cout << "=============================================  " << t << "  =============================================" << std::endl;
    // parcourir all the state indexes 
    for (int i = 0; i < nstates; i++) 
    {
        std::vector<int>& times_vec = times_for_states[i];
        for (auto it = times_vec.begin(); it != times_vec.end(); ) 
        {   
            int s = *it;
            counter = 0 ; // reinitialize counter for each (i, s)
            if (s == 0)
                { // For the moment leave this like that, need a fix in the future to also prune
                    //for (int spdidx = 0; spdidx < nspeeds; spdidx++)
                    //{ // init speed loop
                    //    v_s = initSpeeds.col(spdidx).eval();
                    //    v_t = 2*(p_t - p_s)/(t - s) - v_s; // simple slope rule
                    //    // Quadratic cost for interval [s, t)
                    //    interval_cost = qc.interval_cost(s, t, p_s, p_t, v_s);
                    //    // Candidate cost (DP recurrence)
                    //    candidate = interval_cost;
                    //    if (candidate < pruning_costs(i,s))
                    //    {
                    //    pruning_costs(i,s) = candidate;
                    //    }
                    //}                        
                } // end of cost without changepoint before
                else
                {
                    // Evaluate all ending states
                    p_s = states[s].col(i).eval();
                    v_s = speeds[s].col(i).eval();
                    for(int j = 0; j<nstates; j++)
                    {    
                        p_t = states[t].col(j).eval();
                        interval_cost = qc.interval_cost(s, t, p_s, p_t, v_s);
                        candidate = costs(i, s) + interval_cost;         
                        if (candidate > costs(j,t) * (1 + margin))
                        {counter++;}
                    }
                    
                }
                // Get the time index
            if (counter == nstates) 
            {
                it = times_vec.erase(it); 
            } 
            else 
            {
                ++it;
            }
        }
        times_vec.push_back(t); 
    }
}


void SplineOP::prunev2(int t)
{
    // PRUNING
        std::cout << "===============================================================  " << t << "===============================================================  " << std::endl;
        double min_jt = costs.col(t).minCoeff(); // Best cost at time t
        // parcourir all the state indexes 
        for (int i = 0; i < nstates; i++) 
        {
            std::cout << "Pruning for state idxs " << i << std::endl;
            std::vector<int>& times_vec = times_for_states[i];
            for (auto it = times_vec.begin(); it != times_vec.end(); ) 
            {
                int s = *it; // Get the time index
                if (min_jt < pruning_costs(i, s)) 
                {
                    std::cout << "At time "<< t+1 << " prune time " << s+1 << std::endl;
                    std::cout << " with pruning cost " <<pruning_costs(i,s) << " > min " << min_jt << std::endl;
                    it = times_vec.erase(it); 
                } 
                else 
                {
                    ++it;
                }
            }
            times_vec.push_back(t); 
        }
        pruning_costs.setConstant(std::numeric_limits<double>::infinity());
}

void SplineOP::backtrack_changes()
{
    // Find best final state
    double min_final = std::numeric_limits<double>::infinity();
    int best_final_state = -1;
    for (int j = 0; j < nstates; j++)
    {
        if (costs(j, nobs - 1) < min_final)
        {
        min_final = costs(j, nobs - 1);
        best_final_state = j;
        }
    }

    int j = best_final_state;
    int t = nobs - 1;

    // Backtrack using argmin_s and argmin_i
    while (true)
    {
    int s_prev = argmin_s(j, t);
    int i_prev = argmin_i(j, t);

    if (s_prev == 0)
        break;  // reached the beginning or invalid index

    changepoints.push_back(s_prev);  // record change boundary

    t = s_prev;
    j = i_prev;
    }
    // Reverse the order to chronological (0 → nobs)
    std::reverse(changepoints.begin(), changepoints.end());
}


double SplineOP::get_segment_cost(int s, int t, Eigen::VectorXd p_s, Eigen::VectorXd p_t, Eigen::VectorXd v_s){
    return qc.interval_cost(s, t, p_s, p_t, v_s);
    }

void SplineOP::set_qc(Eigen::MatrixXd& data){
    qc = QuadraticCost(data);  // Initialize the cost function with the provided data
}
void SplineOP::set_states(std::vector<Eigen::MatrixXd> new_states){
    states = new_states;
    nstates = new_states[0].cols(); // assume that all time have the same nb of states
}
void SplineOP::set_initSpeeds(Eigen::MatrixXd new_initSpeeds){
    initSpeeds = new_initSpeeds;
    nspeeds = new_initSpeeds.cols();
}
