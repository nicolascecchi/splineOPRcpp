#include <RcppEigen.h>

#include "SplineOP_constrained.h"
#include "SpeedEstimator.h"
// Continue with initSpeeds that has a problem with the shape and type


// CONSTRUCTOR
SplineOP_constrained::SplineOP_constrained(Eigen::MatrixXd data
                ,size_t nstates
                ,std::vector<int> sp
                ,double data_var
                ,int K
                ,int seed):
     nobs{static_cast<size_t>(data.cols())}
    ,ndims{static_cast<size_t>(data.rows())}
    ,K{K}
    ,sp{sp}
    ,nstates{nstates}   
    ,nspeeds{sp.size()}
    ,speeds(K, nobs, ndims, nstates) // initialize with dim nobs, its elements are Eigen::MatrixXd
    ,costs{K, nstates, data.cols()}//, std::numeric_limits<double>::infinity())
    ,initSpeeds{data.rows(),nspeeds}
    ,states() // default initialization, overwritten in the body of the constructor
    ,argmin_i{K, nstates, data.cols()}
    ,argmin_s{K, nstates, data.cols()}  
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


// Constructor with given speeds

/**
 * @brief Generates a 3D structure of states from observed data.
 * @param nstates Number of generated states per time step (rows).
 * @param data The observed data (ndims x nobs).
 * @param data_var Variance of the Gaussian noise.
 * @param seed Seed for random generation.
 * @return std::vector<Eigen::MatrixXd> The 3D state structure: 
 * vector index = time (t), Matrix = (nstates x ndims).
 */
std::vector<Eigen::MatrixXd> SplineOP_constrained::generate_states(
    size_t nstates,
    Eigen::MatrixXd data, // Input is now MatrixXd
    double data_var, 
    int seed)
{
    //size_t nobs = data.cols(); // Observations/Time (N) already known
    //size_t ndims = data.rows(); // Coordinates/Dimensions (K)
    size_t noise_states = nstates - 1; // Number of states requiring noise
    
    // 3D structure: vector index = time (t), Matrix = (ndims x nstates)
    std::vector<Eigen::MatrixXd> states_3d;
    states_3d.reserve(nobs);

    std::mt19937 gen(seed);
    double std_dev = std::sqrt(data_var);

    // Loop over time (observations)
    for (size_t t = 0; t < nobs; t++)
    {
        // 1. Create the 2D matrix for the current time slice (ndims x nstates)
        Eigen::MatrixXd current_time_states(ndims, nstates);

        // 2. State 0: The observed data point (The Mean)
        // Assign the observed data column (ndims x 1 vector) to the first state column.
        current_time_states.col(0) = data.col(t); 
        Eigen::MatrixXd matrix_of_noise = SplineOP_constrained::generate_matrix_of_noise(gen,std_dev,ndims,noise_states);
        
        // 3. Generate remaining states (1 to nstates-1)
        current_time_states.rightCols(noise_states) = matrix_of_noise.colwise() + data.col(t);
        states_3d.push_back(std::move(current_time_states));
    }
    return states_3d;
}

Eigen::MatrixXd SplineOP_constrained::generate_matrix_of_noise( //MON
    std::mt19937& gen, 
    double std_dev, 
    size_t rows, 
    size_t cols) 
{
    Eigen::MatrixXd matrix_of_noise(rows, cols);
    // Use N(0, std_dev) distribution
    std::normal_distribution<double> normal_dist(0.0, std_dev);
    
    // NOTE: This inner loop is unavoidable with std::normal_distribution.
    for (size_t j = 0; j < cols; ++j) {
        for (size_t i = 0; i < rows; ++i) {
            matrix_of_noise(i, j) = normal_dist(gen);
        }
    }
    return matrix_of_noise;
}

void SplineOP_constrained::predict(int K)
{
    Eigen::VectorXd v_s;
    Eigen::VectorXd v_t;
    Eigen::VectorXd p_s;
    Eigen::VectorXd p_t;
    Eigen::VectorXd best_speed;
    double interval_cost; // Cost of a 
    double candidate; // Cost of candidate value
    double current_MIN = std::numeric_limits<double>::infinity(); // placeholder for lowest cost so far
    int best_i = -1;  // Best previous state index 
    int best_s = -1;  // Best previous time index
    Eigen::Tensor<double, 1> tmp_speed_from_tensor(ndims); // placeholder to get current speeds
    K = K;
    // VOIR ICI .. VER ACA LA REINICIALIZACION
    // EN TOUT CAS, CA NE MARCHE PAS LE PREMIER TOUR NON PLUS
    speeds(K, nobs, ndims, nstates);
    costs{K, nstates, ndims};
    argmin_i{K, nstates, ndims};
    argmin_s{K, nstates, ndims};
    // reinitialize changepoints and costs for new fit (small overhead for first time)
    changepoints = std::vector(1,static_cast<int>(nobs-1)); 
    costs.setConstant(std::numeric_limits<double>::infinity());
    costs.chip(0,0).setConstant(0.0);
    
    for (int k=1; k<K+2; k++) // goes up to K+1 (segments) inclusive
    {   // Loop over data with k segments
        for (size_t t = 1; t < nobs; t++)
        { // current last point
            for (size_t j = 0; j < nstates; j++)
            { // current last state
                p_t = states[t].col(j).eval(); // Fix final position
                current_MIN = std::numeric_limits<double>::infinity();
                best_speed;
                best_i = -1;
                best_s = -1;
                // Find the best solution (state j,time t)
                for (size_t s = 0; s < t; s++)
                { // previous times
                    Rcpp::checkUserInterrupt(); // allow user interruption
                    for (size_t i = 0; i < nstates; i++)
                    { // previous state
                        // Fix start and end position in time and space
                        p_s = states[s].col(i).eval(); // Get starting state position
                        if (s == 0)
                        {
                            for (size_t spdidx = 0; spdidx < nspeeds; spdidx++){ // init speed loop
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
                        }
                        else
                        {
                            // compute speed
                            // sequential chipping. Recall that each cheap removes 1 dimension,
                            // thats why we always chip at dimension 0
                            auto chip_k = speeds.chip(k-1, 0); // get speeds from previous iteration
                            auto chip_state = chip_k.chip(j, 0); // take the ending state slice
                            auto chip_time = chip_state.chip(s,0); // at time s
                            tmp_speed_from_tensor = chip_time.eval();
                            v_s = Eigen::Map<Eigen::VectorXd>(tmp_speed_from_tensor.data(), tmp_speed_from_tensor.size());
                            v_t = 2*(p_t - p_s)/(t - s) - v_s; 
                            // Quadratic cost for interval [s, t)
                            // THIS INTERVAL COST IS BREAKING THE CODE 
                            interval_cost = qc.interval_cost(s, t, p_s, p_t, v_s);
                            // Candidate cost (DP recurrence)
                            candidate = costs(k-1, i, s) + interval_cost;         
                            if (candidate < current_MIN)
                            {
                                current_MIN = candidate;
                                best_speed = v_t;
                                best_i = i;
                                best_s = s;
                            }
                        }
                    }
                }
                // Update tables with best results.
                costs(k, j, t) = current_MIN;
                //update speeds
                auto speeds_chip_k = speeds.chip(k, 0); 
                auto speeds_chip_state = speeds_chip_k.chip(j, 0);
                auto speeds_chip_time = speeds_chip_state.chip(t,0);
                speeds_chip_time = Eigen::TensorMap<Eigen::Tensor<const double, 1>>(best_speed.data(), ndims);
                // update best params
                argmin_i(k, j, t) = best_i;
                argmin_s(k, j, t) = best_s;
            }   
        }     
    }
    SplineOP_constrained::backtrack_changes();
}

void SplineOP_constrained::backtrack_changes()
{
    // Find best final state
    double min_final = std::numeric_limits<double>::infinity();
    int best_final_state = -1;
    int curr_K = K+1;
    for (size_t j = 0; j < nstates; j++)
    {
        if (costs(curr_K, j, nobs - 1) < min_final)
        {
        min_final = costs(curr_K, j, nobs - 1);
        best_final_state = j;
        }
    }

    int j = best_final_state;
    int t = nobs - 1;
    // changepoints.push_back(t); // last point as a changepoint
    // Backtrack using argmin_s and argmin_i
    for (curr_K=K+1; curr_K>1; curr_K--)
    {
        std::cout << "Current K: " << curr_K << std::endl;
        std::cout << "Current best s: " << t << " current best j: " << j << std::endl;
        int s_prev = argmin_s(curr_K, j, t);
        int i_prev = argmin_i(curr_K, j, t);

        changepoints.push_back(s_prev);  // record change boundary

        t = s_prev;
        j = i_prev;
    }
    // Reverse the order to chronological (0 → nobs)
    std::reverse(changepoints.begin(), changepoints.end());
    for (size_t cpt = 0; cpt<changepoints.size(); cpt++)
    {
        changepoints[cpt] += 1;
    }
}


double SplineOP_constrained::get_segment_cost(int s, int t, Eigen::VectorXd p_s, Eigen::VectorXd p_t, Eigen::VectorXd v_s){
    return qc.interval_cost(s, t, p_s, p_t, v_s);
    }

void SplineOP_constrained::set_qc(Eigen::MatrixXd& data){
    qc = QuadraticCost(data);  // Initialize the cost function with the provided data
}

void SplineOP_constrained::set_states(std::vector<Eigen::MatrixXd> new_states){
    states = new_states;
    nstates = new_states[0].cols(); // assume that all time have the same nb of states
}
void SplineOP_constrained::set_initSpeeds(Eigen::MatrixXd new_initSpeeds){
    initSpeeds = new_initSpeeds;
    nspeeds = new_initSpeeds.cols();
}
