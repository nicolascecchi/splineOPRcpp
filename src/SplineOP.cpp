#include "SplineOP.h"
#include <RcppEigen.h>
#include "SpeedEstimator.h"
// Continue with initSpeeds that has a problem with the shape and type


// Define the construction
SplineOP::SplineOP(Eigen::MatrixXd  data
                ,size_t nstates
                ,size_t nspeeds
                ,double data_var
                ,int seed):
     nobs{static_cast<int>(data.cols())}
    ,nstates{nstates}   
    ,nspeeds{nspeeds}

    ,speeds(nobs, Eigen::MatrixXd::Zero(data.rows(), data.cols())) // initialize with dim nobs, its elements are Eigen::MatrixXd
    ,costs{nstates, data.cols()}//, std::numeric_limits<double>::infinity())
    ,initSpeeds{data.rows(),nspeeds}
    ,states() // default initialization, overwritten in the body of the constructor

    ,argmin_i{nstates, data.cols()}
    ,argmin_s{nstates, data.cols()}  
    ,qc{data}

    ,changepoints(1,static_cast<int>(data.size())-1) // place holder, will be superseeded afterwards
    {
        this->costs.setConstant(std::numeric_limits<double>::infinity());
        this->argmin_i.setConstant(-1);
        this->argmin_s.setConstant(-1);

        this->states = generate_states(nstates, data, data_var, seed);
        const std::vector<int> sp{20,40,60};
        this->initSpeeds = EstimateSpeeds(data, sp);
       }


// Constructor with given speeds

/**
 * @brief Generates a 3D structure of states from observed data.
 * * @param nstates Number of generated states per time step (rows).
 * @param data The observed data (ndims x nobs).
 * @param data_var Variance of the Gaussian noise.
 * @param seed Seed for random generation.
 * @return std::vector<Eigen::MatrixXd> The 3D state structure: 
 * vector index = time (t), Matrix = (nstates x ndims).
 */
std::vector<Eigen::MatrixXd> SplineOP::generate_states(
    size_t nstates,
    Eigen::MatrixXd data, // Input is now MatrixXd
    double data_var, 
    int seed)
{
    Eigen::Index nobs = data.cols(); // Observations/Time (N)
    Eigen::Index ndims = data.rows(); // Coordinates/Dimensions (K)
    Eigen::Index noise_states = nstates - 1; // Number of states requiring noise
    
    // 3D structure: vector index = time (t), Matrix = (ndims x nstates)
    std::vector<Eigen::MatrixXd> states_3d;
    states_3d.reserve(nobs);

    std::mt19937 gen(seed);
    double std_dev = std::sqrt(data_var);

    // Loop over time (observations)
    for (Eigen::Index t = 0; t < nobs; t++){
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

Eigen::MatrixXd SplineOP::generate_matrix_of_noise( //MON
    std::mt19937& gen, 
    double std_dev, 
    Eigen::Index rows, 
    Eigen::Index cols) 
{
    Eigen::MatrixXd matrix_of_noise(rows, cols);
    // Use N(0, std_dev) distribution
    std::normal_distribution<double> normal_dist(0.0, std_dev);
    
    // NOTE: This inner loop is unavoidable with std::normal_distribution.
    for (Eigen::Index j = 0; j < cols; ++j) {
        for (Eigen::Index i = 0; i < rows; ++i) {
            matrix_of_noise(i, j) = normal_dist(gen);
        }
    }
    return matrix_of_noise;
}

void SplineOP::predict(double beta){

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
    this->changepoints = std::vector(1,this->nobs-1); 
    this->costs.setConstant(std::numeric_limits<double>::infinity());

    for (size_t j = 0; j < this->nstates; j++){
        this->costs(j, 0) = 0; // random initial cost for 0 data point
    }  
    // Loop over data
    for (size_t t = 1; t < static_cast<size_t>(nobs); t++){ // current last point
        for (size_t j = 0; j < nstates; j++){ // current state
            p_t = states[t].col(j).eval(); // Fix final position
            current_MIN = std::numeric_limits<double>::infinity();
            best_speed;
            best_i = -1;
            best_s = -1;
            // Find the best solution (state j,time t)
            for (size_t s = 0; s < t; s++){ // previous times
                Rcpp::checkUserInterrupt(); // allow user interruption
                for (size_t i = 0; i < nstates; i++){ // previous state
                    // Fix start and end position in time and space
                    p_s = states[s].col(i).eval(); // Get starting state position
                    if (s == 0){
                        //Rcpp::Rcout << "Solution without change" << std::endl;
                        for (size_t spdidx = 0; spdidx < this->nspeeds; spdidx++){ // init speed loop
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
                    else{
                        // compute speed
                        //Rcpp::Rcout << "Changepoint time : " << s << std::endl;
                        v_s = this->speeds[s].col(i).eval();
                        v_t = 2*(p_t - p_s)/(t - s) - v_s; 
                        // Quadratic cost for interval [s, t)
                        // THIS INTERVAL COST IS BREAKING THE CODE 
                        interval_cost = qc.interval_cost(s, t, p_s, p_t, v_s);
                        // Candidate cost (DP recurrence)
                        candidate = costs(i, s) + interval_cost + beta;         
                        if (candidate < current_MIN){
                            current_MIN = candidate;
                            best_speed = v_t;
                            best_i = i;
                            best_s = s;
                              }
                     }
                }
            }

            costs(j, t) = current_MIN;
            this->speeds[t].col(j) = best_speed;
            argmin_i(j, t) = best_i;
            argmin_s(j, t) = best_s;
        }
    }
    SplineOP::backtrack_changes();
}

void SplineOP::backtrack_changes(){
    // Find best final state
    double min_final = std::numeric_limits<double>::infinity();
    int best_final_state = -1;
    for (size_t j = 0; j < this->nstates; j++){
        if (costs(j, this->nobs - 1) < min_final)
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

    if (s_prev < 0)
        break;  // reached the beginning or invalid index

    changepoints.push_back(s_prev);  // record change boundary

    t = s_prev;
    j = i_prev;
    }
    // Reverse the order to chronological (0 → nobs)
    std::reverse(changepoints.begin(), changepoints.end());
}

double SplineOP::get_segment_cost(int s, int t, Eigen::VectorXd p_s, Eigen::VectorXd p_t, Eigen::VectorXd v_s){
    return this->qc.interval_cost(s, t, p_s, p_t, v_s);
    }

void SplineOP::set_qc(Eigen::MatrixXd& data){
    this->qc = QuadraticCost(data);  // Initialize the cost function with the provided data
}
