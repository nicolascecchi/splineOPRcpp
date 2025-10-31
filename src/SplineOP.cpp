#include "SplineOP.h"
#include "Matrix.h"

// Define the construction
SplineOP::SplineOP(std::vector<double>  data
                ,size_t nstates
                ,size_t nspeeds
                ,double data_var
                ,int seed):
    nobs{static_cast<int>(data.size())}
    ,nstates{nstates}   
    ,nspeeds{nspeeds}

    ,speeds(nstates, data.size(), 0.0)
    ,costs(nstates, data.size(), std::numeric_limits<double>::infinity())
    ,initspeeds(nspeeds,size_t{1},-0.017485295)
    ,states(nstates, data.size(), 0.) // place holder for initialization

    ,argmin_i(nstates, data.size(), -1)
    ,argmin_s(nstates, data.size(), -1)  
    ,qc(data)

    ,changepoints(1,static_cast<int>(data.size())) // place holder, will be superseeded afterwards
    {
        this->states = generate_states(nstates, data, data_var, seed);
    }

spop::Matrix<double> SplineOP::generate_states(size_t nstates,const std::vector<double>& data,double data_var, int seed){
    spop::Matrix<double> states(nstates, this->nobs, 0.0);
    //std::random_device rd;
    std::mt19937 gen(seed);

    for (int t = 0; t < nobs; t++)
    {
        // first row = data[t]
        states(0, t) = data[t];

        // remaining rows: Gaussian centered on data[j]
        std::normal_distribution<double> normal_dist(data[t], std::sqrt(data_var));
        for (size_t j = 1; j < nstates; j++)
        {
        states(j, t) = normal_dist(gen);
        }
    }
    return states;
}

void SplineOP::predict(double beta){
    for (size_t j = 0; j < this->nstates; j++){
    this->costs(j, 0) = 0;  // random initial cost for 0 data point
    }
    for (size_t t = 1; t < static_cast<size_t>(nobs); t++)
    {
        for (size_t j = 0; j < nstates; j++)
        { // current state
        double current_MIN = std::numeric_limits<double>::infinity();
        double best_speed;
        int best_i = -1;
        int best_s = -1;

        for (size_t s = 0; s < t; s++){ // previous time
            for (size_t i = 0; i < nstates; i++){ // previous state
            // Fix start and end position in time and space
            double p_s = states(i, s);
            double p_t = states(j, t);
            if (s == 0){
                for (size_t j = 0; j < this->nspeeds; j++){ // init speed loop
                double v_t = 2*(p_t - p_s)/(t - s) - initspeeds(j,0); // simple slope rule
                // Quadratic cost for interval [s, t)
                double interval_cost = qc.quadratic_cost_interval(s, t, p_s, p_t, v_t);
                // Candidate cost (DP recurrence)
                double candidate = interval_cost;

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
                double v_t = 2*(p_t - p_s)/(t - s) - speeds(i, s); 
                
                // Quadratic cost for interval [s, t)
                double interval_cost = qc.quadratic_cost_interval(s, t, p_s, p_t, v_t);
                
                // Candidate cost (DP recurrence)
                double candidate = costs(i, s) + interval_cost + beta;
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
        speeds(j, t) = best_speed;
        argmin_i(j, t) = best_i;
        argmin_s(j, t) = best_s;
        }
    }
    SplineOP::backtrack_changes();
}

void SplineOP::backtrack_changes(){
    std::vector<int> changepoints;
    changepoints.push_back(nobs - 1);  // last time index is always a change point
    // Find best final state
    double min_final = std::numeric_limits<double>::infinity();
    int best_final_state = -1;
    for (size_t j = 0; j < this->nstates; j++){
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

    if (s_prev < 0)
        break;  // reached the beginning or invalid index

    changepoints.push_back(s_prev);  // record change boundary

    t = s_prev;
    j = i_prev;
    }
    // Reverse the order to chronological (0 → nobs)
    std::reverse(changepoints.begin(), changepoints.end());
}

