//--------------------------------------------
//  EMBasins.h
//
//  Created by Jason Prentice on November 13, 2012.
//  Modified by Adrianna Loback.
//  Mexed version checked for local use
//  on December 8, 2018.
//
//--------------------------------------------

#ifndef ____EMBasins__
#define ____EMBasins__

#include <gsl/gsl_rng.h>

#include <vector>
#include <string>
#include <map>

using namespace std;


// ************ RNG ***************
class RNG
{
public:
    RNG();
    ~RNG();
    int discrete(const vector<double>&);
    bool bernoulli(double);
    vector<int> randperm(int);
private:
    gsl_rng* rng_pr;
};

struct State
{
    vector<int> active_constraints;
    vector<int> on_neurons;
    vector<double> P;
    vector<double> weight;
    double freq;
    double pred_prob;

    vector<char> word;
    
    int identifier;
};
// *********************************
typedef map<string,State>::iterator state_iter;
typedef map<string,State>::const_iterator const_state_iter;

// ************ Spike ***************
struct Spike
{
    int bin;
    int neuron_ind;
};
// *********************************
// ************ SpikeComparison ***************
class SpikeComparison
{
public:
    bool operator() (const Spike& lhs, const Spike& rhs) const
    {
        return (lhs.bin < rhs.bin);
    }
};
// *********************************

//typedef priority_queue<Spike, vector<Spike>, SpikeComparison> SpikeHeap;

// ************ EMBasins ***************
class paramsStruct;

template <class BasinT>
class EMBasins
{
public:
    EMBasins(int N, int nbasins);
    EMBasins(vector<vector<double> >& st, double binsize, int nbasins);
    ~EMBasins();
    
    vector<double> train(int niter);
    vector<double> crossval(int niter, int k);      // k-fold cross-validation
    vector<double> test(const vector<vector<double> >& st, double binsize);
    
    int nstates() const {return all_states.size();};
    vector<unsigned long> state_hist() const;
    vector<unsigned long> test_hist() const;
    vector<double> all_prob() const;
    vector<double> test_prob() const;
    vector<double> P() const;    
    vector<paramsStruct> basin_params();
    vector<char> sample(int);
    vector<char> word_list();
    
    vector<double> w;       // 1 x nbasins
    vector<double> m;       // N x nbasins
    
    vector<double> test_logli;
    
protected:
    int nbasins;
    int N;
    double nsamples;
    
    map<string, State> all_states;
    map<string, State> train_states;
    map<string, State> test_states;
    
    vector<BasinT> basins;
    
    RNG* rng;
    
    vector<string> raster;
    
    void update_w();
    double update_P();
    double update_P_test();
    
    double set_state_P(State&);
    vector<Spike> sort_spikes(const vector<vector<double> >&, double) const;
    
};

// *********************************

// ************ HMM ***************
template <class BasinT>
class HMM : public EMBasins<BasinT>
{
public:
    HMM(vector<vector<double> >& st, vector<double>, vector<double>, double binsize, int nbasins);
    
    vector<double> train(int niter);
    vector<int> viterbi(bool);
    
//    vector<char> get_raster();
    vector<double> emiss_prob();
    vector<double> get_forward();
    vector<double> get_backward();
    vector<double> get_P();
    vector<double> get_trans();
    vector<double> stationary_prob();
    pair<vector<double>, vector<double> > pred_prob();
    vector<int> state_v_time();
    
    vector<char> sample(int);
    vector<double> P_indep();
protected:
    int T, tmax, tskip;

   
    vector<double> forward;         // Forward filtering distribution
    vector<double> backward;        // Backward filtering distribution
    //vector<string> words;

    vector<State*> state_list;
    
    void update_forward();
    void update_backward();

    void forward_backward();
    vector<double> trans_at_t(int);
    
    void update_P();

    double logli(bool);
private:
    vector<double> w0;
    vector<double> trans;           // State transition probability matrix
    
    void update_trans();
    vector<double> emiss_obs(bool,int);

};
// *********************************

// ************ Autocorr ***********
template <class BasinT>
class Autocorr : public HMM<BasinT>
{
public:
    Autocorr(vector<vector<double> >& st, double binsize, int nbasins);
    ~Autocorr();
    
    vector<double> train(int niter);
    vector<int> viterbi();
    vector<double> get_basin_trans();
    
private:
    
    void update_forward();
    void update_backward();
    void update_P();
    void update_w();
    void update_basin_trans_indep();
    
    vector<double> trans_at_t(int);
    
    vector<double*> basin_trans;
    double logli();
};




#endif /* defined(____EMBasins__) */
