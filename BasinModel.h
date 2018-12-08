//--------------------------------------------
//  BasinModel.h
//  
//
//  Created by Jason Prentice on 11/22/12.
//
//--------------------------------------------

#ifndef ____BasinModel__
#define ____BasinModel__

#include <iostream>
#include <vector>
#include <map>
#include <string>


using namespace std;

struct State;       // Defined in EMBasins.h
class RNG;

// *********************** myMatrix ****************************
template <class T>
class myMatrix {
public:
    myMatrix();
    myMatrix(vector<T>&, int, int);
    void assign(vector<T>&, int, int);
    int get_N();
    int get_M();
    
    const T& at(const int i, const int j) const;     // Subscripted index
    const T& at(const int i) const;            // Linear index
    
    vector<T>* data();
    
private:
    vector<T> matrix_data;
    int N,M;        // N = #rows, M = #columns
};

// myMatrix definition

template <class T>
myMatrix<T>::myMatrix() : N(0),M(0) {};

template <class T>
myMatrix<T>::myMatrix(vector<T>& _data, int _N, int _M) : N(_N), M(_M) {
    if (_data.size() != _N*_M) {
        cerr << "Matrix dimensions must agree." << endl;
        N = 0; M = 0;
    } else {
        matrix_data = _data;
    }
}

template <class T>
void myMatrix<T>::assign(vector<T>& _data, int _N, int _M) {
    if (_data.size() != _N*_M) {
        cerr << "Matrix dimensions must agree." << endl;
        N = 0; M = 0;
    } else {
        N = _N; M = _M;
        matrix_data = _data;
    }
    return;
}

template <class T>
int myMatrix<T>::get_N() { return N; }

template <class T>
int myMatrix<T>::get_M() { return M; }

template <class T>
const T& myMatrix<T>::at(const int i, const int j) const {
    if (i < 0 || i >= N || j < 0 || j >= M) {
        cerr << "Index exceeds matrix dimensions." << endl;
    } else {
        return matrix_data[j*N + i];
    }
}
template <class T>
const T& myMatrix<T>::at(const int i) const {
    if (i<0 || i>= N*M) {
        cerr << "Index exceeds matrix dimensions." << endl;
    } else {
        return matrix_data[i];
    }
}

template <class T>
vector<T>* myMatrix<T>::data() { return &matrix_data; }
// *****************************************************************

// *********************** paramsStruct ****************************
class paramsStruct {            // Variable-sized structure holding double matrices
public:
    
    paramsStruct();
    int get_nfields();
    void addField(string, myMatrix<double>&);
    const char** fieldNamesArray();
    vector<double>* getFieldData(int);
    int getFieldN(int);
    int getFieldM(int);
    const char* getFieldName(int);
    
private:
    int nfields;
    map<string, myMatrix<double>* > fields;
    vector<const char*> fieldNames;
    
};
// *****************************************************************

// ************************** BasinModel ************************
class BasinModel
{
public:
    BasinModel(int N, int basin_num, RNG* rng) : N(N), basin_num(basin_num), rng(rng) {};

    void reset_stats();
    void increment_stats(const State&);
    void normalize_stats();

    double get_norm() const {return norm;};
    
//    int nparams() const;

protected:
    int N;
    int basin_num;
    double norm;
    
    vector<double> stats;
    RNG* rng;
};
// ***************************************************************

// ************************** IndependentBasin ************************

class IndependentBasin : public BasinModel
{
public:
    IndependentBasin(int,int,RNG*);
    
    static vector<int> get_active_constraints(const State&);    
    
    void doMLE(double);
    double P_state(const State&) const;
    vector<char> sample();
    
    paramsStruct get_params();
//    int nparams() const;
//    vector<double> get_params() const;
private:
//    vector<double> m;
    myMatrix<double> m;
    vector<char> above_thresh_bool;
    vector<int> above_thresh_list;
    
    double prefactor;
    
    void update_thresh_list();
};
// ***************************************************************


#endif /* defined(____BasinModel__) */
