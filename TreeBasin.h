//--------------------------------------------
//  TreeBasin.h
//  
//
//  Created by Jason Prentice on 11/28/12.
//
//--------------------------------------------

#ifndef _TreeBasin_h
#define _TreeBasin_h

#include "BasinModel.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <vector>

class RNG; //Defined in EMBasins.h

// Define a graph type that associates a weight with each
// edge. We store the weights using internal properties as described
// in BGL.
struct EdgeProperty {
    double negI;        // negative mutual information; weight used to find min spanning tree
};


typedef boost::adjacency_list<boost::setS,
boost::vecS,
boost::undirectedS,
boost::no_property,
EdgeProperty>                 Graph;

typedef boost::graph_traits<Graph>::vertex_descriptor       Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor         Edge;

struct TreeEdgeProb {
    TreeEdgeProb(double p10, double p11, double p01, double p00,
                 double r10, double r11, double r01, double r00, int u, int v) :
    source(u), target(v) {cond_prob[0][0]=p00; cond_prob[0][1]=p01; cond_prob[1][0]=p10; cond_prob[1][1]=p11; factor[0][0]=r00; factor[0][1]=r01; factor[1][0]=r10; factor[1][1]=r11;};
    
//    double r10;         // store P(1|0)/P(0|0)
//    double r11;         //  and  P(1|1)/P(0|0)
//    double r01;
//    double r00;

    double cond_prob[2][2]; // P(1|0), P(1|1)
    double factor[2][2];
    int source;
    int target;
};

struct TreeNode {
    TreeNode() : parent(-1) {};
    int parent;
    vector<int> children;
};

class TreeBasin : public BasinModel
{
public:
    TreeBasin(int,int,RNG*);
    static vector<int> get_active_constraints(const State&);
    
    void doMLE(double);
    double P_state(const State&) const;
    vector<char> sample();
    
    paramsStruct get_params();
private:
//    vector<double> P0;
    double P0;
    double alpha;       // regularization parameter
    
    vector<TreeNode> adj_list;
    vector<TreeEdgeProb> edge_list;
    vector<int> below_thresh_list;
    myMatrix<double> J;
    myMatrix<double> m;
//    vector<int> roots;

    Graph G;
    
    double compute_MI(double,double,double);
};


#endif
