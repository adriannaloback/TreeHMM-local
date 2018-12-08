//--------------------------------------------
//  TreeBasin.cpp
//  
//
//  Created by Jason Prentice on 11/28/12.
//
//--------------------------------------------

#include "TreeBasin.h"
#include "EMBasins.h"

#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/foreach.hpp>

#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <queue>
#include <set>


TreeBasin::TreeBasin(int N, int basin_num, RNG* rng) : BasinModel(N,basin_num,rng), G(N), adj_list(N) {
    // Initialize randomly
    // stats 0 to N-1 are <sigma_i>
    int nstats = (N%2==0) ? (N/2)*(N+1) : N*((N+1)/2);
    stats.assign(nstats, 0);
    for (int i=0; i<N; i++) {
        double u = 0.1*((double) rand() / (double) RAND_MAX) + 0.45;
        stats[i] = u;
    }
    // stats N to N(N-1)/2 are <sigma_i sigma_j>; i<j
    for (int i=0; i<N; i++) {
        int ix = (i%2==0) ? (i/2)*(i-1) : i*((i-1)/2);
        for (int j=0; j<i; j++) {
//            double u = 0.6*((double) rand() / (double) RAND_MAX) + 0.2;
//            stats[N + ix + j] = u*min(stats[i],stats[j]);
            stats[N + ix + j] = stats[i]*stats[j];            
        }
    }
    
    // Build graph
    BOOST_FOREACH (Vertex v, vertices(G))
    {
        BOOST_FOREACH (Vertex u, vertices(G))
        {
            if (u != v) {
                EdgeProperty ep;
                ep.negI = 0;                
                add_edge(u, v, ep, G);
            }
        }
    }
    

    //doMLE();
}


vector<int> TreeBasin::get_active_constraints(const State& this_state) {
    vector<int> cons;
    // Constraints 0 to N-1 are the sigma_i
    for (vector<int>::const_iterator it=this_state.on_neurons.begin(); it!=this_state.on_neurons.end(); ++it) {
        cons.push_back(*it);
    }
    
    // Constraints N to N(N-1)/2 are (sigma_i sigma_j); i<j
    for (vector<int>::const_iterator it1=this_state.on_neurons.begin(); it1!=this_state.on_neurons.end(); ++it1) {
        for (vector<int>::const_iterator it2=this_state.on_neurons.begin(); it2!=it1; ++it2) {
            int i = max(*it1, *it2);
            int j = min(*it1, *it2);
            int N = this_state.word.size();
            int ix = (i%2==0) ? (i/2)*(i-1) : i*((i-1)/2);
            cons.push_back(N + ix + j);
        }
    }
    return cons;
    
}

void TreeBasin::doMLE(double alpha) {
    // Add MI to edges of G
    BOOST_FOREACH (Edge e, edges(G)) {
        Vertex u = source(e, G);
        Vertex v = target(e, G);
        int i = max((int) u, (int) v);
        int j = min((int) u, (int) v);
        int ix = (i%2==0) ? (i/2)*(i-1) : i*((i-1)/2);
        double C = stats[N + ix + j];
        if (C > stats[i]*stats[j] + alpha) {
            G[e].negI = -compute_MI(C-alpha, stats[i], stats[j]);
        } else if (C < stats[i]*stats[j] - alpha) {
            G[e].negI = -compute_MI(C+alpha, stats[i], stats[j]);
        } else {
            G[e].negI = 0;
        }
        
    }

    // Find spanning tree (maximizes likelihood over tree topology)
    vector<Edge> tree;
    kruskal_minimum_spanning_tree(G,
                                  back_inserter(tree),
                                  boost::weight_map(get(&EdgeProperty::negI,G)));

    
   
    vector<vector<int> > aux_adj_list (N);
    BOOST_FOREACH(Edge e, tree)
    {
        Vertex u = source(e, G);
        Vertex v = target(e, G);
        aux_adj_list[u].push_back(v);
        aux_adj_list[v].push_back(u);
    }
    
    //cout << "Getting probabilities" << endl;
    // Associate appropriate conditional probabilities with each edge. Need a directionality on the tree, done by BFS.  Vertex 0 will be the root.
    double thresh = 0.1;
//    P0.clear();
    edge_list.clear();
    below_thresh_list.clear();
//    roots.clear();
    for (int i=0; i<N; i++) {
        adj_list[i].parent = -1;
        adj_list[i].children.clear();
    }
    
    vector<char> visited (N,0);

    //cout << "Roots: ";
//    for (int this_root=0; this_root<N; this_root++) {
  //      if (visited[this_root] == 0) {
            //cout << this_root << ", ";
//            roots.push_back(this_root);
//            P0.push_back(1);
            P0 = 1;

            queue<int> to_process;
            to_process.push(0);

            while (!to_process.empty()) {
                int curr_node = to_process.front();
                //cout << curr_node << endl;
                visited[curr_node] = 1;

                for (vector<int>::iterator it=aux_adj_list[curr_node].begin(); it!=aux_adj_list[curr_node].end(); ++it)
                {
                    if (visited[*it]==0) {
                        to_process.push(*it);

                        double m1 = stats[curr_node];
                        double m2 = stats[*it];
                        int i = max(*it, curr_node);
                        int j = min(*it, curr_node);
                        int ix = (i%2==0) ? (i/2)*(i-1) : i*((i-1)/2);
                        double C = stats[N + ix + j];
                        if (C > m1*m2 + alpha) {
                            C -= alpha;
                        } else if (C < m1*m2 - alpha) {
                            C += alpha;
                        } else {
                            C = m1*m2;
                        }
                        if (m1==0 || m2==0) { C = 0; };
                        
                        double p11 = (m1>0) ? C/m1 : 0;
                        double p01 = 1 - p11;
                        double p10 = (m1<1) ? (m2-C)/(1-m1) : 0;
                        double p00 = 1 - p10;
                        
                        double r11=p11;
                        double r01=p01;
                        double r10=p10;
                        double r00=p00;
                        
                        bool below_thresh = false;
                        if (p00 > thresh) {
                            r11 /= p00;
                            r01 /= p00;
                            r10 /= p00;
                            P0 *= p00;
                            r00 = 1;
                        } else {
                            below_thresh = true;
                        }
                        
                        edge_list.push_back(TreeEdgeProb(p10,p11,p01,p00, r10,r11,r01,r00, curr_node, *it));
                        int new_edge = edge_list.size() - 1;
                        if (below_thresh) below_thresh_list.push_back(new_edge);
                        
                        adj_list[curr_node].children.push_back(new_edge);
                        adj_list[*it].parent = new_edge;

                    }
                }
                
                to_process.pop();
            }
//        }
        
//    }
    //cout << endl;
    return;
}

double TreeBasin::P_state(const State& this_state) const {


    if (edge_list.empty()) {
        double P = 1;
        for (int i=0; i<N; i++) {
            P *= (this_state.word[i]==0) ? 1-stats[i] : stats[i];
        }
        return P;
    }

    double P = P0;
    // Factor due to root neuron
//    for (int root_ix=0; root_ix < roots.size(); root_ix++) {
//        P *= P0[root_ix];
        P *= (this_state.word[0]==0) ? 1-stats[0] : stats[0];
//    }
    
    // Factor due to all edges with at least one spike.
    for (vector<int>::const_iterator it=this_state.on_neurons.begin(); it!=this_state.on_neurons.end(); ++it) {
        if (adj_list[*it].parent > -1) {  
            TreeEdgeProb parent = edge_list[adj_list[*it].parent];
            int sigma = this_state.word[parent.source];
            P *= parent.factor[1][sigma];
        }
        for (vector<int>::const_iterator e = adj_list[*it].children.begin(); e!=adj_list[*it].children.end(); ++e) {
            TreeEdgeProb child = edge_list[*e];
            int sigma = this_state.word[child.target];
            if (sigma==0) {          // Prevent double-counting 11 edges
                P *= child.factor[0][1];
            }
        }
    }

    // Factor due to below-threshold edges with no spikes
    for (vector<int>::const_iterator e=below_thresh_list.begin(); e!=below_thresh_list.end(); ++e)
    {
        TreeEdgeProb edge = edge_list[*e];
        if (this_state.word[edge.source] == 0 && this_state.word[edge.target]==0) {
            P *= edge.factor[0][0];
        }
    }
    
    
    return P;
}

double TreeBasin::compute_MI(double Cij, double pi, double pj) {
    double P_joint[4] = {Cij, (pj - Cij), (pi - Cij), (1 - pi - pj + Cij)};
    double S_joint = 0;
    for (int i=0; i<4; i++) {
        S_joint += (P_joint[i]>DBL_EPSILON ? -P_joint[i]*log(P_joint[i]) : 0);
    }
    double Si = (pi>DBL_EPSILON ? -pi*log(pi) : 0) + (pi<(1-DBL_EPSILON) ? (pi-1)*log(1-pi) : 0);
    double Sj = (pj>DBL_EPSILON ? -pj*log(pj) : 0) + (pj<(1-DBL_EPSILON) ? (pj-1)*log(1-pj) : 0);
    
    double I = Si + Sj - S_joint;
    
    return I;
}


vector<char> TreeBasin::sample() {
    vector<char> this_sample (N);
    queue<int> to_process;
    to_process.push(0);
    this_sample[0] = (rng->bernoulli(stats[0])) ? 1 : 0;
    while (!to_process.empty()) {
        int curr_node = to_process.front();
        char sigma = this_sample[curr_node];
        for (vector<int>::const_iterator it=adj_list[curr_node].children.begin(); it!=adj_list[curr_node].children.end(); ++it) {
            int next_node = edge_list[*it].target;
            double p = edge_list[*it].cond_prob[1][sigma];
            this_sample[next_node] = (rng->bernoulli(p)) ? 1 : 0;
            to_process.push(next_node);
        }
        to_process.pop();
    }

    return this_sample;
}

paramsStruct TreeBasin::get_params() {
    vector<double> vec_m (N);
    for (int i=0; i<N; i++) {
        vec_m[i] = stats[i];
    }
    m.assign(vec_m, N,1);
    
    vector<double> vec_J (N*N, 0);
    for (vector<TreeEdgeProb>::iterator it = edge_list.begin(); it!=edge_list.end(); ++it) {
        if (it->cond_prob[1][1] > 0 && it->cond_prob[0][0] > 0 && it->cond_prob[1][0] > 0 && it->cond_prob[0][1] > 0) {
            double this_J = log((it->cond_prob[1][1] * it->cond_prob[0][0]) / (it->cond_prob[1][0] * it->cond_prob[0][1]));
            int i = it->source;
            int j = it->target;
            vec_J[i*N+j] = this_J;
            vec_J[j*N+i] = this_J;
        }
    }
    J.assign(vec_J,N,N);
    
    paramsStruct params;
    params.addField("m",m);
    params.addField("J",J);
    
    return params;
    
}

