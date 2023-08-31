#ifndef tpSCAN_hpp
#define tpSCAN_hpp

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <bitset>

#define _LINUX_

#ifdef _LINUX_
#include <sys/time.h>
#endif

using namespace std;

const int NMAX = 5000000;

class Graph {
   private:
    int n, m, t_max, ts, te, cal_similar_cnt;

    int eps_a2, eps_b2, miu;  // eps_a2/eps_b2 = eps^2
    string eps;

    int* pstart;   // offset of neighbors of nodes   pstart[i]存放 0-i点的度数总和
    int* edges;    // adjacent ids of edges
    int* reverse;  // the position of reverse edge in edges
    int* min_cn;   // minimum common neighbor: -2 means not similar; -1 means similar; 0 means not sure; > 0 means the minimum common neighbor

    int* pa;
    int* rank;  // pa and rank are used for the disjoint-set data structure

    int* cid;  // cluster id

    int* degree;
    int* similar_degree;    // number of adjacent edges with similarity no less than epsilon
    int* effective_degree;  // number of adjacent edges not pruned by similarity

    vector<pair<int, int>> noncore_cluster;
    
    vector<long> t_new_to_old_;
    vector<int> edges_t_idx_;
    vector<pair<int, int>> edges_t_;
    vector<vector<pair<int, int>>> nbr_t_;
    unordered_map<int, vector<vector<pair<int, int>>>>* cn_t_idx_; //index

   public:
    Graph();
    ~Graph();

    void read_graph(const string& path);
    void tpSCAN(string eps_s, int _miu, int t_s, int t_e);
    void cluster_noncore_vertices(int eps_a2, int eps_b2, int mu);
    void load_idx(const string& path);
    void output(const string& path);

   private:
    void get_snapshot(const int t_s, const int t_e);
    int binary_search(const int* array, int b, int e, int val);  // array数组中是否存在val
    int similar_check_OP(int u, int idx, int eps_a, int eps_b);
    int check_common_neighbor(int u, int v, int c);
    int query(int u, int v, int ts, int te, int k);
    int compute_common_neighbor_lowerbound(int u, int v, int eps_a2, int eps_b2);
    void prune_and_cross_link(int eps_a2, int eps_b2, int miu, int* cores, int& cores_e);

    int find_root(int u);
    void my_union(int u, int v);

    void get_eps(const string eps_s);

    int my_check_cn(int u, int v);
    void my_check_index(int u, int v);
};

bool cmp(const pair<int, int>& a, const pair<int, int>& b);

#endif
