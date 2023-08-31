#ifndef Graph_hpp
#define Graph_hpp

#include <string.h>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <list>
#include <queue>
#include <unordered_map>
#include <vector>

#define _LINUX_

#ifdef _LINUX_
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
// #include <sys/resource.h>
#endif

using namespace std;

class Graph {
    unsigned int n_{};  // 最大顶点号 + 1
    unsigned int m_{};  // 边数
    unsigned int effective_m_{};
    long long idx_size_;

    FILE* log_f_;

    vector<long> t_new_to_old_;
    vector<int> edges_idx_;         // 每个时间下标  上的数对应在edges_ 中的下标
    vector<pair<int, int>> edges_;  // 边集 pair<u,v>, 利用上面两个数据得到 每个时间点的边
    vector<vector<pair<int, int>>> nbr_;

    unordered_map<int, int>* nbr_cnt_;  
    unordered_map<int, int>* sub_nbr_cnt_;
    unordered_map<int, vector<int>>* edges_t_;
    unordered_map<int, vector<int>>* cn_;               // common neighbor
    unordered_map<int, vector<pair<int, int>>>* cn_t_;  // common neighbor with time
    unordered_map<int, int>* cn_size_;
    unordered_map<int, vector<int>>* cn_cnt_;
    unordered_map<int, vector<vector<pair<int, int>>>>* cn_t_idx_;  // index

    bool* v_a_;
    bool* v_b_;

    void init_nbr_cnt();
    void init_common_neighbor_time();
    void get_common_neighbors_bl();
    void get_common_neighbors();
    void compute_cn_time_bl(const int t_s);
    void decremental_cn_bl(const int t_s);

    void print_graph_size();
    void print_idx_size();

   public:
    unsigned int t_{};
    Graph();
    ~Graph();
    void load(const string& path);
    void index();
    void index_baseline();

    // void write_idx_txt(const string& path);
    void write_idx(const string& path);
    void load_idx(const string& path);
    void init_log(const string& log_path);

    bool query(int u, int v, int ts, int te, int k);
};

bool cmp1(const pair<int, int>& a, const pair<int, int>& b);
bool cmp2(const pair<int, int>& a, const pair<int, int>& b);
bool cmp3(const pair<int, int>& a, const pair<int, int>& b);

#endif /* Graph_hpp */
