#include "graph.hpp"

Graph::Graph() {
    log_f_ = nullptr;
    nbr_cnt_ = nullptr;
    edges_t_ = nullptr;
    cn_ = nullptr;
    cn_size_ = nullptr;
    sub_nbr_cnt_ = nullptr;
    cn_t_ = nullptr;
    cn_cnt_ = nullptr;
    cn_t_idx_ = nullptr;
    v_a_ = nullptr;
    v_b_ = nullptr;
    idx_size_ = 0;
}

Graph::~Graph() {
    if (cn_ != nullptr) {
        delete[] cn_;
        cn_ = nullptr;
    }
    if (nbr_cnt_ != nullptr) {
        delete[] nbr_cnt_;
        nbr_cnt_ = nullptr;
    }
    if (sub_nbr_cnt_ != nullptr) {
        delete[] sub_nbr_cnt_;
        sub_nbr_cnt_ = nullptr;
    }
    if (edges_t_ != nullptr) {
        delete[] edges_t_;
        edges_t_ = nullptr;
    }
    if (cn_size_ != nullptr) {
        delete[] cn_size_;
        cn_size_ = nullptr;
    }
    if (cn_t_ != nullptr) {
        delete[] cn_t_;
        cn_t_ = nullptr;
    }
    if (cn_cnt_ != nullptr) {
        delete[] cn_cnt_;
        cn_cnt_ = nullptr;
    }
    if (cn_t_idx_ != nullptr) {
        delete[] cn_t_idx_;
        cn_t_idx_ = nullptr;
    }
    if (v_a_ != nullptr) {
        delete[] v_a_;
        v_a_ = nullptr;
    }
    if (v_b_ != nullptr) {
        delete[] v_b_;
        v_b_ = nullptr;
    }
    if (log_f_ != nullptr) {
        fclose(log_f_);
        log_f_ = nullptr;
    }
}

void Graph::load(const string& path) {
    if (log_f_ != nullptr) fprintf(log_f_, "Graph path: %s\n", path.c_str());
    printf("Graph path: %s\n", path.c_str());
    // printf("Loading Graph\n");

    ifstream ifs(path);
    if (!ifs.is_open()) {
        cerr << "open file failed!" << endl;
        exit(-1);
    }

    n_ = 0;
    m_ = 0;
    t_ = 0;

    int u, v;
    long ts, pre_ts = -1;
    while (!ifs.eof()) {
        ifs >> u >> v >> ts;
        if (ifs.fail()) break;  // 防止多读一个空行

        if (u == v) continue;

        if (u > v) swap(u, v);  // 默认 u < v
        edges_.emplace_back(make_pair(u, v));

        // adjust size of neighbor list if necessary.
        if (v + 1 > nbr_.size()) {
            nbr_.resize(v + 1);
        }

        // pre_ts = -1;
        // if (!t_new_to_old_.empty()) {
        //     pre_ts = t_new_to_old_.back();
        // }

        if (ts != pre_ts) {
            // 一起添加的 两个下标一致
            t_new_to_old_.emplace_back(ts);
            edges_idx_.emplace_back(edges_.size() - 1);
            pre_ts = ts;
        }

        int format_t = t_new_to_old_.size() - 1;

        nbr_[u].emplace_back(make_pair(v, format_t));
        // if(nbr_[u].size() > max_deg_) max_deg_ = nbr_[u].size();
        nbr_[v].emplace_back(make_pair(u, format_t));
        // if(nbr_[v].size() > max_deg_) max_deg_ = nbr_[v].size();

        ++m_;  // 统计使用的时间边，不包括u=v
    }
    ifs.close();

    n_ = nbr_.size();
    t_ = t_new_to_old_.size();
    edges_idx_.emplace_back(edges_.size());

    init_nbr_cnt();
    // init_edges_time();
    printf("n = %d, m = %d, t = %d, effective_m_ = %d.\n", n_, m_, t_, effective_m_);

    if (log_f_ != nullptr) {
        fprintf(log_f_, "n = %d, m = %d, t = %d, effective_m_ = %d.\n", n_, m_, t_, effective_m_);
    }

    print_graph_size();
}

void Graph::init_nbr_cnt() {
    effective_m_ = 0;  // 不重复边
    nbr_cnt_ = new unordered_map<int, int>[n_];
    edges_t_ = new unordered_map<int, vector<int>>[n_];
    for (int u = 0; u < n_; ++u) {
        for (auto& i : nbr_[u]) {
            if (u > i.first) continue;  // 统计 u < v

            if (nbr_cnt_[u].find(i.first) != nbr_cnt_[u].end()) {
                ++nbr_cnt_[u][i.first];
                edges_t_[u][i.first].emplace_back(i.second);
            } else {
                nbr_cnt_[u].emplace(i.first, 1);
                edges_t_[u].emplace(i.first, vector<int>{i.second});
            }
        }
        effective_m_ += nbr_cnt_[u].size();
    }
    // effective_m_ /= 2;
}

void Graph::index_baseline() {
#ifdef _LINUX_
    struct timeval t_start, t_end;
    gettimeofday(&t_start, NULL);
#else
    clock_t start = clock();
#endif

    if (cn_t_idx_ == nullptr)
        cn_t_idx_ = new unordered_map<int, vector<vector<pair<int, int>>>>[n_];
    else {
        printf("Index structure exists!\n");
        exit(1);
    }

    if (sub_nbr_cnt_ == nullptr) sub_nbr_cnt_ = new unordered_map<int, int>[n_];
    if (cn_size_ == nullptr) cn_size_ = new unordered_map<int, int>[n_];

    printf("Getting common neighbors...\n");

    get_common_neighbors_bl();  // 计算所有时间的cn

    for (int u = 0; u < n_; ++u) {
        for (auto& i : cn_[u]) {
            // 为 k=0,1,2,3..初始化
            cn_t_idx_[u][i.first].resize((i.second).size());
        }
    }

    for (int t_s = 0; t_s < t_; ++t_s) {
        if (t_s % 100 == 0) printf("t = %d.\n", t_s);  // 记录运行程度

        for (int u = 0; u < n_; ++u) {
            sub_nbr_cnt_[u] = nbr_cnt_[u];
            for (auto& i : cn_[u]) {
                cn_size_[u][i.first] = i.second.size();
            }
        }

        compute_cn_time_bl(t_s);  // 对开始时间为t_s 计算cn_time
        if (t_s == t_ - 1) break;
        decremental_cn_bl(t_s);
    }

#ifdef _LINUX_
    gettimeofday(&t_end, NULL);
    long long t_msec = (t_end.tv_sec - t_start.tv_sec) * 1000 + (t_end.tv_usec - t_start.tv_usec) / 1000;
    printf("Running time (baseline): %lld s, %lld mins\n", t_msec / 1000, t_msec / 1000 / 60);
    if (log_f_ != nullptr) fprintf(log_f_, "Indexing time (Baseline): %lld s, %lld mins\n", t_msec / 1000, t_msec / 1000 / 60);
#else
    clock_t end = clock();
    printf("Running time (baseline): %.2f s, %.2f min\n", (double)(end - start) / CLOCKS_PER_SEC, (double)(end - start) / CLOCKS_PER_SEC / 60);
    if (log_f_ != nullptr) fprintf(log_f_, "Running time (baseline): %.2f s, %.2f min\n", (double)(end - start) / CLOCKS_PER_SEC, (double)(end - start) / CLOCKS_PER_SEC / 60);
#endif

    print_idx_size();
}

// 计算全部时间子图的共同邻居，不带时间边
void Graph::get_common_neighbors_bl() {
    if (cn_ == nullptr) cn_ = new unordered_map<int, vector<int>>[n_];
    if (v_a_ == nullptr) v_a_ = new bool[n_];
    if (v_b_ == nullptr) v_b_ = new bool[n_];
    memset(v_a_, false, sizeof(bool) * n_);
    memset(v_b_, false, sizeof(bool) * n_);

    for (int u = 0; u < n_; ++u) {
        // v_a_  v_b_ 是u的访问标记
        for (auto& i : nbr_[u]) v_a_[i.first] = true;

        for (auto& i : nbr_[u]) {
            int v = i.first;

            if (v < u || v_b_[v]) continue;  // 若nbr[u]中v已经访问过，跳过

            for (auto& item : nbr_[v]) {
                int w = item.first;
                if (v_a_[w] && find(cn_[u][v].begin(), cn_[u][v].end(), w) == cn_[u][v].end()) {
                    cn_[u][v].emplace_back(w);
                }
                v_b_[v] = true;
            }
        }
        memset(v_a_, false, sizeof(bool) * n_);
        memset(v_b_, false, sizeof(bool) * n_);
    }
}

void Graph::compute_cn_time_bl(const int t_s) {
    for (int t_e = t_ - 1; t_e >= t_s; --t_e) {  // 从t-max 减到 t_s
        for (int i = edges_idx_[t_e]; i < edges_idx_[t_e + 1]; ++i) {
            int u = edges_[i].first;
            int v = edges_[i].second;

            if (--sub_nbr_cnt_[u][v] != 0) continue;

            // process (u,v)
            for (int k = 0; k < cn_size_[u][v]; ++k) {
                if (t_s == 0 || cn_t_idx_[u][v][k].back().second < t_e) {  //|| cn_t_idx_[u][v][k].empty()
                    cn_t_idx_[u][v][k].emplace_back(make_pair(t_s, t_e));  // cn_size_[u][v]=0
                }
                int x, y, cs;

                // process (u,w)
                x = u, y = cn_[u][v][k];
                if (x > y) swap(x, y);
                cs = cn_size_[x][y];
                if (t_s == 0 || cn_t_idx_[x][y][cs - 1].back().second < t_e) {  // || empty()
                    cn_t_idx_[x][y][cs - 1].emplace_back(make_pair(t_s, t_e));
                }
                auto it_v = find(cn_[x][y].begin(), cn_[x][y].begin() + cs, v);
                iter_swap(it_v, cn_[x][y].begin() + cs - 1);
                --cn_size_[x][y];

                // process (v,w)
                x = v, y = cn_[u][v][k];
                if (x > y) swap(x, y);
                cs = cn_size_[x][y];
                if (t_s == 0 || cn_t_idx_[x][y][cs - 1].back().second < t_e) {  // || empty()
                    cn_t_idx_[x][y][cs - 1].emplace_back(make_pair(t_s, t_e));
                }
                auto it_u = find(cn_[x][y].begin(), cn_[x][y].begin() + cs, u);
                iter_swap(it_u, cn_[x][y].begin() + cs - 1);
                --cn_size_[x][y];
            }
        }
    }
}

//  删除 t 时间的边
void Graph::decremental_cn_bl(const int t_s) {
    for (int i = edges_idx_[t_s]; i < edges_idx_[t_s + 1]; ++i) {
        int u = edges_[i].first;
        int v = edges_[i].second;

        if (--nbr_cnt_[u][v] != 0) continue;

        // process (u,v)
        for (int k = 0; k < cn_[u][v].size(); ++k) {
            cn_t_idx_[u][v][k].emplace_back(make_pair(t_s + 1, t_));

            int x, y;
            // process (u,w)
            x = u, y = cn_[u][v][k];
            if (x > y) swap(x, y);
            cn_t_idx_[x][y][cn_[x][y].size() - 1].emplace_back(make_pair(t_s + 1, t_));
            auto it_v = find(cn_[x][y].begin(), cn_[x][y].end(), v);
            iter_swap(it_v, cn_[x][y].end() - 1);
            cn_[x][y].pop_back();

            // process (v,w)
            x = v, y = cn_[u][v][k];
            if (x > y) swap(x, y);
            cn_t_idx_[x][y][cn_[x][y].size() - 1].emplace_back(make_pair(t_s + 1, t_));
            auto it_u = find(cn_[x][y].begin(), cn_[x][y].end(), u);
            iter_swap(it_u, cn_[x][y].end() - 1);
            cn_[x][y].pop_back();
        }
        cn_[u][v].clear();
        cn_[u][v].shrink_to_fit();
    }
}

// void Graph::write_idx_txt(const string& path) {
//     ofstream fout(path, ios::out);  // out 文件清空重写
//     for (int u = 0; u < n_; ++u) {
//         for (auto& i : cn_t_idx_[u]) {
//             fout << "\n --------(" << u << "," << i.first << ")-------" << endl;
//             for (int k = 0; k < i.second.size(); ++k) {
//                 fout << "cn = " << k + 1 << endl;
//                 for (auto& t : i.second[k]) {
//                     fout << "[" << t.first + 1 << "," << t.second + 1 << "]" << " ";
//                 }
//                 fout<<endl;
//             }
//         }
//     }
// }

void Graph::write_idx(const string& path) {
    auto fp = fopen(path.c_str(), "wb");
    fwrite(&n_, sizeof(unsigned int), 1, fp);  // u

    for (int u = 0; u < n_; ++u) {
        int u_size = cn_t_idx_[u].size();
        fwrite(&u_size, sizeof(int), 1, fp);  // u.size
        for (auto& nbr : cn_t_idx_[u]) {
            int v = nbr.first;
            int cn_size = nbr.second.size();
            fwrite(&v, sizeof(int), 1, fp);        // v
            fwrite(&cn_size, sizeof(int), 1, fp);  // cn(u,v).size
            for (int k = 0; k < cn_size; ++k) {
                int cnt = cn_t_idx_[u][v][k].size();
                fwrite(&cnt, sizeof(int), 1, fp);
                // for (auto& i : cn_t_idx_[u][v][k]) {
                //     fwrite(&i.first, sizeof(int), 1, fp);
                //     fwrite(&i.second, sizeof(int), 1, fp);
                // }
                fwrite(&cn_t_idx_[u][v][k][0], sizeof(pair<int, int>) * cnt, 1, fp);
            }
        }
    }
    fclose(fp);
    printf("Write index.\n");
}

void Graph::load_idx(const string& path) {
    auto fp = fopen(path.c_str(), "rb");
    if (fp == NULL) {
        printf("cannot open file!\n");
        exit(1);
    }
    printf("Loading index.\n");

    fread(&n_, sizeof(unsigned int), 1, fp);
    t_ = 0;

    if (cn_t_idx_ == nullptr) cn_t_idx_ = new unordered_map<int, vector<vector<pair<int, int>>>>[n_];

    for (int u = 0; u < n_; ++u) {
        int u_size;
        fread(&u_size, sizeof(int), 1, fp);
        // cn_t_idx_[u].reserve(u_size);
        for (int nbr = 0; nbr < u_size; ++nbr) {
            int v, cn_size;
            fread(&v, sizeof(int), 1, fp);
            fread(&cn_size, sizeof(int), 1, fp);
            cn_t_idx_[u][v].reserve(cn_size);
            for (int k = 0; k < cn_size; ++k) {
                int cnt;
                fread(&cnt, sizeof(int), 1, fp);
                cn_t_idx_[u][v][k].resize(cnt);
                fread(&cn_t_idx_[u][v][k][0], cnt * sizeof(pair<int, int>), 1, fp);
                // cn_t_idx_[u][v][k].reserve(cnt);
                // for (int i = 0; i < cnt; ++i) {
                //     int a, b;
                //     fread(&a, sizeof(int), 1, fp);
                //     fread(&b, sizeof(int), 1, fp);
                //     cn_t_idx_[u][v][k].emplace_back(make_pair(a, b));
                // }
                if (t_ < cn_t_idx_[u][v][k].back().second) t_ = cn_t_idx_[u][v][k].back().second;
            }
        }
    }
    fclose(fp);
}

void Graph::index() {
#ifdef _LINUX_
    struct timeval t_start, t_end;
    gettimeofday(&t_start, NULL);
#else
    clock_t start = clock();
#endif

    if (cn_t_idx_ == nullptr) cn_t_idx_ = new unordered_map<int, vector<vector<pair<int, int>>>>[n_];
    if (cn_cnt_ == nullptr) cn_cnt_ = new unordered_map<int, vector<int>>[n_];

    printf("Getting common neighbors with minimum time...\n");
    get_common_neighbors();

    for (int u = 0; u < n_; ++u) {
        for (auto& i : cn_t_[u]) {
            // 为 cn= 0,1,2...初始化，
            cn_t_idx_[u][i.first].resize((i.second).size());
            cn_cnt_[u][i.first].resize((i.second).size());
        }
    }

    printf("Initialize common neighbor time.\n");
    init_common_neighbor_time();

    queue<tuple<int, int, int>> q;
    int* ptr_cn_cnt = nullptr;
    vector<pair<int, int>>* ptr_cn_t = nullptr;

    for (int t_s = 1; t_s < t_; ++t_s) {
        if (t_s % 1000 == 0) printf("t = %d.\n", t_s);                 // record time
        for (int e = edges_idx_[t_s - 1]; e < edges_idx_[t_s]; ++e) {  // delete t-1
            int u = edges_[e].first;
            int v = edges_[e].second;

            if (--nbr_cnt_[u][v] == 0) {  // delete completely

                // process (u,v)
                for (int k = 0; k < cn_t_[u][v].size(); ++k) {
                    if (cn_cnt_[u][v][k] > 0) {
                        q.push(make_tuple(u, v, k));
                        // cn_cnt_[u][v][k] = -1;
                        cn_t_idx_[u][v][k].emplace_back(make_pair(t_s, t_));
                        cn_cnt_[u][v][k] = 0;
                    }
                    int w = cn_t_[u][v][k].first;
                    int old_t = cn_t_[u][v][k].second;

                    // process (u,w)
                    ptr_cn_t = u < w ? &cn_t_[u][w] : &cn_t_[w][u];
                    ptr_cn_cnt = u < w ? &cn_cnt_[u][w][0] : &cn_cnt_[w][u][0];
                    auto it_v_t = lower_bound(ptr_cn_t->begin(), ptr_cn_t->end(), make_pair(v, old_t), cmp2);  // find t = old_t
                    for (int k = it_v_t - ptr_cn_t->begin(); k < ptr_cn_t->size(); ++k) {
                        if (ptr_cn_cnt[k] > 0 && --ptr_cn_cnt[k] < k + 1) {
                            q.push(make_tuple(u, w, k));
                            ptr_cn_cnt[k] = -1;
                        }
                    }
                    while (it_v_t->first != v) ++it_v_t;
                    ptr_cn_t->erase(it_v_t);  // 直接得到时间窗口

                    // process (v,w)
                    ptr_cn_t = v < w ? &cn_t_[v][w] : &cn_t_[w][v];
                    ptr_cn_cnt = v < w ? &cn_cnt_[v][w][0] : &cn_cnt_[w][v][0];
                    auto it_u_t = lower_bound(ptr_cn_t->begin(), ptr_cn_t->end(), make_pair(u, old_t), cmp2);
                    for (int k = it_u_t - ptr_cn_t->begin(); k < ptr_cn_t->size(); ++k) {
                        if (ptr_cn_cnt[k] > 0 && --ptr_cn_cnt[k] < k + 1) {
                            q.push(make_tuple(v, w, k));
                            ptr_cn_cnt[k] = -1;
                        }
                    }
                    while (it_u_t->first != u) ++it_u_t;
                    ptr_cn_t->erase(it_u_t);
                }
                cn_t_[u][v].clear();
                cn_t_[u][v].shrink_to_fit();
            } else {
                int new_t = edges_t_[u][v][edges_t_[u][v].size() - nbr_cnt_[u][v]];
                int cnt = 0;

                for (auto& i : cn_t_[u][v]) {
                    // process (u,v)
                    if (new_t <= i.second) break;
                    int old_t = i.second;
                    i.second = new_t;

                    if (cn_cnt_[u][v][cnt] > 0) {
                        q.push(make_tuple(u, v, cnt));
                        cn_cnt_[u][v][cnt] = -1;
                    }
                    ++cnt;

                    int w = i.first, k;

                    // process(u,w)
                    ptr_cn_t = u < w ? &cn_t_[u][w] : &cn_t_[w][u];
                    ptr_cn_cnt = u < w ? &cn_cnt_[u][w][0] : &cn_cnt_[w][u][0];
                    auto it_v_t = lower_bound(ptr_cn_t->begin(), ptr_cn_t->end(), make_pair(v, old_t), cmp2);
                    for (k = it_v_t - ptr_cn_t->begin(); k < ptr_cn_t->size(); ++k) {
                        if ((*ptr_cn_t)[k].second >= new_t) break;
                        if (ptr_cn_cnt[k] > 0 && --ptr_cn_cnt[k] < k + 1) {
                            q.push(make_tuple(u, w, k));
                            ptr_cn_cnt[k] = -1;
                        }
                    }
                    while (it_v_t->first != v) ++it_v_t;
                    it_v_t->second = new_t;
                    rotate(it_v_t, it_v_t + 1, ptr_cn_t->begin() + k);

                    // process(v,w)
                    ptr_cn_t = v < w ? &cn_t_[v][w] : &cn_t_[w][v];
                    ptr_cn_cnt = v < w ? &cn_cnt_[v][w][0] : &cn_cnt_[w][v][0];
                    auto it_u_t = lower_bound(ptr_cn_t->begin(), ptr_cn_t->end(), make_pair(u, old_t), cmp2);
                    for (k = it_u_t - ptr_cn_t->begin(); k < ptr_cn_t->size(); ++k) {
                        if ((*ptr_cn_t)[k].second >= new_t) break;
                        if (ptr_cn_cnt[k] > 0 && --ptr_cn_cnt[k] < k + 1) {
                            q.push(make_tuple(v, w, k));
                            ptr_cn_cnt[k] = -1;
                        }
                    }
                    while (it_u_t->first != u) ++it_u_t;
                    it_u_t->second = new_t;
                    rotate(it_u_t, it_u_t + 1, ptr_cn_t->begin() + k);
                }
            }
        }
        int u, v, k;
        while (!q.empty()) {  // 考虑unordered_map, 增加cpu缓存命中率
            tie(u, v, k) = q.front();
            q.pop();
            if (u > v) swap(u, v);
            if (cn_t_[u][v].size() < k + 1) {
                cn_t_idx_[u][v][k].emplace_back(make_pair(t_s, t_));
                cn_cnt_[u][v][k] = 0;
            } else {
                int t = cn_t_[u][v][k].second;
                cn_t_idx_[u][v][k].emplace_back(make_pair(t_s, t));
                int i = k + 1;
                while (i < cn_t_[u][v].size() && cn_t_[u][v][i].second == t) ++i;
                cn_cnt_[u][v][k] = i;
            }
        }
    }
    ptr_cn_cnt = nullptr;
    ptr_cn_t = nullptr;

#ifdef _LINUX_
    gettimeofday(&t_end, NULL);
    long long t_msec = (t_end.tv_sec - t_start.tv_sec) * 1000 + (t_end.tv_usec - t_start.tv_usec) / 1000;
    printf("Running time: %lld s, %lld mins\n", t_msec / 1000, t_msec / 1000 / 60);
    if (log_f_ != nullptr) fprintf(log_f_, "Indexing time: %lld s, %lld mins\n", t_msec / 1000, t_msec / 1000 / 60);

#else
    clock_t end = clock();
    printf("Running time: %.2f s, %.2f min\n", (double)(end - start) / CLOCKS_PER_SEC, (double)(end - start) / CLOCKS_PER_SEC / 60);
    if (log_f_ != nullptr) fprintf(log_f_, "Running time: %.2f s, %.2f min\n", (double)(end - start) / CLOCKS_PER_SEC, (double)(end - start) / CLOCKS_PER_SEC / 60);
#endif
    print_idx_size();
}

void Graph::get_common_neighbors() {
    if (cn_t_ == nullptr) cn_t_ = new unordered_map<int, vector<pair<int, int>>>[n_];
    if (v_a_ == nullptr) v_a_ = new bool[n_];
    if (v_b_ == nullptr) v_b_ = new bool[n_];
    memset(v_a_, false, sizeof(bool) * n_);
    memset(v_b_, false, sizeof(bool) * n_);

    for (int u = 0; u < n_; ++u) {
        for (auto& i : nbr_[u]) v_a_[i.first] = true;

        for (auto& i : nbr_[u]) {
            int v = i.first;
            if (u > v || v_b_[v]) continue;  // u < v
            for (auto& j : nbr_[v]) {        // nbr_ (sort by time)
                int w = j.first;
                if (v_a_[w] && find_if(cn_t_[u][v].begin(), cn_t_[u][v].end(), [&w](pair<int, int>& p) { return p.first == w; }) == cn_t_[u][v].end()) {
                    int time_u_w = u < w ? edges_t_[u][w].front() : edges_t_[w][u].front();
                    cn_t_[u][v].emplace_back(make_pair(w, max(time_u_w, max(i.second, j.second))));  // max{ t(u,w),t(u,v),t(v,w) }
                }
            }
            v_b_[v] = true;
        }
        memset(v_a_, false, sizeof(bool) * n_);
        memset(v_b_, false, sizeof(bool) * n_);
    }
}

void Graph::init_common_neighbor_time() {
    for (int u = 0; u < n_; ++u) {
        for (auto& i : cn_t_[u]) {
            int v = i.first, pre_t = -1;
            sort(cn_t_[u][v].begin(), cn_t_[u][v].end(), cmp2);
            for (int k = cn_t_[u][v].size() - 1; k >= 0; --k) {
                int t = cn_t_[u][v][k].second;
                cn_t_idx_[u][v][k].emplace_back(make_pair(0, t));
                if (pre_t == t) {
                    cn_cnt_[u][v][k] = cn_cnt_[u][v][k + 1];
                } else {
                    cn_cnt_[u][v][k] = k + 1;
                    pre_t = t;
                }
            }
        }
    }
}

void Graph::print_graph_size() {
    printf("Graph size: %.2f MB.\n", (float)edges_.size() * 3 * sizeof(int) / 1024 / 1024);
    if (log_f_ != nullptr) fprintf(log_f_, "Graph size: %.2f MB.\n", (float)edges_.size() * 3 * sizeof(int) / 1024 / 1024);
}

void Graph::print_idx_size() {
    if (idx_size_ != 0) printf("Index size: %lld MB.", idx_size_ / 1024 / 1024);

    idx_size_ += sizeof(int);
    idx_size_ += sizeof(int) * n_;  // n_ 个指针

    double average_t = 0;  // 平均时间窗口数量
    long average_t_d = 0;
    int max_t = 0;  // 最多时间窗口数量

    for (int u = 0; u < n_; ++u) {
        idx_size_ += cn_t_idx_[u].size() * sizeof(int);  // key size
        for (const auto& i : cn_t_idx_[u]) {
            int v = i.first;
            for (int k = i.second.size() - 1; k >= 0; --k) {
                idx_size_ += (i.second)[k].size() * 2 * sizeof(int);
                average_t += (i.second)[k].size();
                if ((i.second)[k].size() > max_t) max_t = (i.second)[k].size();
            }
            average_t_d += i.second.size();
        }
    }
    printf("Index size: %.2f MB.\n", (float)idx_size_ / 1024 / 1024);
    printf("Average cn = %.2f, average T = %.2f, max T = %d.\n", double(average_t_d) / effective_m_, average_t / average_t_d, max_t);
    if (log_f_ != nullptr) fprintf(log_f_, "Index size: %.2f MB\n", (float)idx_size_ / 1024 / 1024);
    if (log_f_ != nullptr) fprintf(log_f_, "Average cn = %.2f, average T = %.2f, max T = %d.\n", double(average_t_d) / effective_m_, average_t / average_t_d, max_t);
}

void Graph::init_log(const string& log_path) {
    log_f_ = fopen(log_path.c_str(), "a");  // 追加
    fprintf(log_f_, "\n\n===============================\n");
    time_t now = time(0);
    fprintf(log_f_, "%s\n", ctime(&now));
}

bool Graph::query(int u, int v, int t_s, int t_e, int k) {
    if (u > v) swap(u, v);
    if (cn_t_idx_[u][v].size() < k) return false;
    auto it = upper_bound(cn_t_idx_[u][v][k - 1].begin(), cn_t_idx_[u][v][k - 1].end(), make_pair(t_s, t_e), cmp1);
    --it;
    return it->second <= t_e;
}

bool cmp1(const pair<int, int>& a, const pair<int, int>& b) {
    return a.first < b.first;
}

bool cmp2(const pair<int, int>& a, const pair<int, int>& b) {
    return a.second < b.second;
}

bool cmp3(const pair<int, int>& a, const pair<int, int>& b) {
    if (a.second == b.second) return a.first < b.first;
    return a.second < b.second;
}
