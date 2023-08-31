#include "tpSCAN.hpp"

Graph::Graph() {
    n = m = t_max = ts = te = cal_similar_cnt = 0;
    eps_a2 = eps_b2 = miu = 0;
    pstart = nullptr;
    edges = nullptr;
    reverse = nullptr;
    min_cn = nullptr;
    cid = nullptr;
    degree = nullptr;
    effective_degree = nullptr;
    similar_degree = nullptr;
    pa = nullptr;
    rank = nullptr;

    cn_t_idx_ = nullptr;
}

Graph::~Graph() {
    if (pstart != nullptr) {
        delete[] pstart;
        pstart = nullptr;
    }
    if (edges != nullptr) {
        delete[] edges;
        edges = nullptr;
    }
    if (reverse != nullptr) {
        delete[] reverse;
        reverse = nullptr;
    }
    if (min_cn != nullptr) {
        delete[] min_cn;
        min_cn = nullptr;
    }
    if (cid != nullptr) {
        delete[] cid;
        cid = nullptr;
    }
    if (degree != nullptr) {
        delete[] degree;
        degree = nullptr;
    }
    if (effective_degree != nullptr) {
        delete[] effective_degree;
        effective_degree = nullptr;
    }
    if (similar_degree != nullptr) {
        delete[] similar_degree;
        similar_degree = nullptr;
    }
    if (pa != nullptr) {
        delete[] pa;
        pa = nullptr;
    }
    if (rank != nullptr) {
        delete[] rank;
        rank = nullptr;
    }
    if (cn_t_idx_ != nullptr) {
        delete[] cn_t_idx_;
        cn_t_idx_ = nullptr;
    }
}

void Graph::read_graph(const string& path) {
    printf("Graph path: %s\n", path.c_str());
    // printf("Reading Graph\n");

    ifstream ifs(path);
    if (!ifs.is_open()) {
        cerr << "open file failed!" << endl;
        exit(-1);
    }

    int u, v;
    long t, pre_t = -1;
    while (ifs.good() && !ifs.eof()) {
        ifs >> u >> v >> t;
        if (ifs.fail()) break;

        if (u == v) continue;

        if (u > v) swap(u, v);  // 默认 u < v
        edges_t_.emplace_back(make_pair(u, v));

        // adjust size of neighbor list if necessary.
        if (v + 1 > nbr_t_.size()) {
            nbr_t_.resize(v + 1);
        }

        if (t != pre_t) {
            // 一起添加的 两个下标一致
            t_new_to_old_.emplace_back(t);
            edges_t_idx_.emplace_back(edges_t_.size() - 1);
            pre_t = t;
        }

        int format_t = t_new_to_old_.size() - 1;

        nbr_t_[u].emplace_back(make_pair(v, format_t));
        nbr_t_[v].emplace_back(make_pair(u, format_t));
    }
    ifs.close();

    t_max = t_new_to_old_.size() - 1;
    edges_t_idx_.emplace_back(edges_t_.size());

    // sort neighbor[u] by id
}

void Graph::get_snapshot(const int t_s, const int t_e) {
    if (t_s < 0 || t_s > t_e || t_e > t_max) {
        printf("??? Wrong interval");
        exit(1);
    }
    printf("Interval: [%d, %d].\nt_max:%d.\n", t_s, t_e, t_max);

    for (int i = edges_t_idx_[t_s]; i < edges_t_idx_[t_e + 1]; ++i) {
        if (edges_t_[i].second + 1 > n) n = edges_t_[i].second + 1;  // u < v
    }

    int m_max = 2 * (edges_t_idx_[t_e + 1] - edges_t_idx_[t_s]);
    if (edges == nullptr) edges = new int[m_max];  // max edges
    if (pstart == nullptr) pstart = new int[n + 1];
    if (degree == nullptr) degree = new int[n];

    bitset<NMAX> uset;
    vector<int> temp_arr(n);

    pstart[0] = 0;
    int effective_n = 0;
    for (int u = 0; u < n; ++u) {
        int deg = 0;
        uset.reset();
        for (auto& i : nbr_t_[u]) {
            if (i.second < t_s)
                continue;
            else if (i.second > t_e)
                break;
            else if (uset.test(i.first))
                continue;
            else {
                edges[pstart[u] + deg] = i.first;
            }
        }
        if (deg != 0) ++effective_n;
        pstart[u + 1] = pstart[u] + deg;
        m += deg;
        degree[u] = deg;
        ++degree[u];  // d[u] = {u}
    }

    if (reverse == nullptr) reverse = new int[m];
    if (min_cn == nullptr) min_cn = new int[m];
    memset(min_cn, 0, sizeof(int) * m);

    printf("effective_n: %d (%.2f%%),  n_max:%d,  ", effective_n, (double)(effective_n * 100) / n, n - 1);
    printf("m: %d ,  avg_d(u): %.2f\n", m / 2, (double)m / effective_n);
}

int Graph::binary_search(const int* array, int b, int e, int val) {
    --e;
    if (array[e] < val) return e + 1;
    while (b < e) {
        int mid = b + (e - b) / 2;
        if (array[mid] >= val)
            e = mid;
        else
            b = mid + 1;
    }
    return e;
}

void Graph::tpSCAN(string eps_s, int _miu, int t_s, int t_e) {
    eps = eps_s;
    miu = _miu;
    ts = t_s;
    if (t_e == 0)
        te = t_max;
    else
        te = t_e;
    get_eps(eps_s);

#ifdef _LINUX_
    struct timeval snapshot;
    gettimeofday(&snapshot, NULL);
#endif
    get_snapshot(ts, te);

    if (similar_degree == nullptr) similar_degree = new int[n];  // 相似度初始化为0
    memset(similar_degree, 0, sizeof(int) * n);

    if (effective_degree == nullptr) effective_degree = new int[n];
    for (int i = 0; i < n; i++) effective_degree[i] = degree[i] - 1;  // 有效度初始化为 d[i]-1

    if (pa == nullptr) pa = new int[n];
    if (rank == nullptr) rank = new int[n];
    for (int i = 0; i < n; i++) {
        pa[i] = i;
        rank[i] = 0;
    }

#ifdef _LINUX_
    struct timeval start;
    gettimeofday(&start, NULL);
#else
    clock_t start = clock();
#endif

    int* edge_buf = new int[n];
    int* cores = new int[n];  // 记录核心顶点
    int cores_n = 0;          // 记录核心顶点的数量

    prune_and_cross_link(eps_a2, eps_b2, miu, cores, cores_n);
    // printf("\t*** Finished prune and cross link!\n");

#ifdef _LINUX_
    struct timeval end1;
    gettimeofday(&end1, NULL);
    long long time_use1 = (end1.tv_sec - start.tv_sec) * 1000000 + (end1.tv_usec - start.tv_usec);
#else
    clock_t end1 = clock();
#endif

    int* bin_head = new int[n];                    // 记录第一个ed值为n的点
    int* bin_next = new int[n];                    // bin_next的每一个值指向下一个ed值相同的点，遇到-1结束
    for (int i = 0; i < n; i++) bin_head[i] = -1;  // 初始值为-1

    int max_ed = 0;
    for (int i = 0; i < n; i++)
        if (effective_degree[i] >= miu) {  // 根据每个点的ed大小进行排序，从大到小开始探索
            int ed = effective_degree[i];
            if (ed > max_ed) max_ed = ed;  // 出现新的更大的ed值，标记一下
            bin_next[i] = bin_head[ed];    // 类似链表  bin_head[x] 记录ed=x 的最后一个顶点i（链表），满足ed[i]=x
            bin_head[ed] = i;
        }

    while (true) {  // 退出条件：ed大于miu的点都处理完了
        // 桶排序 按ed非递减顺序访问顶点
        int u = -1;
        if (cores_n)               // 优先处理 已经 确认是core顶点队列，聚类core-core
            u = cores[--cores_n];  // 先处理一遍core顶点(所需cn<=2，直接得到core)
        else {                     // 根据之前的排序找到ed最大的点
            while (max_ed >= miu && u == -1) {
                for (int x = bin_head[max_ed]; x != -1;) {  //-1表示没有这个ed的点
                    int tmp = bin_next[x];
                    int ed = effective_degree[x];  // ed 是真实有效度
                    if (ed == max_ed) {
                        u = x;                           // 找到当前ed最大的结点 退出循环
                        bin_head[max_ed] = bin_next[x];  // head指向倒数第二个点
                        break;
                    } else if (ed >= miu) {
                        bin_next[x] = bin_head[ed];  // 如果x有希望成为core，但ed(x)边小了，重新指向对应head
                        bin_head[ed] = x;
                    }
                    x = tmp;  // 对于真实ed(x)<miu 不是core 从链表中删除结点，同时x指向下一个初始ed为max_ed的结点
                }
                if (u == -1) {  // 初始ed为max_ed的结点 都找完了
                    bin_head[max_ed] = -1;
                    --max_ed;  // 找 --max_ed的链表
                }
            }
        }

        if (u == -1) break;  // ed大于miu的点都处理完了,退出最外层while

        // 否则，遍历 u 的邻居，观察点u是否为core
        int edge_buf_n = 0;
        for (int j = pstart[u]; j < pstart[u + 1]; j++) {  // j表示边 edges[j]表示邻点
            if (min_cn[j] == -2) continue;                 //-2 不相似, j 表示边
            // u，v可能相似
            //  顶点u非core 或者 u是core但uv没抱团，  退出if：u是core，且uv已抱团
            if (similar_degree[u] < miu || find_root(u) != find_root(edges[j])) edge_buf[edge_buf_n++] = j;  // 可能相似
        }

        // 算法4：checkcore!!! 检测 u 是否为 core

        int i = 0;                                                                         // effective_degree[u] >= miu  点u之前被作为邻点处理，可能有ed < miu，从而跳过点 u
        while (similar_degree[u] < miu && effective_degree[u] >= miu && i < edge_buf_n) {  // 还有机会相似的点
            int idx = edge_buf[i];                                                         // idx 是边的索引id
            if (min_cn[idx] != -1) {
                int v = edges[idx];
                // 相反边 min_cn 一样的，这里计算得到 -1或-2
                min_cn[idx] = min_cn[reverse[idx]] = similar_check_OP(u, idx, eps_a2, eps_b2);
                // 若相似 处理u
                if (min_cn[idx] == -1)
                    ++similar_degree[u];
                else
                    --effective_degree[u];

                if (effective_degree[v] >= 0) {  // 若 v还没被探索过  否则ed < 0 表示已被确认是否未core
                    if (min_cn[idx] == -1) {
                        ++similar_degree[v];

                        if (similar_degree[v] == miu) cores[cores_n++] = v;  // 邻点v未被探索，且是core，加入核心序列，优先处理
                    } else
                        --effective_degree[v];
                }
            }

            ++i;  // 退出循环条件，u已经被确认是core还是noncore  若sd[u]>=miu提前退出 可能还有剩余邻点没探索过
        }

        effective_degree[u] = -1;  // ed > 0 表示被访问过了！！！ //标记为已检查过的点

        if (similar_degree[u] < miu) continue;  // u非core 跳过这次大循环

        // 算法5：cluserCore  u是core，聚类邻点也是core的点，core-core

        for (int j = 0; j < edge_buf_n; j++) {
            int idx = edge_buf[j];  // idx是邻边id
            // 若 uv已经相似  并且  邻点是core
            if (min_cn[idx] == -1 && similar_degree[edges[idx]] >= miu)
                my_union(u, edges[idx]);  // 核心抱团
        }
        // 检查剩余邻点 但 相似度未计算
        while (i < edge_buf_n) {
            int idx = edge_buf[i];
            int v = edges[idx];
            // 有一个为true就跳过，1.以确定uv是否相似 2.不确定uv相似但 邻点v非core，这里只聚类core-core  3.v是core且已经抱团，无需计算相似性
            if (min_cn[idx] < 0 || similar_degree[v] < miu || find_root(u) == find_root(v)) {
                ++i;
                continue;
            }
            // 剩下情况是： (u是core)v是core，且uv没抱团，计算相似性
            min_cn[idx] = min_cn[reverse[idx]] = similar_check_OP(u, idx, eps_a2, eps_b2);  // 计算相似情况

            // 说明点v还没有检查过
            if (effective_degree[v] >= 0) {
                if (min_cn[idx] == -1) {
                    ++similar_degree[v];

                    if (similar_degree[v] == miu) cores[cores_n++] = v;  // 邻点v未被探索，且是core，加入核心序列，优先处理
                } else
                    --effective_degree[v];
            }

            if (min_cn[idx] == -1) my_union(u, v);  // 若相似 两个簇合并

            ++i;  // 直到所有邻点都被访问过
        }
        // printf(")\n");
    }
    // printf("\t*** Finished clustering core vertices!\n");

    delete[] edge_buf;
    edge_buf = nullptr;
    delete[] cores;
    cores = nullptr;
    delete[] bin_head;
    bin_head = nullptr;
    delete[] bin_next;
    bin_next = nullptr;

#ifdef _LINUX_
    struct timeval end2;
    gettimeofday(&end2, NULL);
    long long time_use2 = (end2.tv_sec - end1.tv_sec) * 1000000 + (end2.tv_usec - end1.tv_usec);
    long long time_snapshot = (start.tv_sec - snapshot.tv_sec) * 1000000 + (start.tv_usec - snapshot.tv_usec);
    printf("Subgraph time: %.2f s\nPrune time:  %.2f s\nRefine time: %.2f s\n", (double)time_snapshot / 1000000, (double)time_use1 / 1000000, (double)time_use2 / 1000000);
#else
    clock_t end2 = clock();
    printf("Prune time: %.2f s\nSort time: %.2f s\n", (double)(end1 - start) / CLOCKS_PER_SEC, (double)(end2 - end1) / CLOCKS_PER_SEC);
#endif

    cluster_noncore_vertices(eps_a2, eps_b2, miu);
}

void Graph::cluster_noncore_vertices(int eps_a2, int eps_b2, int mu) {
    if (cid == nullptr) cid = new int[n];
    for (int i = 0; i < n; i++) cid[i] = n;

    for (int i = 0; i < n; i++)
        if (similar_degree[i] >= mu) {
            int x = find_root(i);
            if (i < cid[x]) cid[x] = i;  // 如果聚类中有更小的顶点i，那么类别号cid用i表示，所有父结点仍然为x
        }

    noncore_cluster.clear();
    noncore_cluster.reserve(n);
    for (int i = 0; i < n; i++)
        if (similar_degree[i] >= mu) {  // 只有核心才有资格聚类自己的邻居
            for (int j = pstart[i]; j < pstart[i + 1]; j++)
                if (similar_degree[edges[j]] < mu) {  // 核心顶点已经抱团 只需要找非核心的
                    if (min_cn[j] >= 0) {             // 还未确定相似性
                        min_cn[j] = similar_check_OP(i, j, eps_a2, eps_b2);
                        if (reverse[reverse[j]] != j) printf("WA cluster_noncore\n");
                        min_cn[reverse[j]] = min_cn[j];
                        if (min_cn[j] == -1) {
                            ++similar_degree[i];
                            ++similar_degree[edges[j]];
                        }
                    }
                    // 记录 非core 的点属于哪个聚类
                    if (min_cn[j] == -1) noncore_cluster.push_back(make_pair(cid[pa[i]], edges[j]));  // edges[j]是非core
                }
        }
}

void Graph::output(const string& path) {
    printf("\t*** Start write result into disk!\n");
    string out_name = path + "aaa.txt";

    // string out_name = path + "result-" + eps + "-" + to_string(miu) + "-" + to_string(ts) + "-" + to_string(te) + ".txt";
    FILE* fout = fopen(out_name.c_str(), "w");

    if (fout == NULL) {
        printf("Can not open file: %s\n", out_name.c_str());
        exit(1);
    }
    fprintf(fout, "c/n     v_id     c_id\n");

    int mu = miu;
    for (int i = 0; i < n; i++)
        if (similar_degree[i] >= mu) {
            fprintf(fout, " c  %8d %8d\n", i, cid[pa[i]]);  // 聚类内 父结点pa 和 聚类号不同
        }

    sort(noncore_cluster.begin(), noncore_cluster.end());  // 排序是为了去重
    // 去重 重复点移到后面 释放
    noncore_cluster.erase(unique(noncore_cluster.begin(), noncore_cluster.end()), noncore_cluster.end());
    for (int i = 0; i < noncore_cluster.size(); i++) {
        fprintf(fout, " n  %8d %8d\n", noncore_cluster[i].second, noncore_cluster[i].first);
    }

    fclose(fout);
}

void Graph::get_eps(const string eps_s) {
    int i = 0, eps_a = 0, eps_b = 1;
    while (eps_s[i] != '\0' && eps_s[i] != '.') {
        eps_a = eps_a * 10 + (eps_s[i] - '0');
        ++i;
    }

    if (eps_s[i] == '.') {
        ++i;
        while (eps_s[i] != '\0') {
            eps_a = eps_a * 10 + (eps_s[i] - '0');
            eps_b *= 10;
            ++i;
        }
    }

    if (eps_a > eps_b || eps_b > 100 || eps_a <= 0) {
        printf("??? Wrong eps format: %d/%d, %s\n", eps_a, eps_b, eps_s.c_str());
        exit(1);
    }

    eps_a2 = eps_a * eps_a;
    eps_b2 = eps_b * eps_b;
}

void Graph::prune_and_cross_link(int eps_a2, int eps_b2, int miu, int* cores, int& cores_e) {
    for (int i = 0; i < n; i++) {                          // must be iterating from 0 to n-1
        for (int j = pstart[i]; j < pstart[i + 1]; j++) {  // i点所有邻居
            if (edges[j] < i) {                            // 避免重复检查，邻居节点顺序是排好的
                if (min_cn[j] == 0) min_cn[j] = -2;        // 对于 邻边j<i，不考虑(j,i),只考虑(i,j)再反向赋值，直接定定性为-2 不相似
                continue;                                  // this edge has already been checked
            }

            int v = edges[j];                  // v是i的邻居
            int a = degree[i], b = degree[v];  // 此时度是加一的
            if (a > b) swap(a, b);
            // 定理6.1 剪枝规则
            // 近似计算，判断是否没有相似的可能
            if (((long long)a) * eps_b2 < ((long long)b) * eps_a2) {
                min_cn[j] = -2;  // 不相似

                --effective_degree[i];
                --effective_degree[v];
            } else {
                int c = compute_common_neighbor_lowerbound(a, b, eps_a2, eps_b2);  // 相似所需最小的公共邻居

                if (c <= 2) {        // 由于cn(u,v)>=2 若所需 c<=2 直接判定uv相似
                    min_cn[j] = -1;  // 相似

                    ++similar_degree[i];
                    ++similar_degree[v];

                    if (similar_degree[i] == miu) cores[cores_e++] = i;  // cores依次记录核心顶点id
                    if (similar_degree[v] == miu) cores[cores_e++] = v;  // cores_e核顶点数量
                } else {
                    min_cn[j] = c;
                }
            }

            if (min_cn[j] != -2) {  // 若还有机会相似，则将该边的反向边也赋值
                // int r_id = binary_search(edges, pstart[v], pstart[v + 1], i);

                int r_id = pstart[v];
                while(edges[r_id] != i) ++r_id;
                
                reverse[j] = r_id;  // 记录反向边的位置
                reverse[r_id] = j;

                min_cn[r_id] = min_cn[j];  // 最低要求是一样的
            }
        }
    }
}

int Graph::similar_check_OP(int u, int idx, int eps_a2, int eps_b2) {
    ++cal_similar_cnt;
    int v = edges[idx];

    if (min_cn[idx] == 0) {  // 相似度还等于0（应该不存在这种情况）
        int du = degree[u], dv = degree[v];
        int c = compute_common_neighbor_lowerbound(du, dv, eps_a2, eps_b2);

        if (c <= 2) return -1;

        min_cn[idx] = min_cn[reverse[idx]] = c;
    }

    // edges_record_.emplace_back(make_tuple(u, v, min_cn[idx]));  // 记录相似度计算

    // return query(u, v, ts, te, min_cn[idx] - 2);  // tpSCAN

    return check_common_neighbor(u, v, min_cn[idx]);  // pSCAN

    // int a = check_common_neighbor(u, v, min_cn[idx]);
    // int b = query(u, v, ts, te, min_cn[idx] - 2);
    // if (a != b) {
    //     printf("\n\t***Index wrong!!!***\n");
    //     printf("u:%d, v:%d, ts:%d, te:%d, cn:%d, a:%d, b:%d\n", u, v, ts, te, min_cn[idx] - 2, a, b);
    //     exit(-1);
    // }
    // return b;
}

int Graph::check_common_neighbor(int u, int v, int c) {
    int cn = 2;

    if (degree[u] > degree[v]) swap(u, v);

    int du = degree[u] + 1, dv = degree[v] + 1;  // du dv多加2， cn-2 <= min{du-2, dv-2}
    int i = pstart[u], j = pstart[v];
    while (i < pstart[u + 1] && j < pstart[v + 1] && cn < c && du >= c && dv >= c) {  // cn < c <= min{du,dv}
        if (edges[i] < edges[j]) {                                                    // c < min{du,dv} 退出while,min{du,dv}是剩下可能有共同点部分
            --du;
            ++i;
        } else if (edges[i] > edges[j]) {
            --dv;
            ++j;
        } else {
            ++cn;
            ++i;
            ++j;
        }
    }

    if (cn >= c) return -1;
    return -2;
}

int Graph::check_common_neighbor(int u, int v, int c) {
    int cn = 2;

    if (degree[u] > degree[v]) swap(u, v);
    int du = degree[u] + 1, dv = degree[v] + 1;
    int i = pstart[u], j = pstart[v];
}

int Graph::my_check_cn(int u, int v) {
    int cn = 0;
    int i = pstart[u], j = pstart[v];
    while (i < pstart[u + 1] && j < pstart[v + 1]) {
        if (edges[i] < edges[j]) {
            ++i;
        } else if (edges[i] > edges[j]) {
            ++j;
        } else {
            ++cn;
            ++i;
            ++j;
        }
    }
    return cn;
}

void Graph::my_check_index(int u, int v) {
    // if (u > v) swap(u, v);
    int cn = my_check_cn(u, v);
    if (cn == 0) return;
    if (cn_t_idx_[u][v].size() < cn) {
        printf("(%d,%d), cn_t_idx_[u][v].size:%d,  cn:%d\n", u, v, (int)cn_t_idx_[u][v].size(), cn);
        exit(-1);
    }
    auto it = upper_bound(cn_t_idx_[u][v][cn - 1].begin(), cn_t_idx_[u][v][cn - 1].end(), make_pair(ts, te), cmp);
    --it;
    if (it->second > te) {
        printf("(%d,%d), ts:%d, te:%d, idx_ts:%d, idx_te:%d, cn=%d\n", u, v, ts, te, it->first, it->second, cn);
        exit(-2);
    }
}

int Graph::find_root(int u) {
    int x = u;
    while (pa[x] != x) x = pa[x];

    while (pa[u] != x) {
        int tmp = pa[u];
        pa[u] = x;
        u = tmp;
    }

    return x;
}

void Graph::my_union(int u, int v) {
    int ru = find_root(u);
    int rv = find_root(v);

    if (ru == rv) return;

    if (rank[ru] < rank[rv])
        pa[ru] = rv;
    else if (rank[ru] > rank[rv])
        pa[rv] = ru;
    else {
        pa[rv] = ru;
        ++rank[ru];
    }
}

int Graph::compute_common_neighbor_lowerbound(int du, int dv, int eps_a2, int eps_b2) {
    // c 所需最小公共邻居 ceil返回浮点数
    int c = ceil(sqrtl((((long double)du) * ((long double)dv) * eps_a2) / eps_b2));
    // 若没有向上取整，则加1
    // if (((long long)c) * ((long long)c) * eps_b2 < ((long long)du) * ((long long)dv) * eps_a2) ++c;
    return c;
}

void Graph::load_idx(const string& path) {
    auto fp = fopen(path.c_str(), "rb");
    if (fp == NULL) {
        printf("Cannot load index!\n");
        exit(1);
    }
    printf("Loading index...\n");
    unsigned int n_;
    fread(&n_, sizeof(unsigned int), 1, fp);

    if (cn_t_idx_ == nullptr) cn_t_idx_ = new unordered_map<int, vector<vector<pair<int, int>>>>[n_];

    for (int u = 0; u < n_; ++u) {
        int u_size;
        fread(&u_size, sizeof(int), 1, fp);
        // cn_t_idx_[u].reserve(u_size);
        for (int j = 0; j < u_size; ++j) {
            int v, cn_size;
            fread(&v, sizeof(int), 1, fp);
            fread(&cn_size, sizeof(int), 1, fp);
            cn_t_idx_[u][v].resize(cn_size);
            for (int k = 0; k < cn_size; ++k) {
                int cnt;
                fread(&cnt, sizeof(int), 1, fp);
                cn_t_idx_[u][v][k].resize(cnt);
                fread(&cn_t_idx_[u][v][k][0], cnt * sizeof(pair<int, int>), 1, fp);
            }
        }
    }
    fclose(fp);
}

int Graph::query(int u, int v, int t_s, int t_e, int k) {
    if (u > v) swap(u, v);
    if (cn_t_idx_[u][v].size() < k) return -2;
    auto it = upper_bound(cn_t_idx_[u][v][k - 1].begin(), cn_t_idx_[u][v][k - 1].end(), make_pair(t_s, t_e), cmp);
    --it;
    if (it->second <= t_e) return -1;  // 相似
    return -2;                         // 不相似
}

bool cmp(const pair<int, int>& a, const pair<int, int>& b) {
    return a.first < b.first;
}