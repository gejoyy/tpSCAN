#include "graph.hpp"

using namespace std;

int main(int argc, char* argv[]) {

    // 1:-idx  2:grapg_path  3:index_path
    if (strcmp(argv[1], "-idx") == 0) {
        string graph_path(argv[2]);  
        string idx_path(argv[3]);

        Graph G;
        G.init_log("./log.txt");
        G.load(graph_path); 

        G.index();
        G.write_idx(idx_path);
        // G.write_idx_txt(idx_path);
    }
    
    // 1:-idx-bl  2:grapg_path  3:index_path
    if (strcmp(argv[1], "-idx-bl") == 0) {
        string graph_path(argv[2]);
        string idx_path(argv[3]);

        Graph G;
        G.init_log("./log-bl.txt");
        G.load(graph_path); 

        G.index_baseline();
        G.write_idx(idx_path);
        // G.write_idx_txt(idx_path);  
    }

    return 0;
}
