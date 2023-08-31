#include "tpSCAN.hpp"

using namespace std;

void usage() {
	printf("Usage:[0]exe [1]graph-path [2]index-path [3]epsilon [4]mu [5]start-time [6]end-time\n");
}

int main(int argc, char* argv[]) {
	if(argc < 7) {
		usage() ;
		return 0;
	}

    printf("\n\t*** Graph Clustering (eps:%s, mu:%s) ***\n", argv[3], argv[4]);

    Graph G;
	G.load_idx(argv[2]);
    G.read_graph(argv[1]);

#ifdef _LINUX_
	struct timeval start, end;
	gettimeofday(&start, NULL);
#else
	clock_t start, end;
	start = clock();
#endif
    G.tpSCAN(argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));

#ifdef _LINUX_
	gettimeofday(&end, NULL);
	long long use_time = (end.tv_sec - start.tv_sec)*1000000 + (end.tv_usec - start.tv_usec);
	printf("Total time without IO: %.2f s\n", (double)use_time/1000000);
#else
	end = clock();
	printf("Total time without IO: %.2f s\n",(double)(end - start) / CLOCKS_PER_SEC);
#endif

    G.output("./result/");


    // G.tpSCAN(argv[2], atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
    // // G.output("./result/");

    return 0;
}
