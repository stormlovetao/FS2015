//latest version
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <math.h>
#include "graph.h"
#include <getopt.h>
#include "randomGenerator.h"
#include <time.h>
#include <list>

#define EPOC 10
//#define Debug
bool *IsD;
using namespace std;

int subgraphSize = -1, num_random_graphs = 0;
unsigned long long callNautyCount = 0;

//g stores the input graph
Graph *g;

//isRand determines whether the enumeration of random networks is commenced or not
//directed indicates whether the input network is directed or not
bool isRand, directed;

extern unsigned long enumerated_class;

/****************************************************************
****************************************************************/

void print_usage (FILE * stream, int exit_code)
{
	fprintf (stream, "Usage: Kavosh options[inputfile...] \n ");
    fprintf (stream,
		 "\t-h	--help\t\t\tDisplay this usage information. \n"
		 "\t-i	--input filename\tInput filename.\n"
		 "\t-o	--output path\t\tOutput directory.\n"
		 "\t-r 	--random number \tNumber of random graphs (default = 0).\n"
		 "\t-s 	--size motifsize \tMotif size.\n"
		 "\t-u 	--undirected \tUndirected input network\n");
	     
    exit (exit_code);
}

/****************************************************************
****************************************************************/

bool ReadData(const char *path, int numRandomGraphs, const char* Dpath) {
	
	/****************Start: Add on 29th/11/2013 by wangtao****************/
	list<int> Dlist;
	//char dnodefile[] = "./network/D_nodes_exp.txt";
	//char dnodefile[] = "./network/paper_D.txt";
	 FILE * DnodeFile = fopen(Dpath,"r");
    if (DnodeFile == NULL)
    {
        printf("Error opening %s file.\n", Dpath);
        return false;
    }
	
	int dnodelabel;
	while(!feof(DnodeFile))
	{
		fscanf(DnodeFile, "%d\n", &dnodelabel);
		Dlist.push_back(dnodelabel);
	}
	/****************End: Add on 29th/11/2013 by wangtao****************/

	register int i, j;
	int graphSize;
	FILE * inFile = fopen(path, "r");
	
	if (inFile == NULL) {
		printf("Error opening %s file.\n", path);
		return false;
	}
	
	if (!feof(inFile))
		fscanf(inFile, "%d\n", &graphSize);	


	/****************Start: Add on 29th/11/2013 by wangtao****************/
	list<int>::iterator DlistIterator;
	IsD = new bool[graphSize + 1];
	for (i = 1; i <= graphSize; i++)
	{
		DlistIterator = find(Dlist.begin(), Dlist.end(), i);
		if(DlistIterator == Dlist.end())
		{
			IsD[i] = false;
		}
		else
		{
			IsD[i] = true;
		}
	}
	/****************End: Add on 29th/11/2013 by wangtao****************/

	
	g = new Graph(graphSize, subgraphSize);
	while (!feof(inFile)) {
		fscanf(inFile, "%d %d\n", &i, &j);
		if(i == j) continue;
		g->addEdge(i, j);
		if(!directed)
			g->addEdge(j, i);
	}
	
	g->Finalize(directed);
	fclose(inFile);

	return true;
}

/***********************************************************************************
 * This function enumerates the subgraphs related to each vertex of inpur network. *
***********************************************************************************/
void Enumerate() {
	register int v;
	//g->preComC = 0;

	Subgraph* sub = new Subgraph(1, subgraphSize, g->Size());
	for (v = 1; v <= g->Size(); v++)
	{
#ifdef Debug
		printf("+ Exploring Node %d ...\n", v);
#endif Debug

		/*****Add by Wangtao on 29th/11/2013*****/
		if(!IsD[v])
		{
			break;
		}
		/****END***/


		sub->subgraphSize = 1;
		sub->lastVertex = sub->vertices[0] = v;
		
		sub->visited[v] = true;
		g->Nexts(sub, subgraphSize, 0, g->quaT->root);
		sub->visited[v] = false;
	}

	delete sub;
}

/****************************************************************
****************************************************************/

int main(int argc, char *argv[]) {
	double total_random_time = 0 , main_time;
	clock_t start_random_time, end_random_time;
	directed = true;
	register int i, j;
	long long unsigned subgraphCounterMain;
	generator gen;
	int next_option;
	const char *const short_options = "h:i:d:o:r:s:u";
	const struct option long_options[] = {
		{"help",   0, NULL, 'h'},
		{"input",  1, NULL, 'i'},
		{"Dnodes input",  1, NULL, 'd'},
		{"output", 1, NULL, 'o'},
		{"random", 1, NULL, 'r'},
		{"size",   1, NULL, 's'},
		{"undirected",   0, NULL, 'u'},
		{NULL,     0, NULL,  0 }
	};
	
	char *program_name;
    char input_filename[256], output_directory[256], Dpath[256];

    int verbose = 0;
    strcpy(output_directory, "result");

    program_name = argv[0];
    do {
		next_option = getopt_long (argc, argv, short_options, long_options, NULL);
	
		switch (next_option) {
			case 'h':
	    		print_usage (stdout, 0);

			case 'i':
				strcpy(input_filename, optarg);
	    		break;
	    	case 'd':
				strcpy(Dpath, optarg);
	    		break;

			case 'o':
				strcpy(output_directory, optarg);
	    		break;
			
			case 'r':
				 num_random_graphs = atoi(optarg);
	    		break;

			case 's':
				subgraphSize = atoi(optarg);
	    		break;

			case 'u':
				directed = false;
	    		break;
			
			case '?':
	    		print_usage (stdout, 1);
				
			case -1:		/* Done with options. */
			    break;
			
			default:		/* Something else: unexpected. */
                print_usage (stderr, -1);
		}
    } while (next_option != -1);


	if (input_filename == NULL) {
		fprintf(stderr, "Input Argument Error: Please specify an input filename using \"-i <FILE NAME>\".\n");
        print_usage (stderr, -1);
	}
	
	if (subgraphSize == -1) {
		fprintf(stderr, "Input Argument Error: Please specify a motif size using \"-s <MOTIF SIZE>\".\n");
        print_usage (stderr, -1);
	}
	
	printf("Motif Size: %d\n", subgraphSize);
	printf("Input Graph: %s\n", input_filename);
		
	if (!ReadData(input_filename, num_random_graphs, Dpath))
		return 0;

	printf("ReadData Finished!!!\n");

	g->setPath(output_directory);

	clock_t startTime = clock();
	//for main graph
	isRand = false;
	g->subgraphCounter = 0;
	g->notClassSubCounter = 0;
	clock_t mainStartTime = clock();
	Enumerate();
	clock_t mainEndTime = clock();
	g->AllocateCounter();
	//printf("Total Number of Subgraphs: %ld\n", g->subgraphCounter);
	printf("\nTotal Number of Subgraphs: %f\n", (float)g->subgraphCounter);
	g->Extract();
	subgraphCounterMain = g->subgraphCounter;

	//initializing random graphs
	srand(time(NULL));
	isRand = true;
	long long unsigned AvgTotRndSubCounter = 0;
	long long unsigned AvgRndSubCounter = 0;
	printf("\nNumber of Random Graphs: %d\n", num_random_graphs);
	for (i = 1; i <= num_random_graphs; i++) {
		gen.genRandGraph_Edge3(g);
		g->subgraphCounter = 0;
		g->notClassSubCounter = 0;
		clock_t randStartTime = clock();
		Enumerate();

		clock_t randEndTime = clock();
		total_random_time += difftime(randEndTime, randStartTime);

		g->Extract();
    	printf("\nTotal Number of Subgraphs in Random graph %d:      %llu\n", i, g->subgraphCounter);
		printf("Number of Classified Subgraphs in Random graph %d: %llu\n", i, g->subgraphCounter-g->notClassSubCounter);
		AvgTotRndSubCounter += g->subgraphCounter;
		AvgRndSubCounter += g->subgraphCounter-g->notClassSubCounter;
	}

	printf("\nRandoms Completed!!!\n");

	if (0 < num_random_graphs)
		g->calculateZSCORE(num_random_graphs, subgraphCounterMain, output_directory);
	
	delete g;

	clock_t endTime = clock();
	main_time = difftime(mainEndTime, mainStartTime)/(double)CLOCKS_PER_SEC;
	double total_time = difftime(endTime, startTime)/(double)CLOCKS_PER_SEC;

	printf("\n===========RESULTS===========\n");
	printf("\nTopography of Original Network: %d nodes and %d edges\n", g->Size(), directed?g->Edges():g->Edges()/2);
	printf("\nMotif Size: %d\n", subgraphSize);
	printf("\nTotal number of subgraphs in original network: %llu\n", subgraphCounterMain);
	printf("Number of non-isomorphic classes: %lu\n", enumerated_class);
	printf("\nTime Used for Main:      %f\n", main_time); 
	printf("Avg Time Used for Rands: %f\n", total_random_time/(num_random_graphs*CLOCKS_PER_SEC)); 
	printf("\nAverage number of all subgraphs in random networks:        %f\n", AvgTotRndSubCounter/(double)num_random_graphs); 
	printf("Average number of classified subgraphs in random networks: %f\n", AvgRndSubCounter/(double)num_random_graphs); 
	printf("\nTotal Time Used: %f\n", total_time); 
	printf("\n=============================\n");
	printf("Call Nauty %lld\n", callNautyCount );
	
	return 0;
}
