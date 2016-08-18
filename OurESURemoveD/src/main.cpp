//latest version
// modified on Jun 26, 2016
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <math.h>
#include "graph.h"
#include <getopt.h>

#include <time.h>
#include <list>

#include "MyTree.h"
#include <ext/hash_map>

#define EPOC 10
//#define Debug

float frequency_thr = 0;
unsigned long long subgraph_THR = 0xffffffffffff;
int* subgraphDegree;
hash_map<char,int> AsciiToInt;
hash_map<int,char> IntToAscii;
hash_map<std::string,long long int> graphInt;
hash_map<std::string,long long int> treeInt;
using namespace std;

int subgraphSize = -1, num_random_graphs = 0;
float subgraphDensity = 2;
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
	fprintf (stream, "Usage: OurESU options[inputfile...] \n ");
    fprintf (stream,
		 "\t-h	--help\t\t\tDisplay this usage information. \n"
		 "\t-i	--input filename\tInput filename.\n"
		 "\t-o	--output path\t\tOutput directory.\n"
		 "\t-r 	--random number \tNumber of random graphs (default = 0).\n"
		 "\t-s 	--size motifsize \tMotif size.\n"
		 "\t-u 	--undirected \tUndirected input network\n"
		 "\t-l	--limit \tlimitation of Enumerated subgraphs\n"
		 "\t-d	--density \tdensity of subgraph\n"
		 "\t-f	--threshold \tthreshold of frequency\n");
	     
    exit (exit_code);
}

/****************************************************************
****************************************************************/

bool ReadData(const char *path) {
	
	

	register int i, j;
	int graphSize;
	FILE * inFile = fopen(path, "r");
	
	if (inFile == NULL) {
		printf("Error opening %s file.\n", path);
		return false;
	}
	
	if (!feof(inFile))
		fscanf(inFile, "%d\n", &graphSize);	


	subgraphDegree = new int[graphSize + 1];
	
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

	Subgraph* sub = new Subgraph(1, subgraphSize, g->Size());// creat 3 arrays: visited[154],vertices[4], children[154]
	for (v = 1; v <= g->Size(); v++)
	{
#ifdef Debug
		printf("+ Exploring Node %d ...\n", v);
#endif Debug

		// if(g->subgraphCounter >= 100)
		// 	break;//useful

		sub->subgraphSize = 1;
		sub->lastVertex = sub->vertices[0] = v;

		sub->currentEdgeNum = 0;// Tao add
		
		sub->visited[v] = true;
		//cout << "Before Nexts, everything is OK"<<endl;
		//g->Nexts(sub, subgraphSize, 0, g->quaT->root); 
		g->Nexts(sub, subgraphSize, 0);

		sub->visited[v] = false;
		
	}

	delete sub;
}
/****************************************************************
****************************************************************/
std::string graphDegreeSequence(string adj, int subgraphSize)
{

    register int i, j, index = 0;
    //int totalLength = adj.size();
    vector<int> degreeVec;
    for(i = 0; i<subgraphSize; i++)
        degreeVec.push_back(0);



    // for(i = 0; i < subgraphSize; i++)
    // {
    //     for(j = i+1; j<subgraphSize; j++)
    //     {
    //     	index = j*(j-1)/2+i;
    //         if(adj[index] == '1')
    //         {
    //             degreeVec[i] += 1;
    //             degreeVec[j] += 1;
    //         }
    //         //index ++;
    //     }
    // }
    for(i = 0; i < subgraphSize; i++)
    {
        for(j = i+1; j<subgraphSize; j++)
        {
            if(adj[index] == '1')
            {
                degreeVec[i] += 1;
                degreeVec[j] += 1;
            }
            index ++;
        }
    }
   
    std::sort(degreeVec.begin(), degreeVec.end());
   
    string ds = "";
    ds.reserve(subgraphSize);
    for(vector<int>::iterator it = degreeVec.begin(); it != degreeVec.end(); it++)
    {
        ds.push_back(IntToAscii.find(*it-1)->second); 
    }

   
    return ds;
}
/****************************************************************
****************************************************************/
string calculateCam(string graphIn, int subgraphSize)
{
    //cout<<"graphIn"<<graphIn<<endl;
    string tgraph(subgraphSize*subgraphSize,'0');
   

    register int i, j, index;
    // for(i = 0; i < subgraphSize; i++)
    // {
    //     for (j = i+1; j<subgraphSize; j++)
    //     {
    //         index = j*(j-1)/2 +i;
    //         if(graphIn[index] == '1')
    //         {
    //             tgraph[i*subgraphSize + j] = '1';
    //             tgraph[j*subgraphSize + i] = '1';//undirected
    //         }
    //     }
    // }
    index = 0;
    for(i = 0; i < subgraphSize; i++)
    {
        for (j = i+1; j<subgraphSize; j++)
        {
            if(graphIn[index] == '1')
            {
                tgraph[i*subgraphSize + j] = '1';
                tgraph[j*subgraphSize + i] = '1';//undirected
            }
            index++;
        }
    }

    int  n, m, v, k;
    unsigned nCanCode;
    char sCanCode[subgraphSize*subgraphSize];
    set *gv;
    graph g[subgraphSize*subgraphSize];
    graph canong[subgraphSize*subgraphSize];
    nvector lab[subgraphSize],ptn[subgraphSize],orbits[subgraphSize];
    static DEFAULTOPTIONS(options);
    setword workspace[160*subgraphSize];
    /*init for nauty*/
    options.writeautoms = FALSE;
    options.writemarkers = FALSE;
    options.getcanon = TRUE;
    options.defaultptn = TRUE;//wangtao
    options.digraph = TRUE;// FALSE->True,modified by wangtao 4/20th/2015
    statsblk(stats);
    

    n = subgraphSize;
    m = (n+WORDSIZE-1)/WORDSIZE;

    for(v=0; v<n; v++)
    {
        gv=GRAPHROW(g,v,m);
        EMPTYSET(gv,m);
        for(i=0; i<n; i++)
            if(tgraph[v*subgraphSize+i] != '0')//here make sure that the graph is undirected!!!no,depends on tgraph.
                ADDELEMENT(gv,i);
    }//initialize matrix g

   
  
    nauty(g,lab,ptn,NILSET,orbits,&options,&stats,workspace,160*subgraphSize,m,n,canong);
    nCanCode = 0;
    k=0;
    for(i=0; i<n; i++)
    {
        gv=GRAPHROW(canong,i,m);
        for(j=0; j<n; j++)
        {
            nCanCode = nCanCode<<1;
            nCanCode+=ISELEMENT(gv,j);
            sCanCode[k++] = (char)(ISELEMENT(gv,j)+48);//ascii(0) == 48
        }
    }//output cam string
    sCanCode[k]='\0';
    string re = sCanCode;
    
    return re;
}
/****************************************************************
****************************************************************/

int main(int argc, char *argv[]) {

	// for (int ini = 0; ini < 9; ini++)
 //    {
 //        char temp = '1'+ini;
 //        AsciiToInt[temp] = ini;
 //        IntToAscii[ini] = temp;
 //    }
 //    for (int ini = 0; ini < 26; ini++)
 //    {
 //        char temp = 'a'+ini;
 //        AsciiToInt[temp] = 9 + ini;
 //        IntToAscii[9+ini] = temp;
 //    }
 //    for (int ini = 0; ini < 26; ini++)
 //    {
 //        char temp = 'A'+ini;
 //        AsciiToInt[temp] = 35 + ini;
 //        IntToAscii[35+ini] = temp;
 //    }
	for(int ini = 0; ini < 94; ini++)
	{
		char temp = '!' + ini;
		AsciiToInt[temp] = ini;
		IntToAscii[ini] = temp;
	}


	double total_random_time = 0 , main_time;
	clock_t start_random_time, end_random_time;
	directed = true;
	register int i, j;
	long long unsigned subgraphCounterMain;

	int next_option;
	const char *const short_options = "h:i:o:r:s:d:f:l:u";
	const struct option long_options[] = {
		{"help",   0, NULL, 'h'},
		{"input",  1, NULL, 'i'},
		{"output", 1, NULL, 'o'},
		{"random", 1, NULL, 'r'},
		{"size",   1, NULL, 's'},
		{"density", 1, NULL, 'd'},
		{"frequency", 1, NULL, 'f'},
		{"limit", 1, NULL, 'l'},
		{"undirected",   0, NULL, 'u'},
		{NULL,     0, NULL,  0 }
	};
	
	char *program_name;
    char input_filename[256], output_directory[256];

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
	    	
			case 'o':
				strcpy(output_directory, optarg);
	    		break;
			
			case 'r':
				 num_random_graphs = atoi(optarg);
	    		break;

			case 's':
				subgraphSize = atoi(optarg);
	    		break;
	    	case 'd':
	    		subgraphDensity = atof(optarg);
	    		break;
	    	case 'f':
	    		frequency_thr = atof(optarg);
	    		break;
	    	case 'l':
	    		subgraph_THR = atof(optarg);
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
    cout<<"subgraph_THR = " << subgraph_THR<<endl;
    directed = false;//Tao 2016 Jun 25, only work on undirected graph now.

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
	
	if (!ReadData(input_filename))
		return 0;

	printf("ReadData Finished!!!\n");
	cout<<"Nodes:"<<g->nV<<endl;
	cout<<"Edges:"<<g->Edges()<<endl;
	g->setPath(output_directory);

	clock_t startTime = clock();
	//for main graph
	isRand = false;
	g->subgraphCounter = 0;
	g->notClassSubCounter = 0;
	clock_t mainStartTime = clock();
	//printf("before enumerate, everything is OK!\n");
	Enumerate();

	//printf("after enumerate, everything is OK!!\n");
	clock_t mainEndTime = clock();

	//hash_map<std::string, long long int> degreeSeqCount;
	hash_map<string, pair< vector<const string*>, long long int> > degreeSeqPair;
	hash_map<string, pair< vector<const string*>, long long int> >::iterator itor;
	cout<<"graphIntSize = "<<graphInt.size()<<endl;
	for(hash_map<std::string, long long int>::iterator it = graphInt.begin(); it!= graphInt.end(); it++)
	{
		const std::string *graphAdjMatStr = &(it->first);
		//cout<< *graphAdjMatStr<<endl;
		std::string graphDegSeq = graphDegreeSequence(*graphAdjMatStr, subgraphSize);
		//cout<<"No problem"<<endl;
		long long int tempCount = it->second;
		
		itor = degreeSeqPair.find(graphDegSeq);

		if( itor == degreeSeqPair.end())
		{
			vector<const string*> tmpstringvec;
			tmpstringvec.push_back(graphAdjMatStr);
			pair< vector<const string*>, long long int > tmpPair = make_pair(tmpstringvec, tempCount);

			degreeSeqPair[graphDegSeq] = tmpPair;
		}
		else
		{
			(itor->second).first.push_back(graphAdjMatStr);
			(itor->second).second += tempCount;
		}

	}

	//cout<<"still alive!"<<endl;

	hash_map<string, long long int> finalGraph;
	for(itor = degreeSeqPair.begin(); itor != degreeSeqPair.end(); itor++)
	{
		if((itor->second).second <= frequency_thr )
		{
			continue;
		}
		for(int iv = 0; iv < (itor->second).first.size(); iv++ )
		{
			string tempCam = *(((itor->second).first)[iv]);
			string cam = calculateCam(tempCam, subgraphSize);
			callNautyCount += 1;

			hash_map<std::string, long long int>::iterator it2 =  finalGraph.find(cam);
            if (it2 == finalGraph.end())
            {

                finalGraph[cam] = graphInt.find(tempCam)->second;
            }
            else
            {
                (it2->second) += graphInt.find(tempCam)->second;
            }
		}
	}

	//count graph part
	int graph_class = 0, tree_class = 0;
	for(hash_map<std::string, long long int>::iterator it = finalGraph.begin(); it!= finalGraph.end(); it++)
	{
		if (it->second >= frequency_thr)
		{
			graph_class += 1;
		}
	}	

	//count tree part
	hash_map<string, long long int> finalTree;
	for(hash_map<std::string, long long int>::iterator it = treeInt.begin(); it!= treeInt.end(); it++)
	{
		if (it->second >= frequency_thr)
		{
			tree_class += 1;
			//cout<< it->first<<endl;
		}
	}


	subgraphCounterMain = g->subgraphCounter;
	enumerated_class = graph_class + tree_class;

	
	delete g;
	//delete IsD;
	delete subgraphDegree;

	clock_t endTime = clock();
	main_time = difftime(mainEndTime, mainStartTime)/(double)CLOCKS_PER_SEC;
	double total_time = difftime(endTime, startTime)/(double)CLOCKS_PER_SEC;

	printf("\n===========RESULTS===========\n");
	
	printf("\nMotif Size: %d\n", subgraphSize);
	printf("\nTotal number of subgraphs in original network: %llu\n", subgraphCounterMain);
	printf("Number of frequent non-isomorphic classes: %lu\n", enumerated_class);
	printf("Number of frequent non-isomorphic graph classes: %lu\n", graph_class);
	printf("Number of frequent non-isomorphic tree classes: %lu\n", tree_class);
	printf("\nTime Used for Enumerate:      %f\n", main_time); 
	
	printf("\nTotal Time Used: %f\n", total_time); 
	printf("\n=============================\n");
	printf("Call Nauty %lld\n", callNautyCount );
	
	return 0;
}
