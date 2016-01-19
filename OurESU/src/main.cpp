//latest version
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

float THR = 0;
bool *IsD;
int* subgraphDegree;
hash_map<char,int> AsciiToInt;
hash_map<int,char> IntToAscii;
hash_map<std::string,long long int> graphInt;

using namespace std;

int subgraphSize = -1, num_random_graphs = 0;
unsigned long long callNautyCount = 0;

//g stores the input graph
Graph *g;

//isRand determines whether the enumeration of random networks is commenced or not
//directed indicates whether the input network is directed or not
bool isRand, directed;

extern unsigned long enumerated_class;
char ascii[61] = {'1','2','3','4','5','6','7','8','9','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};


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
	subgraphDegree = new int[graphSize + 1];
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

	Subgraph* sub = new Subgraph(1, subgraphSize, g->Size());// creat 3 arrays: visited[154],vertices[4], children[154]
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

    register int i, j, index;
    //int totalLength = adj.size();
    vector<int> degreeVec;
    for(i = 0; i<subgraphSize; i++)
        degreeVec.push_back(0);

    int Dindex = AsciiToInt.find(adj[adj.size()-1])->second;//for only one D node.


    for(i = 0; i < subgraphSize; i++)
    {
        for(j = i+1; j<subgraphSize; j++)
        {
            index = j*(j-1)/2 +i;
            if(adj[index] == '1')
            {
                degreeVec[i] += 1;
                degreeVec[j] += 1;
            }
        }
    }
    int Ddegree = degreeVec[Dindex]; 
    std::sort(degreeVec.begin(), degreeVec.end());
    degreeVec.push_back(Ddegree);// put D's degree at the back of degree sequence.

    string ds = "";
    ds.reserve(subgraphSize + 2);
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
   
    // int edgeN = AsciiToInt.find(graphIn[0])->second +1;
    // //printf("edgeN = %d\n",edgeN );
    // int row,col;
    // int c = 2*edgeN ;
    // for(int h = 1; h <= c; h++)
    // {
    //     row = AsciiToInt.find(graphIn[h])->second;
    //     h++;
    //     col = AsciiToInt.find(graphIn[h])->second;
    //     tgraph[row*subgraphSize+col]='1';
    //     tgraph[col*subgraphSize+row]='1';
    // }
    register int i, j, index;
    for(i = 0; i < subgraphSize; i++)
    {
        for (j = i+1; j<subgraphSize; j++)
        {
            index = j*(j-1)/2 +i;
            if(graphIn[index] == '1')
            {
                tgraph[i*subgraphSize + j] = '1';
                tgraph[j*subgraphSize + i] = '1';
            }
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
    options.defaultptn = FALSE;//wangtao
    options.digraph = TRUE;// FALSE->True,modified by wangtao 4/20th/2015
    statsblk(stats);
    

    n = subgraphSize;
    m = (n+WORDSIZE-1)/WORDSIZE;

    for(v=0; v<n; v++)
    {
        gv=GRAPHROW(g,v,m);
        EMPTYSET(gv,m);
        for(i=0; i<n; i++)
            if(tgraph[v*subgraphSize+i] != '0')//here make sure that the graph is undirected!!!
                ADDELEMENT(gv,i);
    }

    /***Initial ptn and lab Add by wangtao on 1/7/2014***/
    //assert(graphIn.size() == 2*edgeN + subgraphSize + 2);
    //char Dnode1 = graphIn[2*edgeN + subgraphSize +1];
    char Dnode1 = graphIn[graphIn.size()-1];

    int lab_i = 0;
    for(i = 0; i < subgraphSize; i++)
    {
        if(i != AsciiToInt[Dnode1])
        {
            lab[lab_i] = i;
            lab_i ++;
        }
    }
    lab[lab_i] = AsciiToInt[Dnode1];
    for(i = 0; i< subgraphSize-2; i++)
    {
        ptn[i] = 1;
    }
    ptn[subgraphSize - 2] = 0;
    ptn[subgraphSize - 1] = 0;
    /**/

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
            sCanCode[k++] = (char)(ISELEMENT(gv,j)+48);
        }
    }
    sCanCode[k]='\0';
    string re = sCanCode;
    
    return re;
}
/****************************************************************
****************************************************************/

int main(int argc, char *argv[]) {

	for (int ini = 0; ini < 9; ini++)
    {
        char temp = '1'+ini;
        AsciiToInt[temp] = ini;
        IntToAscii[ini] = temp;
    }
    for (int ini = 0; ini < 26; ini++)
    {
        char temp = 'a'+ini;
        AsciiToInt[temp] = 9 + ini;
        IntToAscii[9+ini] = temp;
    }
    for (int ini = 0; ini < 26; ini++)
    {
        char temp = 'A'+ini;
        AsciiToInt[temp] = 35 + ini;
        IntToAscii[35+ini] = temp;
    }



	double total_random_time = 0 , main_time;
	clock_t start_random_time, end_random_time;
	directed = true;
	register int i, j;
	long long unsigned subgraphCounterMain;

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
	printf("before enumerate, everything is OK!\n");
	Enumerate();
	printf("after enumerate, everything is OK!!\n");
	clock_t mainEndTime = clock();

	//hash_map<std::string, long long int> degreeSeqCount;
	hash_map<string, pair< vector<const string*>, long long int> > degreeSeqPair;
	hash_map<string, pair< vector<const string*>, long long int> >::iterator itor;

	for(hash_map<std::string, long long int>::iterator it = graphInt.begin(); it!= graphInt.end(); it++)
	{
		const std::string *graphAdjMatStr = &(it->first);
		std::string graphDegSeq = graphDegreeSequence(*graphAdjMatStr, subgraphSize);
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


	hash_map<string, long long int> finalGraph;
	for(itor = degreeSeqPair.begin(); itor != degreeSeqPair.end(); itor++)
	{
		if((itor->second).second <= THR*(g->subgraphCounter) )
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

	subgraphCounterMain = g->subgraphCounter;
	enumerated_class = finalGraph.size();

	
	delete g;
	delete IsD;
	delete subgraphDegree;

	clock_t endTime = clock();
	main_time = difftime(mainEndTime, mainStartTime)/(double)CLOCKS_PER_SEC;
	double total_time = difftime(endTime, startTime)/(double)CLOCKS_PER_SEC;

	printf("\n===========RESULTS===========\n");
	
	printf("\nMotif Size: %d\n", subgraphSize);
	printf("\nTotal number of subgraphs in original network: %llu\n", subgraphCounterMain);
	printf("Number of non-isomorphic classes: %lu\n", enumerated_class);
	printf("\nTime Used for Enumerate:      %f\n", main_time); 
	
	printf("\nTotal Time Used: %f\n", total_time); 
	printf("\n=============================\n");
	printf("Call Nauty %lld\n", callNautyCount );
	
	return 0;
}
