#include "graph.h"

#include <time.h>

#include <math.h>
#include "MyTree.h"

graph canon[MAXN * MAXM];
extern int* subgraphDegree;
extern bool* IsD;
extern char ascii[];
int* treeChildrenSize;
string* BFSVec;

static DEFAULTOPTIONS(options);
statsblk(stats);
setword workspace[160*MAXM];



int idxID;
unsigned long head;
unsigned long enumerated_class;

FILE * o; 

extern bool isRand;
extern bool* IsD;
extern unsigned long long callNautyCount;
extern hash_map<std::string,long long int> graphInt;
/****************************************************************
****************************************************************/
bool sortCmp(int i, int j) {
	return (j < i);
}

/****************************************************************
****************************************************************/
bool Subgraph::Print() {
	
	for(int i = 0; i < subgraphSize; i++)
		printf("%d ", vertices[i]);
	printf("\n");

	char ch, temp;
	scanf("%c", &ch);
	scanf("%c", &temp);
	return(ch == 'c');
}

/****************************************************************
****************************************************************/
Subgraph::Subgraph(int subgraphSize, int maxSize, int graphSize) {//initialize = 1, 4(subgraphsize), 154(graphSize)
	this->subgraphSize = subgraphSize;
	childCounter = 0;
	visited = new bool[graphSize+1];
	vertices = new unsigned int[maxSize];
	children = new unsigned int[graphSize];


	for(int i = 0; i <= graphSize; i++)
		visited[i] = false;
}

/****************************************************************
****************************************************************/
Subgraph::~Subgraph() {
	delete[] visited;
	delete[] vertices;
}

/****************************************************************
****************************************************************/
void Subgraph::AddChild(int vertex) {
	children[childCounter++] = vertex;
}

/****************************************************************
This function is responsible for enumerating and partially classifying the subgraphs with 'quaT' (Quaternary Tree)
****************************************************************/
void Graph::Nexts(Subgraph *sub, int maxSize, int startChild) {
	int *N;
	N = getNeighbours(sub->lastVertex);

	int addedCounter = 0;
	for(int j = N[0]; j > 0; j--)
	{
		if(N[j] < sub->vertices[0])
			continue;
		/***add by wangtao 2015-7-6***/
		if(IsD[N[j]])
			continue;
		/*end*/
		if(!sub->visited[N[j]]) {
			sub->visited[N[j]] = true;
			sub->AddChild(N[j]);
			addedCounter++;
		}
	}

	for(int c = startChild; c < sub->childCounter; c++)
	{
		sub->lastVertex = sub->vertices[sub->subgraphSize] = sub->children[c];
		sub->subgraphSize++;
		
		if(sub->subgraphSize == maxSize)
		{
			subgraphCounter++;
			
			int subEdgeNum = 0;
			// for(int i = 0; i < maxSize; i++)
			// {
			// 	for(int j = i+1; j < maxSize; j++)
			// 	{
			// 		if(isConnected(sub->vertices[i], sub->vertices[j]))
			// 		{
			// 			subEdgeNum += 1;
			// 		}
			// 	}
			// }
			
			if(subEdgeNum == maxSize -1)// if this subgraph is a tree
			{
				;
			}
			else
			{
				
				std::string adjMatStr = GetAdjMatString(sub->vertices);
				//cout<<adjMatStr;
				hash_map<std::string, long long int>::iterator iter = graphInt.find(adjMatStr);
				if(iter == graphInt.end())
				{
					graphInt[adjMatStr] = 1;
				}
				else
				{
					(iter->second) += 1;
				}
			}

		}
		else
			Nexts(sub, maxSize, c+1);

		sub->subgraphSize--;
		sub->lastVertex = sub->vertices[subgraphSize-1];
	}

	//Removing added children
	for(int i = sub->childCounter - addedCounter; i < sub->childCounter; i++)
		sub->visited[sub->children[i]] = false;
	sub->childCounter -= addedCounter;
}
/****************************************************************
****************************************************************/
int cmp(const void *a, const void *b )
{
	return subgraphDegree[*(int*)a] - subgraphDegree[*(int*)b];
}
string Graph::GetAdjMatString(unsigned int* subVertices)// here, adjMatStr can be transfered into a unsigned interger
{
	register int i,j;
	int tempSubgraph[subgraphSize];
	string adjMatStr = "";
	adjMatStr.reserve(subgraphSize*(subgraphSize-1)/2 + 2);
	for(i = 0; i < subgraphSize; i++)
	{
		subgraphDegree[subVertices[i]] = 0;
		tempSubgraph[i] = subVertices[i];
	}


	for(i= 0; i < subgraphSize; i++)
	{
		for(j = i+1; j < subgraphSize; j++)
		{
			if ( isConnected(tempSubgraph[i], tempSubgraph[j]) )
			{
				subgraphDegree[tempSubgraph[i]] += 1;
				subgraphDegree[tempSubgraph[j]] += 1;
			}
		}
	}

	qsort(tempSubgraph, subgraphSize, sizeof(tempSubgraph[0]), cmp);

	for(i = 1; i < subgraphSize; i++)
	{
		for(j = 0; j < i; j++)
		{
			if ( isConnected(tempSubgraph[i], tempSubgraph[j])  )
			{
				adjMatStr.push_back('1');
			}
			else
			{
				adjMatStr.push_back('0');
			}
		}
	}

	for (i = 0; i < subgraphSize; i++){
		
	   if(IsD[tempSubgraph[i]])
		adjMatStr.push_back(ascii[i]);
	}
	

	return adjMatStr;
}
/****************************************************************
****************************************************************/
Graph::Graph(const int n, int k) {
	register int i, j;
	subgraphSize = k;

	M = ((subgraphSize + WORDSIZE - 1) / WORDSIZE);
	nauty_g = new graph[subgraphSize * MAXM];
	lab = new int[subgraphSize];
	ptn = new int[subgraphSize];
	orbits = new int[subgraphSize];


	printingTime = 0;
	head = 1;
    nV = n;
	nE = 0;
	nEd = 0;
	E_temp.resize(nV+1);
	h = sizeof(Entry) << 3;

	options.writeautoms = FALSE;
	options.getcanon = TRUE;
	options.defaultptn = FALSE;//wangtao 2015-7-5
	options.digraph = TRUE;
	
	rowSize = (int)ceil((float)nV/8);
	adjMat = new char[rowSize*(nV+1)+1];
	
	for(i = 1; i <= nV; i++) {
		for(j = 0; j <= (nV >> 3); j++) {
			adjMat[i*rowSize + j] = 0;
		}
	}	
}

/****************************************************************
****************************************************************/

Graph::~Graph() {
    printf("Destroying graph ...\n");
	delete [] adjMat;
    for(int i = 1; i <= nV; i++) {
        delete [] E[i];    
    }
    delete [] E;
    delete [] lab;
    delete [] ptn	;
    delete [] orbits;
 	delete [] nauty_g;


    delete [] degree;

    fclose(am);
}


/****************************************************************
****************************************************************/

void Graph::setPath(char *path) { 
    char file[256];
    sprintf(file, "%s/adjMatrix.txt", path);
    printf("Opening %s ...\n", file); 
    fflush(stdout);
	am = fopen(file, "w+");
    if(!am) {
        printf("Can't open %s\n", path);
        sprintf(file, "result/adjMatrix.txt", path);
        printf("Writing to %s instead ...\n", file); 
    	am = fopen(file, "w+");
        if (!am) {
            printf("Error again ... Sorry!\n");
            exit(-1);
        }
    }
}


/****************************************************************
****************************************************************/

void Graph::Print() {
	register int i, j;
	
	fflush(stdout);
	for(i = 1; i <= nV; i++) {
		printf("Node %d: (%d) ", i, E[i][0]);
		for(j = 1; j <= E[i][0]; j++) {
			printf("%d ", E[i][j]);
		}
		printf("\n");
	}
	printf("---------------------------------------------------------------\n");
	
}

/****************************************************************
****************************************************************/

void Graph::addEdgeAdjMat(vertex u, vertex v) {
	adjMat[(rowSize*u) + (v>>3)] |= (1 << (v % 8));
}




/****************************************************************
****************************************************************/

void Graph::addEdge(vertex u, vertex v) {
	nE++;
	E_temp[u].push_back(v);
	E_temp[v].push_back(u);
	adjMat[(rowSize*u) + (v>>3)] |= (1 << (v % 8));
}

/****************************************************************
****************************************************************/

bool Graph::isConnected(vertex u, vertex v) {
	return adjMat[(rowSize*u) + (v>>3)] & (1 << (v % 8));
}

/****************************************************************
****************************************************************/

int* Graph::getNeighbours(vertex v) {
	return E[v];
}



/****************************************************************
****************************************************************/
void Graph::Finalize(bool directed) {
	register int i, j, max = 0;
	vector<int>::iterator it;
	vector<int> degs;
	
	E = new int*[nV+1];
	degree = new int[nE];
	int degInd = 0;
	
	for(i = 1; i <= nV; i++) {
		sort(E_temp[i].begin(), E_temp[i].end(), sortCmp);
		it = unique(E_temp[i].begin(), E_temp[i].end());		
		E_temp[i].resize(it - E_temp[i].begin());
		
		degs.push_back(E_temp[i].size());

		if(max < E_temp[i].size())
			max = E_temp[i].size();
		
		E[i] = new int[E_temp[i].size() + 1];
		E[i][0] = E_temp[i].size();
		for(j = 0; j < E_temp[i].size(); j++) {
			E[i][j+1] = E_temp[i][j];
			if(isConnected(i, E_temp[i][j])) {
				degree[degInd] = i;
				degInd++;
				nEd++;
			}
		}
		E_temp[i].resize(0);
		E_temp[i].clear();
	}
	maxDegree = max;
	
	int *temp = new int[max+1];
	for(int jj = 0; jj < max; jj++) {
		temp[jj] = 0;	
	}
	
	for(int jj =0; jj < degs.size(); jj++) {
		temp[degs[jj]]++;
	}

    delete []temp;

	printf("Number of Nodes: %d\n", nV);
	printf("Number of Edges: %d\n", nE);
	printf("Maximum Degree: %d\n", directed?maxDegree:maxDegree/2);
}

