#include "graph.h"

#include <time.h>

#include <math.h>
#include "MyTree.h"

graph canon[MAXN * MAXM];
extern int* subgraphDegree;
//extern bool* IsD;
extern unsigned long long subgraph_THR;
extern float subgraphDensity;
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
//extern bool* IsD;
extern unsigned long long callNautyCount;
extern hash_map<std::string,long long int> graphInt;
extern hash_map<std::string,long long int> treeInt;
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
	currentEdgeNum = 0;// Tao added 2016-1-19

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
void Graph::Nexts(Subgraph *sub, int maxSize, int startChild) {//g->Nexts(sub, subgraphSize, 0);

	//cout<<subgraphCounter<<endl;
	if ((subgraph_THR != 0)&&(subgraphCounter >= subgraph_THR))
	{
		return;
	}//useful
	int *N;
	
	N = getNeighbours(sub->lastVertex);//sub->lastVertex = sub->vertices[0] = v; for v in 1:154

	int addedCounter = 0;
	
	for(int j = N[0]; j > 0; j--)
	{
		if(N[j] < sub->vertices[0])
			continue;
	
		if(!sub->visited[N[j]]) {
			sub->visited[N[j]] = true;
			sub->AddChild(N[j]);//3 arrays: visited[154],vertices[4], children[154]
			addedCounter++;
		}
	}

	for(int c = startChild; c < sub->childCounter; c++)
	{
		
		sub->lastVertex = sub->vertices[sub->subgraphSize] = sub->children[c];
		sub->subgraphSize++;

		if((sub->subgraphSize == maxSize)&&(subgraphCounter < subgraph_THR))// useful
		{

			subgraphCounter++;
			
			int subEdgeNum = 0;
			for(int i = 0; i < maxSize; i++)
			{
				for(int j = i+1; j < maxSize; j++)
				{
					if(isConnected(sub->vertices[i], sub->vertices[j]))
					{
						subEdgeNum += 1;
					}
				}
			}

			if(subEdgeNum == maxSize -1)//if(subEdgeNum == maxSize -1)// if this subgraph is a tree
			{
				std::string treeStr = GetTreeString(sub->vertices);
				// if (treeStr == "1<11<111>")
				// {
				// 	for (int i = 0; i < maxSize; ++i)
				// 	{
				// 		cout<< sub->vertices[i]<<" ";
				// 	}
				// 	cout<<endl;

				// }
				hash_map<std::string, long long int>::iterator iter = treeInt.find(treeStr);
				if(iter == graphInt.end())
				{
					treeInt[treeStr] = 1;
				}
				else
				{
					(iter->second) += 1;
				}
				//cout<<"treeStr: "<< treeStr<<endl;
			}
			else if (float(subEdgeNum)/maxSize < subgraphDensity)//check graph density
			{
				//cout<<"I'm here, I'm alive!"<<endl;
				std::string adjMatStr = GetAdjMatString(sub->vertices);
				//cout<<adjMatStr<<endl;
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
string Graph::GetTreeString(unsigned int* subVertices)
{
	// cout<< "Entering: GetTreeString, subVertices: "<<endl;
	// for (int i = 0; i < subgraphSize; ++i)
	// {
	// 	cout<<subVertices[i]<<" ";
	// }
	// cout<<endl;
	//Find the root, o
	
	register int i,j;
	int* tempSubgraph = new int[subgraphSize];
	for (int i = 0; i < subgraphSize; ++i)
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
				subgraphDegree[tempSubgraph[j]] += 1;//undirected!
			}
		}
	}
	int* subVerticesDegree =  new int[subgraphSize];
	for (int i = 0; i < subgraphSize; ++i)
	{
		subVerticesDegree[i] = subgraphDegree[tempSubgraph[i]];
		//cout<<subVerticesDegree[i]<<" ";
	}
	//cout<<endl;
	int remainNum = subgraphSize;
	int toDelNode;
	while(remainNum>2)
	{
		std::vector<int> deleted_nodes;
		for (int i = 0; i < subgraphSize; ++i)
		{
			if (subgraphDegree[tempSubgraph[i]] == 1 && tempSubgraph[i]!=-1 ) // leaf node with degree ==1
			{
				// remove tempSubgraph[i]
				toDelNode = tempSubgraph[i];
				tempSubgraph[i] = -1;
				remainNum -= 1;
				deleted_nodes.push_back(toDelNode);
				subgraphDegree[toDelNode] = 0;
			}
		
		}
		//update the degree of remaining nodes
		std::vector<int>::iterator it;
		for (it = deleted_nodes.begin(); it < deleted_nodes.end(); ++it)
		{
			int deleted_node = *it;
			//cout<<"d"<<deleted_node;
			for (int i = 0; i < subgraphSize; ++i)
			{
				if (tempSubgraph[i]!=-1 && isConnected(tempSubgraph[i], deleted_node))
				{
					subgraphDegree[tempSubgraph[i]] -= 1;
				}
			}

		}
		//cout<<endl;
		// for (int i = 1; i < subgraphSize+1; ++i)
		// {
		// 	cout<<subgraphDegree[i]<<"s ";
		// 	/* code */
		// }
		// cout<<endl;
	}
	int root1 = -1, root2 = -1;
	int root1Ind = -1, root2Ind = -1;
	//cout<<"remainNum: "<<remainNum<<endl;
	if (remainNum ==1 )
	{
		for (int i = 0; i < subgraphSize; ++i)
		{
			if (tempSubgraph[i]!=-1)
			{
				root1 = tempSubgraph[i];
				root1Ind = i;
			}
		}
	}
	else if (remainNum ==2)
	{
		for (int i = 0; i < subgraphSize; ++i)
		{
			if (tempSubgraph[i]!=-1 && root1 == -1)
			{
				root1 = tempSubgraph[i];
				root1Ind = i;
			}
			else if (tempSubgraph[i]!= -1 && root1 != -1)
			{
				root2 = tempSubgraph[i];
				root2Ind = i;
			}
		}
		if (subVerticesDegree[root1Ind] > subVerticesDegree[root2Ind]) // select root which has more children
		{
			root2 = -1;
		}
		else if (subVerticesDegree[root2Ind] > subVerticesDegree[root1Ind])
		{
			root1 = -1;
		}

	}


	// cout<<"root1: " <<root1<<" root2: "<<root2<<endl;
	// cout<<"root1Ind: " <<root1Ind<<" root2Ind: "<<root2Ind<<endl;

	string outStr1 = ">", outStr2 = ">";
	//Call Tree part, input a layered matrix, o as the root
	if (root1 == -1 && root2 == -1)
	{
		cout<<"ERROR: root1 == -1 and root2 == -1"<<endl;
		exit(1);
	}
	if (root1 != -1)
	{
		outStr1 = subGetTreeString(subVertices, root1, root1Ind);
		
	}
	if (root2 != -1)
	{
		outStr2 = subGetTreeString(subVertices, root2, root2Ind);
	}

	delete [] tempSubgraph;

	if (outStr1 <= outStr2)
	{
		return outStr1;
	}
	else
	{
		return outStr2;
	}
	
}

string Graph::subGetTreeString(unsigned int* subVertices, int root, int rootInd )
{
		
		MyTree tree(subgraphSize);
	
		bool used[subgraphSize];
		for(int i = 0; i< subgraphSize; i++)
		{
			used[i] = false;
		}
		list<int> queue;
  		queue.push_back(1);//always transform the root as 1
  		int lookup[subgraphSize+1];
  		lookup[0] = 0;//no sense
  		lookup[1] = root;
  		used[rootInd] = true;// mark root as used, root's index is root1Ind
  		int ind = 2; // other nodes coding from 2

  		while(!queue.empty())
  		{
  			int k = queue.front();
  			queue.pop_front();
  			for (int i = 0; i < subgraphSize; ++i)
  			{
  				if (isConnected(lookup[k], subVertices[i]) && used[i]==false)
  				{
  					tree.addEdge(k, ind);
  					queue.push_back(ind);
  					used[i] = true;
  					lookup[ind] = subVertices[i];
  					ind++;
  					
  				}
  			}

  		}

  		treeChildrenSize = new int[subgraphSize+1];
  		for(int i = 1; i< subgraphSize +1; i++)
		{
		   treeChildrenSize[i] = tree.ChildrenNum(i);
		}

		BFSVec = new string[subgraphSize+1];
		
		string outStr = tree.TreeCamGen(subgraphSize); //tree1.TreeCamGen(int** subgraphTree, int level);
		
		delete [] BFSVec;
		delete [] treeChildrenSize;
		return outStr;
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
	adjMatStr.reserve(subgraphSize*(subgraphSize-1)/2 );
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
				subgraphDegree[tempSubgraph[j]] += 1;//undirected!
			}
		}
	}

	qsort(tempSubgraph, subgraphSize, sizeof(tempSubgraph[0]), cmp);//smaller degree -> bigger degree
	
	// for(i = 1; i < subgraphSize; i++)
	// {
	// 	for(j = 0; j < i; j++)
	// 	{
	// 		if ( isConnected(tempSubgraph[i], tempSubgraph[j])  )
	// 		{
	// 			adjMatStr.push_back('1');
	// 		}
	// 		else
	// 		{
	// 			adjMatStr.push_back('0');
	// 		}
	// 	}
	// }
	for(i = 0; i < subgraphSize; i++)
	{
		for(j = i+1; j < subgraphSize; j++)
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



	return adjMatStr;
}
/****************************************************************
****************************************************************/
Graph::Graph(const int n, int k) {
	register int i, j;
	subgraphSize = k;

	M = ((subgraphSize + WORDSIZE - 1) / WORDSIZE);
	nauty_g = new graph[subgraphSize * MAXM];//graph=int, A 'graph' consists of n contiguous sets.  The i-th set represents 
											//   the vertices adjacent to vertex i, for i = 0,1,...,n-1.
	lab = new int[subgraphSize];
	ptn = new int[subgraphSize];
	orbits = new int[subgraphSize];


	printingTime = 0;
	head = 1;
    nV = n;
	nE = 0;
	nEd = 0;
	E_temp.resize(nV+1);
	h = sizeof(Entry) << 3;// what is h for?

	options.writeautoms = FALSE;
	options.getcanon = TRUE;
	options.defaultptn = TRUE;//wangtao 2015-12-19
	options.digraph = TRUE;
	
	// a char array for storing adjacent matrix
	// rowSize indicates how many neighbours does a node maximaly have
	rowSize = (int)ceil((float)nV/8);
	adjMat = new char[rowSize*(nV+1)+1];
	
	// initialize each node's neighbour
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
	
	// int *temp = new int[max+1];
	// for(int jj = 0; jj < max; jj++) {
	// 	temp[jj] = 0;	
	// }
	
	// for(int jj =0; jj < degs.size(); jj++) {
	// 	temp[degs[jj]]++;
	// }

 //    delete []temp;

	printf("Number of Nodes: %d\n", nV);
	printf("Number of Edges: %d\n", nE);
	printf("Maximum Degree: %d\n", directed?maxDegree:maxDegree/2);
}

