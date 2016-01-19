#include "graph.h"
#include "QuaternaryTree.h"
#include <time.h>

graph canon[MAXN * MAXM];


static DEFAULTOPTIONS(options);
statsblk(stats);
setword workspace[160*MAXM];

long long unsigned * C_main;
long long unsigned * C_rand;
double * mean;
double * var;
double* Score;
long long unsigned int* ID;
int idxID;
unsigned long head;
unsigned long enumerated_class;

FILE * o; 

extern bool isRand;
extern bool* IsD;
extern unsigned long long callNautyCount;
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
Subgraph::Subgraph(int subgraphSize, int maxSize, int graphSize) {
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
void Graph::Nexts(Subgraph *sub, int maxSize, int startChild, QNode* curNode) {
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
		
		// quaT->cur_node = curNode;
		// for(int i = 0; i < sub->subgraphSize-1; i++)
		// {
		// 	if(isConnected(sub->vertices[i], sub->lastVertex))
		// 	{
		// 		if(isConnected(sub->lastVertex, sub->vertices[i]))//two way connection
		// 		{
		// 			if(quaT->cur_node->two == NULL)
		// 				quaT->cur_node->two = quaT->add_node();
		// 			quaT->cur_node = quaT->cur_node->two;
		// 		}
		// 		else
		// 		{
		// 			if(quaT->cur_node->pone == NULL)
		// 				quaT->cur_node->pone = quaT->add_node();
		// 			quaT->cur_node = quaT->cur_node->pone;
		// 		}
		// 	}
		// 	else if(isConnected(sub->lastVertex, sub->vertices[i]))
		// 	{
		// 		if(quaT->cur_node->mone == NULL)
		// 			quaT->cur_node->mone = quaT->add_node();
		// 		quaT->cur_node = quaT->cur_node->mone;
		// 	}
		// 	else
		// 	{
		// 		if(quaT->cur_node->zero == NULL)
		// 			quaT->cur_node->zero = quaT->add_node();
		// 		quaT->cur_node = quaT->cur_node->zero;
		// 	}
		// }

		if(sub->subgraphSize == maxSize)
		{
			subgraphCounter++;
			// quaT->Check();

			//calling Nauty if the leaf is new
			// if(quaT->cur_node->leaf == NULL)
			// {
			// 	quaT->cur_node->leaf = Classify(sub);
			// 	if(quaT->cur_node->leaf->tag)//for random networks and subgraphs that do not exist in the original network
			// 		notClassSubCounter++;
			// }
		}
		else
			Nexts(sub, maxSize, c+1, quaT->cur_node);

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
Graph::Graph(const int n, int k) {
	register int i, j;
	subgraphSize = k;
	T = new tree(k);
	M = ((subgraphSize + WORDSIZE - 1) / WORDSIZE);
	nauty_g = new graph[subgraphSize * MAXM];
	lab = new int[subgraphSize];
	ptn = new int[subgraphSize];
	orbits = new int[subgraphSize];
	quaT = new QTree();

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
    delete [] ID;
    delete [] mean;
    delete [] var;
    delete [] Score;
    delete [] C_rand;
    delete [] C_main;
    delete [] degree;
    delete T;
	delete quaT;
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

void Graph::deleteEdgeAdjMat(vertex u, vertex v) {
	adjMat[(rowSize*u) + (v>>3)] &= (~(1 << (v % 8)));
}

/****************************************************************
****************************************************************/

void Graph::swapEdge(vertex v, int ind, vertex u) {	
	if(u < E[v][ind]) {
		ind++;
		while(ind <= E[v][0] && u < E[v][ind]) {
			E[v][ind-1] = E[v][ind];
			ind++;
		}
		if(ind <= E[v][0] && u == E[v][ind]) {
			while(ind <= E[v][0]) {
				E[v][ind-1] = E[v][ind];
				ind++;
			}
			E[v][0]--;
		}
		else
			E[v][ind-1] = u;
	}
	else {
		ind--;
		while(ind > 0 && u > E[v][ind]) {
			E[v][ind+1] = E[v][ind];
			ind--;
		}
		if(ind > 0 && u == E[v][ind]) {
			ind++;
			ind++;
			while(ind <= E[v][0]) {
				E[v][ind-1] = E[v][ind];
				ind++;
			}
			E[v][0]--;
		}
		else
			E[v][ind+1] = u;
	}	
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
int Graph::get_vertex() {
	int ind = rand() % nEd;
	return degree[ind];
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

/****************************************************************
This funciton calls Nauty and uses the returned canonical labet to
traverse the Binary Tree

For subgraphs of random networks the Binary Tree is not modified,
it is just checked to see if the subgraph is classified before in
the original network.
****************************************************************/
Leaf* Graph::Classify(Subgraph *sub) {
	register int i = 0, j, l, k;
    set *gv;
	int subgraphSize = sub->subgraphSize;
	
	for (i = 0; i < subgraphSize; i++) {
		gv = GRAPHROW(nauty_g, i, M);
		EMPTYSET(gv, M);		
		for(j = 0; j < subgraphSize; j++) {
			if(i == j)
				continue;
			if (isConnected(sub->vertices[i], sub->vertices[j])) 
			{
				ADDELEMENT(gv, j);
			}
		}
	}
	
	/***Initial ptn and lab Add by wangtao on 1/8th/2014***/
   
    for(i = 0; i < subgraphSize-1; i++)
    {
        lab[i] = i+1;
    }
    lab[subgraphSize-1] = 0;
    for(i = 0; i< subgraphSize-2; i++)
    {
        ptn[i] = 1;
    }
    ptn[subgraphSize - 2] = 0;
    ptn[subgraphSize - 1] = 0;
    /***/

	nauty(nauty_g, lab, ptn, NULL, orbits, &options, &stats, 
		  workspace, 160*MAXM, M, subgraphSize, canon);
	callNautyCount += 1;
	
	T->init_cur_node();

	if (!isRand) {
		for (i = 0; i < subgraphSize-1; i++) {	
			for(j = 0; j < subgraphSize; j++) {
				if(i == j)
					continue;
				if(isConnected(sub->vertices[lab[i]], sub->vertices[lab[j]])) 
					T->insert_one_main();
				else
					T->insert_zero_main();
			}
		}
		
		for(j = 0; j < subgraphSize-2; j++) {
			if(isConnected(sub->vertices[lab[i]], sub->vertices[lab[j]])) 
				T->insert_one_main();
			else
				T->insert_zero_main();
		}
		
		if(isConnected(sub->vertices[lab[i]], sub->vertices[lab[j]])) 
			return T->update_one_main(1);
		else 
			return T->update_zero_main(1);
	}
	else {
		for (i = 0; i < subgraphSize-1; i++) {	
			for(j = 0; j < subgraphSize; j++) {
				if(i == j)
					continue;
				if(isConnected(sub->vertices[lab[i]], sub->vertices[lab[j]])) 
					if(!T->insert_one_rand())
						return T->trueLeaf;
				else
					if(!T->insert_zero_rand())
						return T->trueLeaf;
			}
		}
		
		for(j = 0; j < subgraphSize-2; j++) {
			if(isConnected(sub->vertices[lab[i]], sub->vertices[lab[j]])) 
				if(!T->insert_one_rand())
					return T->trueLeaf;
			else
				if(!T->insert_zero_rand())
					return T->trueLeaf;
		}
		
		if(isConnected(sub->vertices[lab[i]], sub->vertices[lab[j]])) 
			return T->update_one_rand(1);
		else 
			return T->update_zero_rand(1);
	}
}

/****************************************************************
Memory Allocation
****************************************************************/
void Graph::AllocateCounter() {
	int class_num = T->get_leafnum();
	C_main = new long long unsigned[class_num + 1];
	C_rand = new long long unsigned[class_num];
	C_main[0] = 0;

	mean = new double[class_num];
	var = new double[class_num];
	Score = new double[class_num];
	
	register int i;
	for(i = 0; i < class_num; i++) {
		mean[i] = 0.0;
		var[i] = 0.0;
		C_rand[i] = 0;
	}
	
	idxID = 0;
	ID = new long long unsigned int[class_num];
	for(int i = 0; i < class_num; i++)
		ID[i] = 0;
}

/****************************************************************
****************************************************************/
void Graph::DFS(Node * cur) {
	if(!cur->left && !cur->right) {
		Leaf * leaf = (Leaf *)cur;
		if(leaf->count > 0) {
			C_rand[head] = leaf->count;
			head++;
			leaf->count = 0;
		}
		return;
	}
	if(cur->left)
		DFS(cur->left);
	if(cur->right)
		DFS(cur->right);
} 

/****************************************************************
****************************************************************/
void Graph::print_adjMatrix(char * str) {
	register int i, j;
	int l = 0;
	int index = 0;
	int maxpow = subgraphSize * subgraphSize -1;

	for(i = 0; i < subgraphSize; i++) {
		for(j = 0; j < subgraphSize; j++) {
			if(i == j) {
				l++;
				fprintf(am,"0");
			}
			else
			{
				index = i*(subgraphSize)+(j-l);
				fprintf(am, "%d", str[index]);
				if (str[index] == 1)
					ID[idxID] += (long long unsigned int) (pow(2, (maxpow - (i*subgraphSize+j))));
			}
		}
		fprintf(am,"\n");
	}
	fprintf(am,"ID: %llu", ID[idxID]);
	fprintf(am,"\n\n");
	idxID++;
}

/****************************************************************
****************************************************************/

void Graph::DFSmain(Node * cur, char * str, int lev) {
	if(!cur->left && !cur->right) {
		//printf("checkpoint1\n");
		clock_t tempTime = clock();
		print_adjMatrix(str);
		printingTime += difftime(clock(), tempTime)/(double)CLOCKS_PER_SEC;
		Leaf * leaf = (Leaf *)cur;
		head++;
		C_main[head] = leaf->count;
		leaf->count = 0;
		C_main[0]++;
		return;
	}
	if(cur->left) {
		str[lev] = 0;
		DFSmain(cur->left, str, lev+1);
	}
	if(cur->right) {
		str[lev] = 1;
		DFSmain(cur->right, str, lev+1);
	}
}

/****************************************************************
Extracting the data required for calculating Z-Scores
****************************************************************/
void Graph::Extract() {
	
	register int i;
	int class_num = T->get_leafnum();
	char * adj_str = new char[subgraphSize*(subgraphSize-1)];
	Node * current = T->return_root();
	
	head = 0;
	if(isRand){
		DFS(current);
		for(i = 0; i < class_num; i++) {
			mean[i] += C_rand[i];
			var[i] += (C_rand[i]*C_rand[i]);
		}
	}	
	else {
		DFSmain(current, adj_str, 0);
		enumerated_class = C_main[0];	
		printf("Number of Non-isomorphic Classes: %lu\n", enumerated_class);
	}
}

/****************************************************************
****************************************************************/
void Graph::calculateZSCORE(int RAND, int subgraphCounter, char *path) {
	FILE * cm;
	int i;
	for (i = 0; i < T->get_leafnum(); i++) {
		mean[i] = mean[i]/RAND;
		var[i] = sqrt((var[i]-(RAND*(mean[i]*mean[i])))/RAND);
		if(var[i] != 0)
			Score[i] = ((double)C_main[i+1] - mean[i])/var[i];
		else
			Score[i] = -1.0;
	}
    
    char file[256];
    sprintf(file, "%s/ZScore.txt", path);
    printf("Writing ZScores to %s ...\n", file); 
	cm = fopen(file, "w+");
    if(!cm) {
        printf("Can't open %s\n", path);
        sprintf(file, "result/ZScore.txt", path);
        printf("Writing to %s instead ...\n", file); 
    	cm = fopen(file, "w+");
        if (!cm) {
            printf("Error again ... Sorry!\n");
            exit(-1);
        }
    }

	fprintf(cm, "TOTAL NUMBER OF CLASSES:: %lu\n\n", enumerated_class);
	fprintf(cm, "ID \t\t NUM IN REAL \t\t MEAN IN RANDOM \t VAR IN RANDOM \t\t ZSCORE\n");
	for (i = 0; i < T->get_leafnum(); i++) {
		if (var[i] != 0)
		{
			fprintf(cm, "%llu \t\t ", ID[i]);
			fprintf(cm, "%llu \t\t ", C_main[i+1]);
			fprintf(cm, "%f \t\t ", mean[i]);
			fprintf(cm, "%f \t\t ", var[i]);
			fprintf(cm, "%f\n", Score[i]);
		}
		if (var[i] == 0) 
		{
			fprintf(cm, "**%llu \t\t ", ID[i]);
			fprintf(cm, "%llu \t\t ", C_main[i+1]);
			fprintf(cm, "%f \t\t ", mean[i]);
			fprintf(cm, "%f \t\t ", Score[i]);
			fprintf(cm, "%f\n", (double)C_main[i+1] - mean[i]);
		}
	}

	fclose (cm);
}
