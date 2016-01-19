// This main.cpp is used to test the speed between MyTree and Nauty
// ./FSMining -i ./network/V4_net_exp.txt -d ./network/D_nodes_exp.txt -s 9 -f 0/1
#include <ext/hash_map>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <map>
#include <math.h>
#include "graph.h"
#include <getopt.h>
#include <string>
#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <sstream>
#include <exception>
#include <assert.h>
//#include <iostream>
#include "randomGenerator.h"
#include "MyTree.h"
// #include "Tree.h"
// #include "Vertex.h"
#include "Shared.h"
#define EPOC 10
//#define Debug
#define OneD //one D node
//#define NoD //no D node
#define OUTGastonFormat
using namespace std;
using namespace   __gnu_cxx;

void GEN ( int n, int k, int root, int level, int reminder, int m);
void NEG ( int n, int k, int root, int level, int reminder, int m);

//1st dimension is the layer, while the 2nd is the list of selected nodes in each layer. subgraph[k][0] is the number of nodes in layer k.
//vertex subgraph[subgraphSize][subgraphSize+1];
vertex **subgraph;
int subgraphSize = 11, num_random_graphs = 0;
int MM = subgraphSize*subgraphSize;
double THR = 0.0000;
double seqTHR = 0.0000;
long long int countGraph =0;
long long int countGraph2 = 0;
double DENSITY = 0.5*subgraphSize*(subgraphSize-1);
bool *Visited;
bool *IsD;
//bool *IsC; //check if a node is c node
bool containC; //contain a C node or not
int Dcount = 0;
hash_map<std::string,long long int> treeInt;
hash_map<std::string,long long int> graphInt;
hash_map<char,int> AsciiToInt;
hash_map<int,char> IntToAscii;
// This array collects the valid children at each layer.
unsigned int **childSet;

// This array maps the childSet to the subgraph[level]
unsigned int **Index;

//total number of enumerated subgraphs.
unsigned long long subgraphCounter;

//g stores the input graph
Graph *g;

long num;

bool isRand;
int treeflag;


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
             "\t-s 	--size motifsize \tMotif size.\n" );

    exit (exit_code);
}

/****************************************************************
****************************************************************/


ofstream outfile;


bool ReadData(const char *path, const char *Dpath)
{
    register int i, j;
    int graphSize;
    int dnodelabel;
    int cnodelabel;
    list<int> Dlist;//list for D nodes set
    
   
    FILE * DnodeFile = fopen(Dpath,"r");
    if (DnodeFile == NULL)
    {
        printf("Error opening %s file.\n", Dpath);
        return false;
    }
    
    while (!feof(DnodeFile))
    {
        fscanf(DnodeFile, "%d\n", &dnodelabel);
        Dlist.push_back(dnodelabel);
        //cout<<"dnodelabel: "<<dnodelabel<<endl;
    }
    
    list<int>::iterator DlistIterator;
    
    FILE * inFile = fopen(path, "r");

    if (inFile == NULL)
    {
        printf("Error opening %s file.\n", path);
        return false;
    }

    if (!feof(inFile))
        fscanf(inFile, "%d\n", &graphSize);
    IsD = new bool[graphSize+1];
    
    for(i = 1; i <= graphSize; i++)
    {
        DlistIterator = find (Dlist.begin(), Dlist.end(),i);
        if(DlistIterator == Dlist.end())
            IsD[i] = false;
        else
        {
            IsD[i] = true;
            //cout<<"is D node.............: "<<i<<endl;
        }
        
    }
    Visited = new bool[graphSize+1];
    for(i = 1; i <= graphSize; i++)
        Visited[i] = false;


    g = new Graph(graphSize, subgraphSize);
    while (!feof(inFile))
    {
        fscanf(inFile, "%d %d\n", &i, &j);
        if(i == j) continue;
        g->addEdge(i, j);
    }

    g->Finalize();

    //g->Print();

    fclose(inFile);
    fclose(DnodeFile);

    childSet = new unsigned int*[subgraphSize];
    Index    = new unsigned int*[subgraphSize];

    for(i = 0; i < subgraphSize; i++)
    {
        childSet[i] = new unsigned int[g->MaxDegree() * subgraphSize + 1];
        Index[i] = new unsigned int[g->MaxDegree() * subgraphSize + 1];
    }

    return true;
}

/***************************************************************
 * This function finds the valid children in each given level. *
***************************************************************/

void initChildSet(int root, int level)
{
    register int *N;
    register int i, j;
    const int *parents = subgraph[level-1];

    childSet[level][0] = 0;
    for(i = 1; i <= parents[0]; i++)
    {
        N = g->getNeighbours(parents[i]);
        for(j = 1; j <= N[0] && root <= N[j]; j++)
        {
#ifdef OneD
            if((!Visited[N[j]])&&(!IsD[N[j]]))
#endif
#ifdef NoD
                if(!Visited[N[j]])
#endif
                {
                    Visited[N[j]] = true;
                    childSet[level][0]++;
                    childSet[level][childSet[level][0]] = N[j];
                }
        }
    }
}

// bool treeIso(std::string si1, std::string si2)
// {

//     char* t1 = new char[si1.size()];
//     char* t2 = new char[si2.size()];
//     //cout<<"si1:"<<si1<<endl;
//     string s1 = si1.substr(0,si1.length()-1);
//     //cout<<"si2:"<<si2<<endl;
//     string s2 = si2.substr(0,si2.length()-1);

//     string Dnode_id_s1 = si1.substr(si1.length()-1,1);
//     string Dnode_id_s2 = si2.substr(si2.length()-1,1);
//     strcpy(t1,s1.c_str());
//     strcpy(t2,s2.c_str());
// //printf("t1: %s  \n",t1);
// //  printf("t2: %s  \n",t2);
//     Tree * tr1 = new Tree(t1,subgraphSize,AsciiToInt);
//     Tree * tr2 = new Tree(t2,subgraphSize,AsciiToInt);

//     /*find D node, make the D node as root*/
//     int Did=999,Did2=999;
//     for(int i = 0; i < tr1->c_vertices.size(); i++)
//     {
//         if(atoi(tr1->c_vertices[i]->id.c_str()) == AsciiToInt[Dnode_id_s1.at(0)])
//         {
//             Did = i;
//             break;
//         }
//         //cout<<Dnode_id_s1.at(0)<<"xxx"<<AsciiToInt[Dnode_id_s1.at(0)]<<endl;
//         //cout<<atoi(tr1->c_vertices[i]->id.c_str())<<endl;
//     }
//     if(Did > 998)
//         cout<<"error! No D node!"<<endl;
//     //cout<<Did<<endl;


// //  std::vector<Vertex *> centers = tr1->findCenter();
//     /*set the D node as root node*/
//     // tr1->reRoot(centers[0]);
//     tr1->reRoot(tr1->c_vertices[Did]);

//     for(int i = 0; i < tr2->c_vertices.size(); i++)
//     {
//         if(atoi(tr2->c_vertices[i]->id.c_str()) == AsciiToInt[Dnode_id_s2.at(0)])
//         {
//             Did2 = i;
//             break;
//         }
//         //cout<<Dnode_id_s1.at(0)<<"xxx"<<AsciiToInt[Dnode_id_s1.at(0)]<<endl;
//         //cout<<atoi(tr1->c_vertices[i]->id.c_str())<<endl;
//     }
//     if(Did2 > 998)
//         cout<<"error! No D node!"<<endl;
//     // cout<<Did2<<endl;

//     //centers = tr2->findCenter();
//     //tr2->reRoot(centers[0]);
//     tr2->reRoot(tr2->c_vertices[Did2]);
//     std::string mp;
//     mp = tr1->rootedIsomorphic(tr2);
//     /*
//      if (mp.compare("")==0) {
//     	if (centers.size()==2) {
//     		tr2->reRoot(centers[1]);
//     		mp = tr1->rootedIsomorphic(tr2);
//     		}
//     	}
//       */
//     delete[] t1;
//     delete[] t2;
//     /*
//     for(unsigned int i = 0; i < tr1->c_vertices.size(); i++){
//        delete tr1->c_vertices[i];
//     }
//     for(unsigned int i = 0; i < tr2->c_vertices.size(); i++){
//        delete tr2->c_vertices[i];
//     }
//     */
//     delete tr1;
//     delete tr2;

//     //cout<<"mp"<<mp<<endl;
//     if (mp.compare("")!=0)
//     {
//         return true;
//     }
//     else
//     {
//         return false;
//     }
// }
/****************************************************************************************
 * This function Explores the root vertex and generates subgraphs containing that node. *
****************************************************************************************/

void Explore(vertex root, int level, int reminder)
{
    register int i, j, k; // k is the number of selected nodes in the current level.
    /*  if (subgraphCounter >10000) {
    return;
    }  */

#ifdef Debug
    printf("************************************\n");
    printf("*****   Exploring level %3d  *******\n", level);
    printf("************************************\n");
#endif
#ifdef Density_D
    if ((double)(subgraph[level-1][subgraphSize+2]+(subgraphSize-subgraph[level-1][subgraphSize+1]))/DENSITY > 0.25)
    {
        return;
    }
#endif
    if (reminder == 0)   //reminder == 0 assures level <= subgraphSize
    {

        i = 0;

#ifdef Density_D
        int tempSubgraph[subgraphSize];
        for (int l = 0; l < level; l++)
        {
            for(k = 1; k <= subgraph[l][0]; k++)
            {
                tempSubgraph[i++] = subgraph[l][k];
            }
        }
        int edgeCount =0;
        for ( int x = 0; x < subgraphSize; x++)
        {
            for (int y = x+1; y < subgraphSize; y++)
            {
                if (g->isConnected(tempSubgraph[x],tempSubgraph[y]) || g->isConnected(tempSubgraph[y],tempSubgraph[x]))
                {
                    edgeCount++;
                }
            }
        }
#endif



        /*888888888888888888888888888888888888888888888888888888888888888888888*/

      

                if (subgraph[level-1][subgraphSize+2] < subgraphSize)//the graph is a tree actually
                {
                    if(treeflag==1)
                    {
                        string treeCam;
                        treeCam = g->gettree_cam(subgraph, level);
                    }
                    
                
                    

                    else
                        g->Classify(subgraph, level);

                 }
              
                else
                {
                    ;
                }
           
            
               
        
        /*888888888888888888888888888888888888888888888888888888888888888888888*/
        return;
    }


    initChildSet(root, level);

#ifdef Debug
    printf("Valid Children in level %d:\nN = { ", level);
    for(k = 1; k <= childSet[level][0]; k++)
    {
        printf("%d ", childSet[level][k]);
    }
    printf("}\n");
#endif

    for(k = 1; k <= reminder; k++)
    {
        if( childSet[level][0] < k )   //There is not enough child to choose m from.
        {
            for(i = 1; i <= childSet[level][0]; i++)
            {
                Visited[childSet[level][i]] = false;
            }
            return;
        }

#ifdef Debug
        printf("Selecting %d node(s) from level %d\n", k, level);
        printf("Initial Selection = { ");
        for(i = 1; i <= k; i++)
        {
            printf("%d ", childSet[level][i]);
        }
        printf("}\n");
#endif
        subgraph[level][0] = k;
        int addEdge = 0;
        for(i = 1; i <= k; i++)
        {
            register int *tempN;
            subgraph[level][i] = childSet[level][i];
            Index[level][i] = i;
//#ifdef Density_D
            tempN = g->getNeighbours(subgraph[level][i]);
            if (tempN[0]==1)
            {
                addEdge++;
                continue;
            }

            for (int x = 1; x <= subgraph[level-1][0]; x++)
            {
                if (g->isConnected(subgraph[level][i],subgraph[level-1][x]) || g->isConnected(subgraph[level-1][x],subgraph[level][i]))
                {
                    addEdge++;
                }
            }


        }
        for (int y = 1; y <= k; y++)
        {
            for (int z=y+1; z <=k; z++)
            {

                if (g->isConnected(subgraph[level][y],subgraph[level][z]) || g->isConnected(subgraph[level][z],subgraph[level][y]))
                {
                    addEdge++;
                }
            }
//#endif
        }
//#ifdef Density_D
        subgraph[level][subgraphSize+1] = subgraph[level-1][subgraphSize+1]+k;
        subgraph[level][subgraphSize+2] =  subgraph[level-1][subgraphSize+2]+addEdge;
//#endif


        Explore(root, level + 1, reminder - k);
        GEN( childSet[level][0], k, root, level, reminder, k);

#ifdef Debug
        printf("************************************\n");
        printf("*****    Back to level %3d   *******\n", level);
        printf("************************************\n");
#endif
    }

    for(i = 1; i <= childSet[level][0]; i++)
    {
        Visited[childSet[level][i]] = false;
    }
    subgraph[level][0] = 0;
    return;
}

/***************************************************************************************************
 * The following                                                                                *   ee functions generate all C(n, k) in Gray Code order, Adopted from Rusky code. *
***************************************************************************************************/

void swap( int i, int j, int root, int level, int reminder, int m)
{
#ifdef Debug
    printf("Switch %d with %d in level %d\n", childSet[level][j], childSet[level][i], level);
#endif
    int addEdge=0;
    Index[level][i] = Index[level][j];
    subgraph[level][Index[level][i]] = childSet[level][i];

    for(int i2 = 1; i2 <= subgraph[level][0]; i2++)
    {
        int *tempN;

//#ifdef Density_D
        tempN = g->getNeighbours(subgraph[level][i2]);
        if (tempN[0]==1)
        {
            addEdge++;
            continue;
        }

        for (int x = 1; x <= subgraph[level-1][0]; x++)
        {
            if (g->isConnected(subgraph[level][i2],subgraph[level-1][x]) || g->isConnected(subgraph[level-1][x],subgraph[level][i2]))
            {
                addEdge++;
            }
        }
    }
    for (int y = 1; y <= subgraph[level][0]; y++)
    {
        for (int z=y+1; z <=subgraph[level][0]; z++)
        {

            if (g->isConnected(subgraph[level][y],subgraph[level][z]) || g->isConnected(subgraph[level][z],subgraph[level][y]))
            {
                addEdge++;
            }
        }
//#endif
    }
    subgraph[level][subgraphSize+1] = subgraph[level-1][subgraphSize+1]+subgraph[level][0];
    subgraph[level][subgraphSize+2] =  subgraph[level-1][subgraphSize+2]+addEdge;
    Explore(root, level + 1, reminder - m);
}

/****************************************************************
****************************************************************/

void GEN( int n, int k, int root, int level, int reminder, int m)
{
    /* if (subgraphCounter >10000) {
                   return;
                }
      */
    if (k > 0 && k < n)
    {
        GEN( n-1, k, root, level, reminder, m);
        if (k == 1)
            swap( n, n-1, root, level, reminder, m);
        else
            swap( n, k-1, root, level, reminder, m);

        NEG( n-1, k-1, root, level, reminder, m);
    }
}

/****************************************************************
****************************************************************/

void NEG( int n, int k, int root, int level, int reminder, int m)
{
    if (k > 0 && k < n)
    {
        GEN( n-1, k-1, root, level, reminder, m);

        if (k == 1)
            swap( n-1, n, root, level, reminder, m);
        else
            swap( k-1, n, root, level, reminder, m);

        NEG( n-1, k, root, level, reminder, m);
    }

}

/***********************************************************************************
 * This function enumerates the subgraphs related to each vertex of inpur network. *
***********************************************************************************/
void Enumerate()
{
    register int v;

    for (v = 1; v <= g->Size(); v++)
    {
#ifdef Debug
        printf("+ Exploring Node %d ...\n", v);
#endif
        if(!IsD[v])
            break;
        /* if (subgraphCounter >10000) {
            break;
         }
        */
        //cout<<"v:"<<v<<endl;
        subgraph[0][0] = 1;
        subgraph[0][1] = v;
        /*add at Aug. 8 2012*/
        subgraph[0][subgraphSize+1] = 1; //node count
        subgraph[0][subgraphSize+2] = 0; //edge count
        /**/
        Visited[v] = true;
        Explore(v, 1, subgraphSize - 1);
        Visited[v] = false;
    }
}
/****************************************************************
****************************************************************/
std::string graphDegreeSequence(string adj, int subgraphSize, hash_map<char,int> AsciiToInt)
{
    //cout<<"adj"<<adj<<endl;
    //int totalLength = adj.size();
    //int n = subgraphSize*subgraphSize;  edgeNumber*2+subgraphSize+1
    int numEdge = AsciiToInt.find(adj[0])->second+1;
    int totalLength = numEdge*2 + subgraphSize +1;//adj.size();
    //cout<<"..."<<adj.substr(0,totalLength)<<endl;
    char Dnode = adj[totalLength];
    // cout<<"Dnode"<<Dnode<<endl;
    // cout<<adj<<endl;
    string ds(subgraphSize,'0');
    string ds2(numEdge,'0');
    vector<char> vec;
    char tempDegree;
    //vector<char> vec2;
    for (int i = 2*numEdge+1; i < totalLength; i++)
    {
        if(AsciiToInt[Dnode] == (i - (2*numEdge+1)))
        {
            tempDegree = adj[i];
            continue;
        }
        vec.push_back(adj[i]);
        //std::cout<<adj[i*subgraphSize + i]<<" ";
    }
    std::sort(vec.begin(),vec.end());
    /*put the d node in the last position*/
    vec.push_back(tempDegree);
    for (int i = 0; i < subgraphSize; i++)
    {
        ds[i] = vec[i];
    }
    /*degree sequence 2*/
    /*
    int node1, node2;
    for (int i = n; i < totalLength;i++) {
       node1 = AsciiToInt[adj[i]];
       i++;
       node2 = AsciiToInt[adj[i]];
       vec2.push_back(adj[node1*subgraphSize+node1] + adj[node2*subgraphSize+node2]-2);

    }

    std::sort(vec2.begin(),vec2.end());
    for (int i = 0; i < numEdge;i++) {
       ds2[i] = vec2[i];

    }
    */
    // std::cout<<std::endl;
    vec.clear();
    // cout<<"ds"<<ds<<endl;
    return ds;
}
/****************************************************************
****************************************************************/
std::string graphDegreeSequence(string adj, int subgraphSize)
{

    int totalLength = adj.size();
    int n = subgraphSize*subgraphSize;
    int numEdge = (totalLength - n)/2;
    string ds(subgraphSize,'0');
    string ds2(numEdge,'0');
    vector<char> vec;
    vector<char> vec2;
    for (int i = 0; i < subgraphSize; i++)
    {
        vec.push_back(adj[i*subgraphSize + i]);
        //std::cout<<adj[i*subgraphSize + i]<<" ";
    }

    std::sort(vec.begin(),vec.end());
    for (int i = 0; i < subgraphSize; i++)
    {
        ds[i] = vec[i];
    }
    /*degree sequence 2*/
    /*
    int node1, node2;
    for (int i = n; i < totalLength;i++) {
       node1 = AsciiToInt[adj[i]];
       i++;
       node2 = AsciiToInt[adj[i]];
       vec2.push_back(adj[node1*subgraphSize+node1] + adj[node2*subgraphSize+node2]-2);

    }

    std::sort(vec2.begin(),vec2.end());
    for (int i = 0; i < numEdge;i++) {
       ds2[i] = vec2[i];

    }
    */
    // std::cout<<std::endl;
    return ds;
}


string calculateCam(string graphIn, int subgraphSize, hash_map<char,int> AsciiToInt)
{
    //cout<<"graphIn"<<graphIn<<endl;
    string tgraph(subgraphSize*subgraphSize,'0');
   
    int edgeN = AsciiToInt.find(graphIn[0])->second +1;
    //printf("edgeN = %d\n",edgeN );
    int row,col;
    int c = 2*edgeN ;
    for(int h = 1; h <= c; h++)
    {
        row = AsciiToInt.find(graphIn[h])->second;
        h++;
        col = AsciiToInt.find(graphIn[h])->second;
        tgraph[row*subgraphSize+col]='1';
        tgraph[col*subgraphSize+row]='1';
    }
    int i, j, n, m, v, k;
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
    assert(graphIn.size() == 2*edgeN + subgraphSize + 2);
    char Dnode1 = graphIn[2*edgeN + subgraphSize +1];
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


string calculateCam(string graphIn, int subgraphSize)
{
    // cout<<graphIn<<endl;
    // string tgraph(subgraphSize*subgraphSize,'0');
    string tgraph = graphIn.substr(0,subgraphSize*subgraphSize);
    for(int h = 0; h < subgraphSize; h++)
    {
        tgraph[h*subgraphSize+h]='0';
    }
    int i, j, n, m, v, k;
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
    options.defaultptn = TRUE;
    options.digraph = TRUE;
    statsblk(stats);


    n = subgraphSize;
    m = (n+WORDSIZE-1)/WORDSIZE;

    for(v=0; v<n; v++)
    {
        gv=GRAPHROW(g,v,m);
        EMPTYSET(gv,m);
        for(i=0; i<n; i++)
            if(tgraph[v*subgraphSize+i] != '0')
                ADDELEMENT(gv,i);
    }

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

std::vector<std::string> split(std::string const &input)
{
    std::stringstream buffer(input);
    std::vector<std::string> ret;

    std::copy(std::istream_iterator<std::string>(buffer),
              std::istream_iterator<std::string>(),
              std::back_inserter(ret));
    return ret;
}
std::string intToString(int i)
{

    std::string s;
    std::stringstream out;
    out << i;
    s = out.str();
    return s;

}

double realtime()
{
    struct timeval tp;
    struct timezone tzp;
    gettimeofday(&tp, &tzp);
    return tp.tv_sec + tp.tv_usec * 1e-6;
}
/****************************************************************
****************************************************************/

int main(int argc, char *argv[])
{
    

   outfile.open("./tree.txt");

    printf("THR = %f\n", THR);
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
    //cout<<AsciiToInt.find('a')->first<<endl;
    //cout<<AsciiToInt.find('a')->second<<endl;
    double total_random_time = 0;
    double random_time , main_time;
    clock_t start_random_time, end_random_time;
    timeval tv, tv1, tvPartition;
    timeval treeover, DSeqover, Camover;
    register int i, j;
    long long subgraphCounterMain;
    generator gen;
      int next_option;
  const char *const short_options = "hi:d:o:r:s:f:";
  const struct option long_options[] = {
        {"help",   0, NULL, 'h'},
        {"input",  1, NULL, 'i'},
        {"Dnodes input",  1, NULL, 'd'},
        {"output", 1, NULL, 'o'},
        {"random", 1, NULL, 'r'},
        {"size",   1, NULL, 's'},
        {"treeflag",   1, NULL, 'f'},       
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

            case 'f':
                treeflag = atoi(optarg);
                break;
            
            case '?':
                print_usage (stdout, 1);
                
            case -1:        /* Done with options. */
                break;
            
            default:        /* Something else: unexpected. */
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
    
    num_random_graphs = 0;
   
    subgraph = new int*[subgraphSize];
    for (int i = 0; i < subgraphSize; i++)
        subgraph[i] = new int[subgraphSize+3];   //add two cells on subgraph save the number of nodes, edges

    num = 0;
    printf("Graphlets Size: %d\n", subgraphSize);
    printf("Input Graph: %s\n", input_filename);

    if (!ReadData(input_filename, Dpath))
        return 0;

    g->setPath(output_directory);

    double realStart_t = realtime();

    clock_t startTime = clock();
    gettimeofday(&tv, 0);
    //for main graph
    isRand = false;
    subgraphCounter = 0;
    Enumerate();
    

    gettimeofday(&tvPartition, 0);

    cout<<"partitionTime"<<(tvPartition.tv_sec - tv.tv_sec + (double)(tvPartition.tv_usec - tv.tv_usec) / CLOCKS_PER_SEC)<<endl;
    exit(1);

    int deleted=0;
    vector<std::string> frequentG;
    
    vector< pair<std::string, long long int> > finalTree;
    hash_map<std::string, vector< const std::string* > > degreeSeq;
    hash_map<std::string, long long int> degreeSeqCount;
    /*******************Tree isomorphic*****************************/

    /*
    int treeSize = treeInt.size();
 
    int flag=0;
    
    int countTree =0;
  
    for (hash_map<std::string, long long int>::iterator it = treeInt.begin(); it != treeInt.end(); it++)
    {
        std::string tree1 = it->first;
        cout<<"Tree1 = "<<tree1<<endl;
        long long int tempCount = it->second;

        flag = 0;
    
        countTree++;
        
        for (int vi = 0; vi < finalTree.size(); vi++)
        {

            std::string tree2 = finalTree[vi].first;


            if(treeIso(tree1,tree2))
            {
                finalTree[vi].second = tempCount + finalTree[vi].second;
                flag = 1;
                break;
            }
        }
        if(flag == 0)
        {
            pair<string, long long int> tempPair = make_pair(tree1, tempCount);
            finalTree.push_back(tempPair);
        }
        
    }
    */


 
    /******************write non-isomorphic trees into file***********/
    int frequentTree = 0;
    ofstream fileTree;

    fileTree.open(("resultTree_"+intToString(subgraphSize)+".txt").c_str());

    for (hash_map<std::string, long long int>::iterator it = treeInt.begin(); it != treeInt.end(); it++)
    {
        if(it->second > (THR*subgraphCounter))
        {
            frequentTree++;
            //fileTree<< it->first<<" "<<it->second<<"\n";
        }
    }
    
    fileTree.close();

    printf("OutTree: %d\n", frequentTree);

    treeInt.clear();
    finalTree.clear();
    gettimeofday(&treeover, 0);

    cout<<"Time of dealing with tree after Partition before graph: "<<(treeover.tv_sec - tvPartition.tv_sec + (double)(treeover.tv_usec - tvPartition.tv_usec) / CLOCKS_PER_SEC)<<endl;
    
    /****************************************************************/
    /*****************graph isomorphic********************************/
    int graphSize = graphInt.size();
    //cout<<"graphSize:"<<graphSize<<endl;
    vector<bool> CountedGraph(graphSize, true);

    int flag = 0;
    int flag2 = 0;
    //degree sequence
    /*1210 2013**/
    for (hash_map<std::string, long long int>::iterator it = graphInt.begin(); it != graphInt.end(); it++)
    {
        const std::string *graph1 = &(it->first);

        //cout<<"Graph1 = " << *graph1 <<endl;

        std::string gds1 = graphDegreeSequence(*graph1,subgraphSize,AsciiToInt);
        long long int tempCount = it->second;
        //vector<const std::string*> tempGraphs;
        //tempGraphs.push_back(graph1);

        if(degreeSeqCount.find(gds1)==degreeSeqCount.end())
        {
            degreeSeqCount[gds1] = tempCount;
            vector<const std::string*> tempGraphsTemp;
            tempGraphsTemp.push_back(graph1);
            degreeSeq[gds1] = tempGraphsTemp;
        }
        else
        {
            degreeSeqCount.find(gds1)->second = degreeSeqCount.find(gds1)->second + tempCount;
            degreeSeq.find(gds1)->second.push_back(graph1);
        }

    }

    gettimeofday(&DSeqover, 0);

    cout<<"Time of dealing with Degree sequence: "<<(DSeqover.tv_sec - treeover.tv_sec  + (double)(DSeqover.tv_usec - treeover.tv_usec ) / CLOCKS_PER_SEC)<<endl;

    hash_map<std::string, long long int> finalGraph;
    int countCam = 0;
    for (hash_map<std::string, vector<const std::string*> >::iterator it = degreeSeq.begin(); it != degreeSeq.end(); it++)
    {
        vector<const std::string*> tempV = it->second;//graph strings' address with the same degreeSeq 

        if (degreeSeqCount.find(it->first)->second <= (seqTHR*subgraphCounter))
        {
            continue;
        }

        for (int iv = 0; iv < tempV.size(); iv++)
        {
            string tempCam = *(tempV[iv]);
            //cout<<"tempCam = "<<tempCam<<endl;
            countCam++;
            string cam =  calculateCam(tempCam,subgraphSize,AsciiToInt);
            //cout<<"cam = "<< cam<<endl;

            hash_map<std::string, long long int>::iterator it2 =  finalGraph.find(cam);
            if (it2 == finalGraph.end())
            {

                finalGraph[cam] = graphInt.find(tempCam)->second;
            }
            else
            {
                (it2->second) =  (it2->second) + graphInt.find(tempCam)->second;
            }
        }
        tempV.clear();
    }


    gettimeofday(&Camover, 0);

    cout<<"Time of dealing with Cam: "<<(Camover.tv_sec - DSeqover.tv_sec   + (double)(Camover.tv_usec - DSeqover.tv_usec) / CLOCKS_PER_SEC)<<endl;

    // cout<<"here 3333"<<endl;
    ofstream fileGraph;
    fileGraph.open(("resultGraph_"+intToString(subgraphSize)+".txt").c_str());
    //cout<<"Number of finalGraph : "<<finalGraph.size()<<endl;
    int frequentGraph = 0;
    for (hash_map<std::string, long long int> ::iterator it = finalGraph.begin(); it != finalGraph.end(); it++)
    {
        if (it->second > THR*subgraphCounter)
        {
            frequentGraph ++;
            //fileGraph<< it->first<<" "<<it->second<<"\n";
        }
    }
    printf("OutGraph: %d\n", frequentGraph);
    fileGraph.close();

    /****************************************************************/
    //g->AllocateCounter();
    //std::cout<<"countCam"<<countCam<<std::endl;
    //std::cout<<"degreeSeq:"<<degreeSeq.size()<<std::endl;
    //std::cout<<"finalGraph:"<<finalGraph.size()<<std::endl;
    //std::cout<<"deleted:"<<deleted<<std::endl;
    //std::cout<<"frequentGCount:"<<frequentG.size()<<std::endl;
    
    //g->Extract();

    clock_t end_main_time = clock();
    main_time = difftime(end_main_time, startTime)/(double)CLOCKS_PER_SEC;
    //printf("Time Used for main graph: %f\n", main_time);

    subgraphCounterMain = subgraphCounter;
    //for random graphs
    srand(time(NULL));
    isRand = true;
    long long boz = 0;
    //printf("Number of Random Graphs: %d\n", num_random_graphs);
    for (i = 1; i <= num_random_graphs; i++)
    {
        
        gen.genRandGraph_Edge(g);
        subgraphCounter = 0;
        Enumerate();
        g->Extract();
    }

    if (0 < num_random_graphs)
        g->calculateZSCORE(num_random_graphs, subgraphCounterMain, output_directory);

    for(i = 0; i < subgraphSize; i++)
    {
        delete [] Index[i];
        delete [] childSet[i];
    }

    for (int i = 0; i < subgraphSize; i++)
        delete [] subgraph[i];

    delete [] subgraph;
    delete [] Index;
    delete [] childSet;
    delete [] Visited;
    delete g;
    gettimeofday(&tv1, 0);
    clock_t endTime = clock();
    double total_time = difftime(endTime, startTime)/(double)CLOCKS_PER_SEC;
    printf("Time Used: %f\n", total_time);
    cout<<(tv1.tv_sec - tv.tv_sec + (double)(tv1.tv_usec - tv.tv_usec) / CLOCKS_PER_SEC)<<endl;

    double realEnd_t = realtime();
    cout<<"Total real Time: "<<realEnd_t - realStart_t <<endl;
    printf("=============================\n");
    //printf("graphInt: %d\n", graphInt.size());

    //printf("finalGraph: %d\n", finalGraph.size());
    printf("Tree Cam Count: %d\n", frequentTree);
    printf("Graph Cam Count: %d\n", finalGraph.size());
    printf("Total Cam Count: %d\n", finalGraph.size() + frequentTree);
    printf("Total Number of Subgraphs: %lld\n", subgraphCounter);


    return 0;
}
