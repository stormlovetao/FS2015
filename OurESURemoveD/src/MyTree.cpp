#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>
#include <stdint.h>

#include "MyTree.h"
using namespace std;

#define BFSLENGTH 64

extern int* treeChildrenSize;
extern string* BFSVec;
// int treeChildrenSize[] = {2,2,0,0,0};
// string BFSVec[] = {"","","","",""};


/***** 构造函数 *****/
MyTree::MyTree(int V)
{
  this->V = V;
  adj = new vector<int>[V+1];   // 初始化V+1条链表
}

MyTree::~MyTree()
{
  delete [] adj;
}
  
/* 添加边，构造邻接表 */
void MyTree::addEdge(int v, int w)
{
  adj[v].push_back(w);     // 将w加到v的list
}

int MyTree::ChildrenNum(int v)
{
  return adj[v].size();
}  


string MyTree::BFString(int v, bool rootIsD )
{
  string re = "";
  re.reserve(BFSLENGTH);
  list<int> queue;
  
  queue.push_back(v);
  if(rootIsD)
  {
	re.push_back('0');
  }
  else
  {
	re.push_back('1');
  }
	
  re.push_back('<');
  
  vector<int>::iterator i;
  
  while(!queue.empty())
  {
	// 出队
   int k = queue.front();
	
	queue.pop_front();

   
	for(i = adj[k].begin(); i!=adj[k].end(); i++)
	{
	  queue.push_back(*i);
	   
	  re.push_back('1');
	  
	}
	re.push_back('<');
  
  }
  while(re[re.size()-1] == '<')
  {
	re.resize(re.size()-1);
  }
  re.push_back('>');
  return re;
}

string MyTree::BFString(int v )
{
  string re = "";
  re.reserve(BFSLENGTH);
  list<int> queue;
  
  queue.push_back(v);
  
  re.push_back('1');
  
	
  re.push_back('<');
  
  vector<int>::iterator i;
  
  while(!queue.empty())
  {
	// 出队
   int k = queue.front();
	
	queue.pop_front();

   
	for(i = adj[k].begin(); i!=adj[k].end(); i++)
	{
	  queue.push_back(*i);
	   
	  re.push_back('1');
	  
	}
	re.push_back('<');
  
  }
  while(re[re.size()-1] == '<')
  {
	re.resize(re.size()-1);
  }
  re.push_back('>');
  return re;
}


bool sortCmp1(int a, int b )
{
  return treeChildrenSize[a] > treeChildrenSize[b];
}

bool sortCmp2(int a, int b)
{
  return BFSVec[a] < BFSVec[b];
}
string MyTree::TreeCamGen(int subgraphSize)
{
	
	// cout<<"Entering TreeCamGen, subgraphSize= "<<subgraphSize<<endl; 
	// PrintTree(5);
	int** subgraphtree;
	subgraphtree = new int*[subgraphSize];
	for (int i = 0; i < subgraphSize; ++i)
	{
		subgraphtree[i] = new int[subgraphSize+1];
	}
	// root node is labeled as 1
	int root = 1;
	
	subgraphtree[0][0] = 1;// first level only contain the root, size = 1
	subgraphtree[0][1] = root;
	int level = 1;
	int nextLevelCount = adj[root].size();
	int* nextLevelNodes =  new int[subgraphSize];
	int ind = 0;
	for (vector<int>::iterator it = adj[root].begin(); it != adj[root].end(); it++)
	{
		nextLevelNodes[ind] = *it;
		ind++;
	}
	
	// for (int i = 0; i < nextLevelCount; ++i)
	// {
	// 	cout<<nextLevelNodes[i]<<" ";
	// }
	// cout<<endl;
	while(nextLevelCount != 0)
	{
		subgraphtree[level][0] = nextLevelCount;
		for (int i = 0; i < nextLevelCount; ++i)
		{
			subgraphtree[level][i+1] = nextLevelNodes[i]; // store nodes on LEVEL level+1, level starts from 0
		}
		
		//update nextLevelNodes
		nextLevelCount = 0;
		for (int i = 0; i < subgraphtree[level][0]; ++i)
		{
			int tmpNode = subgraphtree[level][i+1];

			for (vector<int>::iterator it = adj[tmpNode].begin(); it != adj[tmpNode].end(); it++)
			{
				nextLevelNodes[nextLevelCount] = *it;
				nextLevelCount++;
			}
		}
		level++;
		
	}
	//cout<<"level = "<<level<<endl;
	/**Tao add above on Aug14 2016**/



	int i,j,k,curlevel,uplevel,parentNode;
	for(curlevel = level-2; curlevel >= 1; curlevel--)
	{
		uplevel = curlevel-1;
		for(i = 1; i<=subgraphtree[uplevel][0]; i++)
		{
		  	parentNode = subgraphtree[uplevel][i];
			if (adj[parentNode].size() <= 1)
			{
				continue;
			}

			vector<int> childrenNumList;
			int flag = 1;
			for(vector<int>::iterator j = adj[parentNode].begin(); j != adj[parentNode].end(); j++)
			{
				int children_num = adj[*j].size();
				if (find(childrenNumList.begin(), childrenNumList.end(), children_num) != childrenNumList.end())
				{
					flag = 0;
					break;
				}
				else
					childrenNumList.push_back(children_num);
			}
		  
			if (flag==1)
			{ 
				sort(adj[parentNode].begin(), adj[parentNode].end(), sortCmp1);
				continue;
			}
			else
			{

				//cout<< "dddd "<<parentNode<<endl;
				TreeAdjust(parentNode);
			} 
		}

	}
	int rootNode = subgraphtree[0][1];
	string outstr = BFString(rootNode);
	for (int i = 0; i < subgraphSize; ++i)
	{
		delete [] subgraphtree[i];
	}

  	return outstr;
}

string MyTree::TreeCamGen(int** subgraphTree, int level)
{
  int i,j,k,curlevel,uplevel,parentNode;
  for(curlevel = level-2; curlevel >= 1; curlevel--)
  {
	uplevel = curlevel-1;
	for(i = 1; i<=subgraphTree[uplevel][0]; i++)
	{
	  parentNode = subgraphTree[uplevel][i];
	  if (adj[parentNode].size() <= 1)
	  {
		continue;
	  }

	  vector<int> childrenNumList;
	  int flag = 1;
	  for(vector<int>::iterator j = adj[parentNode].begin(); j != adj[parentNode].end(); j++)
	  {
		int children_num = adj[*j].size();
		if (find(childrenNumList.begin(), childrenNumList.end(), children_num) != childrenNumList.end())
		{
		   flag = 0;
		   break;
		}
		else
		  childrenNumList.push_back(children_num);
	  }
	  
	  if (flag==1)
	  { 

		sort(adj[parentNode].begin(), adj[parentNode].end(), sortCmp1);
		continue;
	  }
	  else
	  {
		//cout<< "dddd"<<parentNode<<endl;
		TreeAdjust(parentNode);
	  } 
	}

  }
  int rootNode = subgraphTree[0][1];
  return BFString(rootNode, true);
}



void MyTree::TreeAdjust(int parentNode)
{
  

  for (vector<int>::iterator iter = adj[parentNode].begin(); iter != adj[parentNode].end(); iter++)
  {
	BFSVec[*iter] = BFString(*iter);
  }

  sort(adj[parentNode].begin(), adj[parentNode].end(), sortCmp2);
  
}
int MyTree::TreeNodeNum()
{
	return this->V;
} 

void MyTree::PrintTree(int v)
{
  int i;
  std::vector<int>::iterator iter;

  for (i = 1; i<= v;i++)
  {
	cout<<i<<": ";
	for (iter = adj[i].begin(); iter!=adj[i].end(); iter++)
	{
	  cout<<*iter<<" ";
	}
   
	cout<<endl;
  }
}