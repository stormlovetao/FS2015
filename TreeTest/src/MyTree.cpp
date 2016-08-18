#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>
#include <stdint.h>

#include "MyTree.h"
using namespace std;

#define BFSLENGTH 64
int treeChildrenSize[]={2,0,2,0,0};
std::string BFSVec[]={"","","","",""};

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
/* 从顶点v出发广度优先搜索 */
/*
void MyTree::BFS(int v)
{
  // BFS辅助队列
  list<int> queue;
  
  
  queue.push_back(v);
  
  list<int>::iterator i;
  
  while(!queue.empty())
  {
	// 出队
	v = queue.front();
	cout << v << " ";
	queue.pop_front();
  
	// 检测已出队的顶点s的所有邻接顶点
	
	for(i = adj[v].begin(); i!=adj[v].end(); i++)
	{
	  
		queue.push_back(*i);
	   
	}
  }
}*/



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
  return BFSVec[a] > BFSVec[b];
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
  return BFString(rootNode);
}



void MyTree::TreeAdjust(int parentNode)
{
  

  for (vector<int>::iterator iter = adj[parentNode].begin(); iter != adj[parentNode].end(); iter++)
  {
	BFSVec[*iter] = BFString(*iter);
  }

  sort(adj[parentNode].begin(), adj[parentNode].end(), sortCmp2);
  
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