#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>
#include <stdint.h>

#include "MyTree.h"
using namespace std;
  
extern char* ascii;


uintx_t const_1 = 1ULL;
uintx_t const_2 = 2ULL;
uintx_t const_3 = 3ULL;

/***** 构造函数 *****/
MyTree::MyTree(int V)
{
  this->V = V;
  adj = new list<int>[V+1];   // 初始化V+1条链表
  
}

MyTree::~MyTree()
{
  // for (int i = 0; i<V+1; i++)
  // {
  //   delete adj[i];
  // }
  delete [] adj;
}
  
/* 添加边，构造邻接表 */
void MyTree::addEdge(int v, int w)
{
  adj[v].push_back(w);     // 将w加到v的list
}
  
/* 从顶点v出发广度优先搜索 */
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
}

string MyTree::BFString(int v, bool rootIsD )
{
  string re = "";
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
  
  list<int>::iterator i;
  
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


uintx_t MyTree::BFString(int v )
{
  
  uintx_t re = 0;
  
  int offset = sizeof(uintx_t)*8;

  list<int> queue;
  
  queue.push_back(v);
  
  //re.push_back('1');
  offset -= 2;
  re |= (const_1<<offset);
  offset -= 2;
  //re.push_back('<');
  re |= (const_2<<offset);
  
  list<int>::iterator i;
  
  while(!queue.empty())
  {
    // 出队
   int k = queue.front();
    
    queue.pop_front();

   
    for(i = adj[k].begin(); i!=adj[k].end(); i++)
    {
      queue.push_back(*i);
       
      //re.push_back('1');
      offset -= 2;
      if (offset < 0)
      {
        cout<< "In function MyTree::BFString, offset error!"<<endl;
        exit(1);
      }
      re |= (const_1<<offset);

      
    }
    //re.push_back('<');
    offset -= 2;
      if (offset < 0)
      {
        cout<< "In function MyTree::BFString, offset error!"<<endl;
        exit(1);
      }
    re |= (const_2<<offset);
  
  }

  while( int((re>>offset)&const_3) == int(const_2))
  {
    re &= (~const_3<<offset);
    offset += 2;
  }
  offset -= 2;
  re |= (const_3<<offset);

  
  return re;
} 

string MyTree::TreeCamGen(int** subgraphTree, int level)
{
  int i,j,k,curlevel,uplevel,parentNode;
  for(curlevel = level-2; curlevel >= 1; curlevel--)// we don't need to adjust leaf nodes in the lowest level
  {
    uplevel = curlevel-1;
    for(i = 1; i<=subgraphTree[uplevel][0]; i++)
    {
      parentNode = subgraphTree[uplevel][i];
      if (adj[parentNode].size() <= 1)
      {
        continue;
      }
      else
      {
        TreeAdjust(parentNode);
      } 
    }

  }
  int rootNode = subgraphTree[0][1];
  
  return BFString(rootNode, true);
}

bool sortCmp(pair<int, uintx_t> a, pair<int, uintx_t> b)
{
  return a.second < b.second;
}

void MyTree::TreeAdjust( int parentNode )
{
  list<int> nodes_list(adj[parentNode]);//copy elements of adj[parentNode] into nodes_list
  vector< pair<int, uintx_t> > tmp;
  for (list<int>::iterator iter = nodes_list.begin(); iter != nodes_list.end(); iter++)
  {
    uintx_t re = BFString(*iter);
    pair<int, uintx_t>apair = make_pair(*iter, re);
    tmp.push_back(apair);
  }

  sort(tmp.begin(), tmp.end(), sortCmp);
  list<int> list_tmp;
  for(vector< pair<int, uintx_t> >::iterator iter = tmp.begin(); iter!= tmp.end(); iter++)
  {
    list_tmp.push_back(iter->first);
  }
  adj[parentNode].swap(list_tmp);




}





void MyTree::PrintTree(int v)
{
  int i;
  std::list<int>::iterator iter;

  for (i = 1; i<= v;i++)
  {
    cout<<i<<": ";
    for (iter = adj[i].begin(); iter!=adj[i].end(); iter++)
    {
      cout<<*iter;
    }
   
    cout<<endl;
  }
}