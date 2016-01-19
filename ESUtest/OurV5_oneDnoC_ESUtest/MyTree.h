
#ifndef MYTREE_H
#define MYTREE_H
#include <list>
#include <iostream>
#include <stdint.h>
#include <ext/hash_map>
#include <string>
using namespace   __gnu_cxx;
namespace   __gnu_cxx
{
        template<> struct hash< std::string >
        {
                size_t operator()( const std::string& x ) const
                {
                        return hash< const char* >()( x.c_str() );
                }
        };
}

typedef  uint64_t uintx_t;
extern int* treeChildrenSize;
extern std::string* BFSVec;

/* 邻接表存储有向图 */
class MyTree
{
  int V;            // 顶点的数量
  std::vector<int> *adj;       // 邻接表
public:
  MyTree(int V);        // 构造函数
  ~MyTree();
  void addEdge(int v, int w); // 向图中添加一条边
  //void BFS(int v);       // BFS遍历
  std::string BFString(int v, bool isD);      
  std::string BFString(int v);
  void PrintTree(int);
  std::string TreeCamGen(int** subgraphTree, int level);
  void TreeAdjust( int parentNode);
  int ChildrenNum(int v);

};

#endif //MYTREE_H