#include <vector>
#include "MyTree.h"
#include <stdlib.h>

using namespace std;

int main(int argc, char const *argv[])
{
	MyTree tree(5);
	tree.addEdge(1,2);
	tree.addEdge(1,3);
	tree.addEdge(3,4);
	tree.addEdge(3,5);
	int** subgraphtree;
	subgraphtree = new int*[5];
	for (int i = 0; i < 5; ++i)
	{
		subgraphtree[i] = new int[5+1];
	}

	tree.PrintTree(5);
	subgraphtree[0][0]=1, subgraphtree[0][1]=1;
	subgraphtree[1][0]=2, subgraphtree[1][1]=2,subgraphtree[1][2]=3;
	subgraphtree[2][0]=2, subgraphtree[2][1]=4,subgraphtree[2][2]=5;

	string cam = tree.TreeCamGen(subgraphtree, 3);
	cout<<cam<<endl;
	return 0;
}