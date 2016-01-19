#ifndef	QUATERNARYTREE_H
#define QUATERNARYTREE_H

#include <string>
#include <math.h>
#include "ZeroOneTree.h"

using namespace std;


class QNode {
public:
	QNode *mone, *zero, *pone, *two;
	Leaf* leaf;
	QNode() {leaf = NULL; mone = zero = pone = two = NULL;}
	//~TerNode() {if(mone) delete mone; if(zero) delete zero; if(pone) delete pone;}
};

class QTree {
public:
	QNode* root;
	QNode * node;
	int node_ptr;
	QNode* cur_node;
	tree* t;

	QNode * add_node();
	void allocate_node();
	QTree() {root = cur_node = new QNode; t = new tree(1); allocate_node();}
	~QTree() { delete root; delete [] node; }
	void Check();
};
#endif //QUATERNARYTREE_H
