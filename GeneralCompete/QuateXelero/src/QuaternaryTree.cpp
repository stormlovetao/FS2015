#include "QuaternaryTree.h"

#define MAX 1000

void QTree::Check()
{
		if(cur_node->leaf != NULL)
			if(!cur_node->leaf->tag)
				t->update_leaf(cur_node->leaf, 1);
}

void QTree::allocate_node() {
	node = new QNode [MAX];
	node_ptr = 0;
}

QNode * QTree::add_node() {
	QNode * n = &node[node_ptr];
	node_ptr++;
	if(node_ptr == MAX)
		allocate_node();
	return n;
}