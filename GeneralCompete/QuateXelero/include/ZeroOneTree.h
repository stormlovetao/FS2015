#ifndef	ZEROONETREE_H
#define ZEROONETREE_H

#include <string>
#include <math.h>

using namespace std;


class Node {
public:
	bool visited;
	Node * left;
	Node * right;
	Node() {	left = NULL; right = NULL; visited = false;} //printf("node counter: %d\n", ++otherNodeCounter);}
	//~Node() { /*printf("Destroying ... \n"); fflush(stdout); }*/ if(right != NULL) delete right; if(left != NULL) delete left; }
};

class Leaf:public Node {
	public:
		long long unsigned count;
		bool tag; //indicates whether this leaf is a real leaf or an internal node in the binary tree (false = leaf, true = internal node)
		Leaf() {count = 0; tag = false;}
		Leaf(bool t) {tag = t;}
		//~Leaf() {delete subgraphIndexes;}
};

struct value {
	double * C;
	char ** adj_str;
};

class tree {
	protected:
		Node * root;
		Node * node;
		Node * cur_node;
		int node_ptr;
		int leaf_ptr;
		unsigned long leaf_num;
		double * C;
		int head;
		value v;
		int subgraphSize;
	public:
		Leaf * leaf;
		Leaf * trueLeaf;////
		tree(int);
		Node * add_node();
		Node * return_root();
		Leaf * add_leaf();
		void update_leaf(Leaf * l, int val);
		void allocate_node();
		void allocate_leaf();
		unsigned long  get_leafnum();
		void init_cur_node();
		void insert_one_main();
		void insert_zero_main();
		bool insert_one_rand();
		bool insert_zero_rand();
		Leaf* update_one_main(int val);
		Leaf* update_zero_main(int val);
		Leaf* update_one_rand(int val);
		Leaf* update_zero_rand(int val);
		double * extract();
		void DFS(Node * cur);
		value destroy();
		void DFS_value(Node * cur, char * str, int lev);
        ~tree();
};
#endif //ZEROONETREE_H
