#include "randomGenerator.h"

generator::generator() {
	
}

int generator::binarySearch(int *E, int v, int w, int l, int h) {
	int m;
	while(h >= l) {
		m = (l+h)/2;
		if(E[m] == w)
			return m;
		else {
			if(E[m] < w)
				h = m - 1;
			else
				l = m + 1;
		}
	}
	return -1;
}

/////////////////// NOT USED ////////////////////
void generator::genRandGraph_Edge(Graph * g) {
	register int i, j;
	int len = g->Size();
	int p, q, a, b, c, d, c_ind, d_ind;
	int *Na, *Nb;
	
	for(j = 0; j < numOFexchange; j++) {
		for (a = 1; a <= len; a++) {
			Na = g->getNeighbours(a);
			do {
				b = g->get_vertex();	
			} while (a == b);
			
			Nb = g->getNeighbours(b);
			fflush(stdout);
			
			p = q = -1;
			
			for (i = 0; i < numOFexchange; i++) {
				c_ind = (rand() % Na[0]) + 1;
				if(g->isConnected(a, Na[c_ind]) && !g->isConnected(b, Na[c_ind]) && Na[c_ind] != b && !g->isConnected(Na[c_ind], a)) {
					c = Na[c_ind];
					p = c_ind;
					break;
				}
			}
			
			if (p == -1)
				continue;
			
			for (i = 0; i < numOFexchange; i++) {
				fflush(stdout);
				d_ind = (rand() % Nb[0]) + 1;
				if(g->isConnected(b, Nb[d_ind])	&& !g->isConnected(a, Nb[d_ind]) && Nb[d_ind] != a && !g->isConnected(Nb[d_ind], b)) {
					d = Nb[d_ind];
					q = d_ind;
					break;
				}
			}
			
			if (q == -1)
				continue;
			
			g->deleteEdgeAdjMat(a,c);
			g->deleteEdgeAdjMat(b,d);
			g->addEdgeAdjMat(a,d);
			g->addEdgeAdjMat(b,c);
			
			g->swapEdge(a, p, d);
			g->swapEdge(b, q, c);
		
			p = binarySearch(g->getNeighbours(d),d,b,1,g->getNeighbours(d)[0]);
			g->swapEdge(d, p, a);
			p = binarySearch(g->getNeighbours(c),c,a,1,g->getNeighbours(c)[0]);
			g->swapEdge(c, p, b);			
		}
	}
}

/////////////////// NOT USED ////////////////////
void generator::genRandGraph_Edge2(Graph * g) {
	register int i, j, k;
	int len = g->Size();
	int p, q, a, b, c, d, c_ind, d_ind;
	
	int *Na, *Nb, *outs;

	for(j = 0; j < numOFexchange; j++)
	{
		for (a = 1; a <= len; a++)
		{
			Na = g->getNeighbours(a);
			for(k = 0; k < Na[0]; k++)
			{
				p = k+1;
				c = Na[p];
				if(g->isConnected(c, a))
					continue;

				for (i = 0; i < numOFtrial; i++)
				{
					b = (rand() % len) + 1;
					if(a == b || c == b || g->isConnected(b, c))
						continue;

					Nb = g->getNeighbours(b);
					outs = new int[Nb[0]];
					int size = 0;
					for(d_ind = 1; d_ind <= Nb[0]; d_ind++)
						if(g->isConnected(b, Nb[d_ind]))
							outs[size++] = d_ind;

					if(size == 0)
						continue;

					d_ind = rand() % size;
					q = outs[d_ind];
					d = Nb[q];
					if(g->isConnected(d, b))
						continue;
					if(a == d || c == d || g->isConnected(a, d))
						continue;

					break;
				}
				if(i < numOFtrial)
				{
					g->deleteEdgeAdjMat(a,c);
					g->deleteEdgeAdjMat(b,d);
					g->addEdgeAdjMat(a,d);
					g->addEdgeAdjMat(b,c);

					g->swapEdge(a, p, d);
					g->swapEdge(b, q, c);

					p = binarySearch(g->getNeighbours(d),d,b,1,g->getNeighbours(d)[0]);
					g->swapEdge(d, p, a);
					p = binarySearch(g->getNeighbours(c),c,a,1,g->getNeighbours(c)[0]);
					g->swapEdge(c, p, b);			

				}//if
			}//k
		}//a
	}//j
}

/////////////////// USED ///////////////////
void generator::genRandGraph_Edge3(Graph * g) {
	register int i, j, k;
	int len = g->Size();
	int p, q, a, b, c, d, s, d_ind;
	
	int *Na, *tNa, *Nb, *outs;
	Na = new int[g->MaxDegree()+1];
	bool bidir;


	for(j = 0; j < numOFexchange; j++)
	{
		for (a = 1; a <= len; a++)
		{
			tNa = g->getNeighbours(a);
			s = 0;
			for(k = 0; k < tNa[0]; k++)
				if(g->isConnected(a, tNa[k+1]))
					Na[s++] = tNa[k+1];

			for(k = 0; k < s; k++)
			{
				c = Na[k];

				bidir = false;
				if(g->isConnected(c, a))
					bidir = true;

				for (i = 0; i < numOFtrial; i++)
				{
					b = (rand() % len) + 1;
					if(a == b || c == b || g->isConnected(b, c) || g->isConnected(c, b))
						continue;

					Nb = g->getNeighbours(b);
					outs = new int[Nb[0]];
					int size = 0;
					for(d_ind = 1; d_ind <= Nb[0]; d_ind++)
						if(g->isConnected(b, Nb[d_ind]))
							outs[size++] = d_ind;

					if(size == 0)
					{
						delete outs;
						continue;
					}

					d_ind = rand() % size;
					q = outs[d_ind];
					d = Nb[q];
					if((!bidir && g->isConnected(d, b)) || (bidir && !g->isConnected(d, b)))
					{
						delete outs;
						continue;
					}
					if(a == d || c == d || g->isConnected(a, d) || g->isConnected(d, a))
					{
						delete outs;
						continue;
					}

					break;
				}
				if(i < numOFtrial)
				{

					g->deleteEdgeAdjMat(a,c);
					g->deleteEdgeAdjMat(b,d);
					g->addEdgeAdjMat(a,d);
					g->addEdgeAdjMat(b,c);

					if(bidir)
					{
						g->deleteEdgeAdjMat(c,a);
						g->deleteEdgeAdjMat(d,b);
						g->addEdgeAdjMat(d,a);
						g->addEdgeAdjMat(c,b);
					}

					p = binarySearch(g->getNeighbours(a),a,c,1,g->getNeighbours(a)[0]);
					if(p==-1)
					{
						printf("-1!!!\n");
						int t;
						scanf("%d", &t);
					}
					g->swapEdge(a, p, d);
					g->swapEdge(b, q, c);
					p = binarySearch(g->getNeighbours(d),d,b,1,g->getNeighbours(d)[0]);
					if(p==-1)
					{
						printf("-1!!!\n");
						int t;
						scanf("%d", &t);
					}
					g->swapEdge(d, p, a);
					p = binarySearch(g->getNeighbours(c),c,a,1,g->getNeighbours(c)[0]);
					if(p==-1)
					{
						printf("-1!!!\n");
						int t;
						scanf("%d", &t);
					}
					g->swapEdge(c, p, b);
				}//if
			}//k
		}//a
	}//j
}
