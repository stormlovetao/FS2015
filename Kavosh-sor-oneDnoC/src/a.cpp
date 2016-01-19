int* GetPtn(int nV, int Dnum){

	int tempGraph[nV];
	for(int i = 0; i <nV; i++)
	{
		tempGraph[i] = i+1;
	}

	int m = (nV + WORDSIZE -1)/WORDSIZE;

	graph * g = new graph[nV * MAXM];

	set *gv;
	int i,j;
	for( i = 0; i < nV; i++)
	{
		gv = GRAPHROW(g, i, m);
		EMPTYSET(gv, m);
		for(j = 0; j < nV; j++)
		{
			if(i == j)
				continue;
			if(isConnected(tempGraph[i],tempGraph[j]))
			{
				ADDELEMENT(gv, j);
			}
		}
	}
	nvector lab[nV],ptn[nV],orbits[nV];
	static DEFAULTOPTIONS(options);
    setword workspace[160*subgraphSize];
    /*init for nauty*/
    options.writeautoms = FALSE;
    options.writemarkers = FALSE;
    options.getcanon = TRUE;
    options.defaultptn = FALSE;/
    options.digraph = FALSE;
    statsblk(stats);

	for(i = 0; i< nV; i++)
		lab[i] = i;
	for(i = 0; i < Dnum-1; i++)
		ptn[i] = 1;
	ptn[Dnum-1] = 0;

	for(i = Dnum; i < nV-1; i++)
		ptn[i] = 1;
	ptn[nV -1] = 0;

	graph canong[nV*nV];
	nauty(g,lab,ptn,NULL,orbits,&options,&stats,workspace,160*MAXM,m,nV,canon);

	return ptn;
}