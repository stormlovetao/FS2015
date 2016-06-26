#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include <stdlib.h>
int main(int argcc, char *argv[])
{
	const char* const short_options = "h:i:";
	const struct option long_options[] = 
	{
		{"help", 1, NULL, 'h'},
		{"input",1, NULL, 'i'}	
	};
	int test;
	int next_option;
	do{
		next_option = getopt_long(argcc, argv, short_options, long_options, NULL);
		switch(next_option)
		{
			case 'h':
				test = atoi(optarg);
				break;
			case -1:
				break;
		}

	}while(next_option != -1);
	printf("%d", test);
	return 0;
}
