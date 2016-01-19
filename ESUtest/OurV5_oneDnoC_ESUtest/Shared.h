#ifndef SHARED_H_
#define SHARED_H_

inline 		int 	max ( unsigned int a, unsigned int b ){
		return a > b ? a : b;
}

enum bfsState {FRESH, OPEN, CLOSED};

#endif /*SHARED_H_*/
