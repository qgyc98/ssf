#ifndef __MYPAIR_H
#define __MYPAIR_H

struct MyPair
{
	long node1;
	long node2;
	MyPair(){}
	~MyPair(){}


};

struct Less
{
		// porovnava se pouze podle jednoho bodu
	bool operator()(const MyPair & p1,const MyPair & p2) const { return  p1.node1 < p2.node1;  }
};

#endif


