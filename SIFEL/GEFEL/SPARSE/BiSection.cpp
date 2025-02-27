// BiSection.cpp 

#include "BiSection.h"

DSS_NAMESPASE_BEGIN

CMcKee::CMcKee()
{
	p_node_level = p_order = nodes = NULL;
	size=domA=domB=0;
}

void CMcKee::Init(SparseConectivityMtxII* mtx)
{
	this->mtx = mtx;
	n=mtx->N();

	p_node_level = new long[n];
	memset(p_node_level,0,n*sizeof(long));

	p_order = new long[n];

	nodes = NULL; 
	size = 0;
}

CMcKee::~CMcKee()
{
	delete [] p_node_level;
	delete [] p_order;
}

BOOL CMcKee::IsAvailable(int v)
{
	return (p_node_level[v]==0);
}

void CMcKee::PrepareValid()
{
	for (long v=0; v<n; v++)
		p_node_level[v]=-1;

	for (long i=0; i<size;i++)
		p_node_level[nodes[i]]=0;
}

long CMcKee::FindFirstNode()
{
	long r=-1,cl,min = n+1;
	for (long i=0; i<size; i++)        
		if (IsAvailable(nodes[i]) && (cl=mtx->ColumnLength(nodes[i]))<min) {min = cl;r=nodes[i];}
	return r;
}


void CMcKee::ComputeLevels()
{
	PrepareValid();

	long l=1;
	long k=0;
	long last_level_start = 0;
	long next_level_start = 0;
	long last_level_count = 0;
	long next_level_count = 0;

	while(k<size)
	{
		if (last_level_count == 0)
		{
			long r=FindFirstNode();
			last_level_start = k;

			if (r==-1)// no available node 
				break;
			p_node_level[r] = l;
			p_order[k++] = r;last_level_count=1;
		}
		next_level_start = k;
		next_level_count = 0;
		// mark new level
		for (long ilr = last_level_start; ilr<next_level_start; ilr++)
		{
			long lr = p_order[ilr];
			//foreach (long neig in ColumnsIndexes[lr])
			IntArrayList& al= *mtx->ColumnsIndexes[lr];
			for (long idx = al.Count-1; idx>=0; idx--)
			{
				long neig = al.Items[idx];
				if (!IsAvailable(neig)) continue;
				p_node_level[neig] = l+1;
				p_order[k++] = neig;next_level_count++;
			}
		}
		l++;
		last_level_count = next_level_count;
		last_level_start = next_level_start;
	}
}

void CMcKee::DivideByMidLevel()
{
	long last_level_start = size/2;
	int midlevel = p_node_level[p_order[last_level_start]];
	while(last_level_start>0 && p_node_level[p_order[last_level_start-1]]==midlevel)
		last_level_start--;

	long next_level_start = last_level_start;

	while(next_level_start<size && p_node_level[p_order[next_level_start]]==midlevel)
		next_level_start++;

	domA = last_level_start;
	domB = size - next_level_start;

	Array::Reverse(p_order+last_level_start,size-last_level_start);
	Array::Copy(p_order,nodes,size);
}

CBiSection::CBiSection(SparseConectivityMtxII* mtx)
{
	this->mtx = mtx;
	mck.Init(mtx);
}

void CBiSection::RecurBiSectOrder(IntArrayList* order)
{
	RecurBiSect(order->Items,order->Count);
}

void CBiSection::RecurBiSect(long* nodes,long size)
{
	if (size<=3)
		return;

	long dA,dB;
	BiSect(nodes,size,dA,dB);
	RecurBiSect(nodes,dA);
	RecurBiSect(nodes+dA,dB);
}

void CBiSection::BiSect(long* nodes,long size, long& domA, long& domB)
{
	mck.nodes = nodes;
	mck.size = size;
	mck.ComputeLevels();
	mck.DivideByMidLevel();

	domA = mck.domA;
	domB = mck.domB;
}



DSS_NAMESPASE_END
