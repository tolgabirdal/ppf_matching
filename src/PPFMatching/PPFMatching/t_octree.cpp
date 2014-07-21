#include "t_octree.h"
#include <cstdlib>

using namespace std;

int t_octree_insert(TOctreeNode* octree, float* data)
{
	int octNew = -1;
	if (octree->children[0]==0)
	{
		if (!octree->data)
			octree->data = data;
		else
		{
			int i, octOld;
			for(i=0; i!=8; ++i)
			{
				octree->children[i] = (TOctreeNode*)calloc(1, sizeof(TOctreeNode));
				octree->children[i]->center[0] = octree->center[0] + octree->halfDim[0] * (i&4 ? .5f : -.5f);
				octree->children[i]->center[1] = octree->center[1] + octree->halfDim[1] * (i&2 ? .5f : -.5f);
				octree->children[i]->center[2] = octree->center[2] + octree->halfDim[2] * (i&1 ? .5f : -.5f);

				octree->children[i]->halfDim[0] = octree->halfDim[0]*0.5f;
				octree->children[i]->halfDim[1] = octree->halfDim[1]*0.5f;
				octree->children[i]->halfDim[2] = octree->halfDim[2]*0.5f;
			}

			octOld = get_octant(octree->data, octree->center);
			octNew = get_octant(data, octree->center);

			t_octree_insert(octree->children[octOld], octree->data);
			t_octree_insert(octree->children[octNew], data);
		}
	}
	else
	{
		octNew = get_octant(data, octree->center);
		t_octree_insert(octree->children[octNew], data);
	}

	return octNew;
}

void t_octree_init(TOctreeNode *root, float center[3], float halfDim[3])
{
	int i;
	root->data=0;
	
	for (i=0; i!=8; i++)
		root->children[i]=0;

	root->center[0]=center[0];
	root->center[1]=center[1];
	root->center[2]=center[2];
	root->halfDim[0]=halfDim[0];
	root->halfDim[1]=halfDim[1];
	root->halfDim[2]=halfDim[2];
}

void t_octree_destroy(TOctreeNode *octree)
{
	if (octree)
	{
		if (octree->children)
		{
			int i;
			for(i=0; i!=8; ++i) 
			{
				t_octree_destroy(octree->children[i]);
				octree->children[i]=0;
			}			
		}
		octree->data=0;
		free(octree);
		octree=0;
	}
}

void t_octree_query_in_bbox(TOctreeNode *octree, float xrange[2], float yrange[2], float zrange[2], std::vector<float*>& results)
{
	if (octree->children[0]==0)
	{
		if (octree->data)
		{
			const float* p = octree->data;

			if(p[0]>xrange[1] || p[1]>yrange[1] || p[2]>zrange[1]) 
				return;
			if(p[0]<xrange[0] || p[1]<yrange[0] || p[2]<zrange[0]) 
				return;

			results.push_back(octree->data);
		}
	}
	else
	{
		int i;
		for(i=0; i!=8; ++i) 
		{
			float cmax[3], cmin[3];
			cmax[0] = octree->children[i]->center[0] + octree->children[i]->halfDim[0];
			cmax[1] = octree->children[i]->center[1] + octree->children[i]->halfDim[1];
			cmax[2] = octree->children[i]->center[2] + octree->children[i]->halfDim[2];

			if(cmax[0]<xrange[0] || cmax[1]<yrange[0] || cmax[2]<zrange[0]) 
				continue;

			cmin[0] = octree->children[i]->center[0] - octree->children[i]->halfDim[0];
			cmin[1] = octree->children[i]->center[1] - octree->children[i]->halfDim[1];
			cmin[2] = octree->children[i]->center[2] - octree->children[i]->halfDim[2];

			if(cmin[0]>xrange[1] || cmin[1]>yrange[1] || cmin[2]>zrange[1]) 
				continue;

			t_octree_query_in_bbox(octree->children[i], xrange, yrange, zrange, results);
		} 
	}
}