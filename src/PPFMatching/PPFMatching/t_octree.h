
#ifndef _T_OCTREE_H_
#define _T_OCTREE_H_

#include <vector>
//using namespace std;

typedef struct TOctreeNode
{
	float center[3], halfDim[3];
	float* data;
	struct TOctreeNode *children[8];
}TOctreeNode;

/*

// designed for linear octree
typedef struct TOctree
{
	int maxLevels, numNodes;
	TOctreeNode* allNodes;
	TOctreeNode* root;
}TOctree;
*/

#if defined(__cplusplus)
extern "C" {
#endif 

	int t_octree_insert(TOctreeNode* octree, float* data);
	void t_octree_init(TOctreeNode* root, float center[3], float halfDim[3]);
	void t_octree_destroy(TOctreeNode *octree);
	void t_octree_query_in_bbox(TOctreeNode *octree, float xrange[2], float yrange[2], float zrange[2], std::vector<float*>& results);

	__inline static int get_octant(const float* point, const float root[3])
	{
		int oct = 0;
		oct |= 4 * (point[0] >= root[0]);
		oct |= 2 * (point[1] >= root[1]);
		oct |= 1 * (point[2] >= root[2]);
		return oct;
	}

#if defined(__cplusplus)
}
#endif 

#endif 