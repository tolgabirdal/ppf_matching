
#ifndef __T_PLY__H__
#define __T_PLY__H__

typedef struct
{
	float v[3][3];			// 3 Vertices of coordinates of triangle
	float n[3][3];			// corresponding normals
	unsigned int voxel;			// voxel index
} TFace3D; // this structure isn't used for now. Instead TTriangle3D is used

typedef struct
{
	int v1, v2, v3;		// indices into vertex array
	int n1, n2, n3;		// indices into normals array
	int t1, t2, t3;		// texture coordinate indices
	int c1, c2, c3;		// Indices into the colors array. For each of 3 vertices there is an index into color array
	int texID;				// opengl ID (used in opengl texture mapping)
	int groupID;			// group ID (in obj for example, for multiple objects)
} TTriangle3D;


/*  Here I am defining different structures for Vertex, Normal, TexCoords etc,
even though they look the same. This is because it's more friendly when coding
All the coordinates are in float precision. However, I would recommend doing any 
important computation in doubles, and casting the result to float if float is referring to
single precision.
*/

typedef struct
{
	float vx, vy, vz;
} TVertex3D;

typedef struct
{
	float nx, ny, nz;
} TNormal3D;

typedef struct
{
	float u, v, t;
} TTexCoords3D;

typedef struct
{
	float r, g, b;
} TRGBColor;

typedef struct TMesh_
{
	//int magic;
	TTriangle3D* triangles;		// triangle array of a mesh
	TVertex3D* vertices;			// vertex array of the mesh
	TNormal3D* normals;			// normals array of the mesh
	TTexCoords3D* texCoords;		// texture array of the mesh
	TRGBColor* colors;				// color array of the mesh
	int NumVertices;			// total number of vertices
	int NumTriangles;			// total number of triangles
	int isNormalsComputed;		// set this to 0 when vertices are updated
	int isTextureSet;			// have we loaded a texture?
	int** AdjFaceList;
	int *CountAdjFaceList;
	int isAdjFaceListComputed;
	int numConnected;
} TMesh;


/*

magic: A value identifying that this structure is a mesh. NOT USED

Vertices: 3D Vertices that make up the mesh
Triangle: 3 vertex indices.   mesh.vertices[t.v1] is first vertex of the triangle and so on (t is the triangle)
Normals:  Normals per all triangles or vertices. 
Colors:   Color array for triangles. It's defined for all vertices.

A color can be sampled like this: 

For first vertex:
cv[0] = mesh->colors[mesh->triangles[j].c1].r;
cv[1] = mesh->colors[mesh->triangles[j].c1].g;
cv[2] = mesh->colors[mesh->triangles[j].c1].b;
cv[3] = 1.0;

For second vertex:
cv[0] = mesh->colors[mesh->triangles[j].c2].r;
cv[1] = mesh->colors[mesh->triangles[j].c2].g;
cv[2] = mesh->colors[mesh->triangles[j].c2].b;
cv[3] = 1.0;

... and so on.



Texture Coordinates: Are (usually 2D) coordinates of refering to image coordinates. 
For each 3D vertex, there exists an image coordinate.


isNormalsComputed:  Set to 0 whenever vertices or faces changes. Set to 1 whenever normals are computed. 
Don't change it if an operation doesn't alter the surface normals.

isTextureSet:		Set this to 0 whenever texture changes. Set to 1 whenever texture is loaded, or its
vertices are recomputed.

*/

typedef struct
{
	TMesh** meshPtr;
	long counterVertex, counterTri;
	int eol;
}PlyArg;

#if defined __cplusplus
extern "C" {
#endif

	int read_mesh_ply(TMesh** Mesh3D, const char* FileName);
	int t_write_mesh_ply(TMesh* mesh, const char* FileName);
	void destroy_mesh(TMesh** mesh);
	TMesh* clone_mesh(const TMesh* mesh);
	void transform_mesh(TMesh* Mesh, double Pose[16]);
	TMesh* transform_mesh_new(TMesh* Mesh, double Pose[16]);
	int get_mesh_vertices(TMesh* Mesh, unsigned char* verticesPC, int stride, int withNormals, int flipViewpoint=0);

#if defined __cplusplus
}
#endif

#endif