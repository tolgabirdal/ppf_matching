/*M///////////////////////////////////////////////////////////////////////////////////////
 //
 //  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
 //
 //  By downloading, copying, installing or using the software you agree to this license.
 //  If you do not agree to this license, do not download, install,
 //  copy or use the software.
 //
 //
 //                           License Agreement
 //                For Open Source Computer Vision Library
 //
 // Copyright (C) 2013, OpenCV Foundation, all rights reserved.
 // Third party copyrights are property of their respective owners.
 //
 // Redistribution and use in source and binary forms, with or without modification,
 // are permitted provided that the following conditions are met:
 //
 //   * Redistribution's of source code must retain the above copyright notice,
 //     this list of conditions and the following disclaimer.
 //
 //   * Redistribution's in binary form must reproduce the above copyright notice,
 //     this list of conditions and the following disclaimer in the documentation
 //     and/or other materials provided with the distribution.
 //
 //   * The name of the copyright holders may not be used to endorse or promote products
 //     derived from this software without specific prior written permission.
 //
 // This software is provided by the copyright holders and contributors "as is" and
 // any express or implied warranties, including, but not limited to, the implied
 // warranties of merchantability and fitness for a particular purpose are disclaimed.
 // In no event shall the Intel Corporation or contributors be liable for any direct,
 // indirect, incidental, special, exemplary, or consequential damages
 // (including, but not limited to, procurement of substitute goods or services;
 // loss of use, data, or profits; or business interruption) however caused
 // and on any theory of liability, whether in contract, strict liability,
 // or tort (including negligence or otherwise) arising in any way out of
 // the use of this software, even if advised of the possibility of such damage.
 //
 //M*/

// Author : Tolga Birdal. Please see the notice in tmesh.h

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <string.h> 
#include "c_utils.h"
#include "tmesh.h"
#include "rply.h"

#if defined T_OPENMP
#include <omp.h>
#endif

static int vertex_cb(p_ply_argument argument) 
{
	PlyArg* paPtr=0;
	long counterVertex=0, counterTri=0;
	TMesh** meshPtr=0;
	TMesh* mesh;
    long eol=0;
	double val;

	ply_get_argument_user_data(argument, (void**)(&paPtr), &eol);
    
	mesh=*((TMesh**)paPtr->meshPtr);

	counterVertex=(paPtr->counterVertex);

	val=ply_get_argument_value(argument);
	switch (eol)
	{
	case 0:
		mesh->vertices[counterVertex].vx=val;
		break;
	case 1:
		mesh->vertices[counterVertex].vy=val;
		break;
	case 2:
		mesh->vertices[counterVertex].vz=val;
		break;
	case 3:
		mesh->normals[counterVertex].nx=val;
		break;
	case 4:
		mesh->normals[counterVertex].ny=val;
		break;
	case 5:
		mesh->normals[counterVertex].nz=val;
	case 6:
		mesh->colors[counterVertex].r=(float)val/255.0f;
		break;
	case 7:
		mesh->colors[counterVertex].g=(float)val/255.0f;
		break;
	case 8:
		mesh->colors[counterVertex].b=(float)val/255.0f;
		break;
	};

    //if (eol == argument->length*3-1)
	if (eol == paPtr->eol)
		paPtr->counterVertex++;

    return 1;
}

static int face_cb(p_ply_argument argument) 
{
	PlyArg* paPtr=0;
	long counterVertex=0, counterTri=0;
	TMesh** meshPtr=0;
	TMesh* mesh;
    long length, value_index,eol=0;
	double val;
	ply_get_argument_user_data(argument, (void**)(&paPtr), &eol);

	mesh=*((TMesh**)paPtr->meshPtr);

    ply_get_argument_property(argument, NULL, &length, &value_index);
	val= ply_get_argument_value(argument);
	switch (value_index) {
		case -1:
			return 1;
			break;
		case 0:
			mesh->triangles[paPtr->counterTri].v1=val;
			mesh->triangles[paPtr->counterTri].c1=val;
			mesh->triangles[paPtr->counterTri].n1=val;
			break;
		case 1: 
			mesh->triangles[paPtr->counterTri].v2=val;
			mesh->triangles[paPtr->counterTri].c2=val;
			mesh->triangles[paPtr->counterTri].n2=val;
			break;
		case 2:
			mesh->triangles[paPtr->counterTri].v3=val;
			mesh->triangles[paPtr->counterTri].c3=val;
			mesh->triangles[paPtr->counterTri].n3=val;
			paPtr->counterTri++;
			break;
		default: 
			break;
	}

    return 1;
}


void create_mesh_3D(TMesh** Mesh, int NumVertices, int NumTriangles)
{
	TMesh** mesh=(TMesh**)Mesh;
	(*mesh)= (TMesh*)calloc(1, sizeof(TMesh)); //new TMesh();
	(*mesh)->vertices=0;
	(*mesh)->texCoords=0;
	(*mesh)->triangles=0;
	(*mesh)->normals=0;
	(*mesh)->colors=0;

	if (NumVertices)
	{
		(*mesh)->vertices= (TVertex3D*)calloc(NumVertices, sizeof(TVertex3D)); //new TVertex3D[NumVertices];
		(*mesh)->texCoords= (TTexCoords3D*)calloc(NumVertices, sizeof(TTexCoords3D)); //new TexCoords[NumVertices];
	}
	if (NumTriangles)
	{
		(*mesh)->triangles= (TTriangle3D*)calloc(NumTriangles, sizeof(TTriangle3D)); //new TTriangle3D[NumTriangles];
		(*mesh)->colors= (TRGBColor*)calloc(NumTriangles, sizeof(TRGBColor)); //new TRGBColor[NumTriangles];
	}

	(*mesh)->normals= (TNormal3D*)calloc(MAX(NumTriangles, NumVertices), sizeof(TNormal3D)); //new TNormal3D[NumTriangles];
	
	(*mesh)->NumTriangles=NumTriangles;
	(*mesh)->NumVertices=NumVertices;
	(*mesh)->isNormalsComputed=0;
	(*mesh)->isTextureSet=0;
	(*mesh)->numConnected=0;
}

int read_mesh_ply(TMesh** Mesh3D, const char* FileName) 
{
    long nvertices, ntriangles;
	PlyArg pa;
	long counterVertex=0, counterTri=0;
	int nelements;

	//TMesh* mesh=0;
	
    p_ply ply = ply_open(FileName, NULL, 0, NULL);

    if (!ply) 
	{
		printf("Cannot open PLY file\n");
		return -1;
	}
    
	if (!ply_read_header(ply)) 
	{
		printf("Cannot read ply header. Possibly the input is not in PLY format\n");
		return -2;
	}

	pa.counterTri=0;
	pa.counterVertex=0;
	pa.meshPtr=&(*Mesh3D);
    nvertices = ply_set_read_cb(ply, "vertex", "x", vertex_cb, ((void*)(&pa)), 0);
    ply_set_read_cb(ply, "vertex","y", vertex_cb, ((void*)(&pa)), 1);
    ply_set_read_cb(ply, "vertex", "z", vertex_cb, ((void*)(&pa)), 2);

	pa.eol = 3;
	nelements = ply_get_num_elements(ply);
	if (nelements>(long)1)
	{		
		ply_set_read_cb(ply, "vertex", "nx", vertex_cb, ((void*)(&pa)), 3);
		ply_set_read_cb(ply, "vertex", "ny", vertex_cb, ((void*)(&pa)), 4);
		ply_set_read_cb(ply, "vertex", "nz", vertex_cb, ((void*)(&pa)), 5);
		pa.eol = 5;
	}

	if (nelements>(long)2)
	{
		ply_set_read_cb(ply, "vertex", "red", vertex_cb, ((void*)(&pa)), 3);
		ply_set_read_cb(ply, "vertex", "green", vertex_cb, ((void*)(&pa)), 4);
		ply_set_read_cb(ply, "vertex", "blue", vertex_cb, ((void*)(&pa)), 5);
		pa.eol = 8;
	}

    ntriangles = ply_set_read_cb(ply, "face", "vertex_indices", face_cb, ((void*)(&pa)), 0);

	//create_mesh_3D(&(*Mesh3D), nvertices, ntriangles);
	create_mesh_3D(Mesh3D, nvertices, ntriangles);
	((TMesh*)(*Mesh3D))->numConnected=1;

	if (nelements>(long)1)
		(*Mesh3D)->isNormalsComputed = 1;

   // printf("%ld\n%ld\n", nvertices, ntriangles);
    if (!ply_read(ply)) 
	{
		ply_close(ply);
		return -3;
	}

    return 0;
}


int t_write_mesh_ply(TMesh* mesh, const char* FileName)
{
	int i,j;

	// Create .ply file in ascii format.
	p_ply oply = ply_create(FileName, PLY_ASCII, 0, 0, NULL);
	if (!oply) 
		return -1;

	// Add "vertex" element.
	if (!ply_add_element(oply, "vertex", mesh->NumVertices)) 
		return -1;

	// Add "vertex" properties. if the type parameter is not PLY_LIST, two last
	// parameters are ignored. So, PLY_LIST is passed for two last parameters
	// as the most unsuitable value.
	if (!ply_add_property(oply, "x", PLY_FLOAT, PLY_LIST, PLY_LIST)) 
		return -1;

	if (!ply_add_property(oply, "y", PLY_FLOAT, PLY_LIST, PLY_LIST)) 
		return -1;

	if (!ply_add_property(oply, "z", PLY_FLOAT, PLY_LIST, PLY_LIST)) 
		return -1;

	if (mesh->isNormalsComputed)
	{
		if (!ply_add_property(oply, "nx", PLY_FLOAT, PLY_LIST, PLY_LIST)) 
			return -1;

		if (!ply_add_property(oply, "ny", PLY_FLOAT, PLY_LIST, PLY_LIST)) 
			return -1;

		if (!ply_add_property(oply, "nz", PLY_FLOAT, PLY_LIST, PLY_LIST)) 
			return -1;
	}

	// Add "face" element.
	if (!ply_add_element(oply, "face", mesh->NumTriangles)) 
		return -1;
	// Add "face" only property. It is a list of vertex indices. 
	if (!ply_add_property(oply, "vertex_indices", PLY_LIST, PLY_UCHAR, PLY_UINT)) 
		return -1;

	// Add a comment and an obj_info.
	if (!ply_add_comment(oply, "This PLY File is generated by Tolga")) 
		return -1;
	if (!ply_add_obj_info(oply, "Tolga's PLY Mesh Dump")) 
		return -1;

	// Write .ply header.
	if (!ply_write_header(oply)) 
		return -1;

	if (mesh->isNormalsComputed)
	{
		for (i=0; i<mesh->NumVertices; i++)
		{
			ply_write(oply, mesh->vertices[i].vx);
			ply_write(oply, mesh->vertices[i].vy);
			ply_write(oply, mesh->vertices[i].vz);
			ply_write(oply, mesh->normals[i].nx);
			ply_write(oply, mesh->normals[i].ny);
			ply_write(oply, mesh->normals[i].nz);
		}
	}
	else
	{
		for (i=0; i<mesh->NumVertices; i++)
		{
			ply_write(oply, mesh->vertices[i].vx);
			ply_write(oply, mesh->vertices[i].vy);
			ply_write(oply, mesh->vertices[i].vz);
		}
	}

	for (i=0; i<mesh->NumTriangles; i++)
		{
			ply_write(oply, 3);
			ply_write(oply, mesh->triangles[i].v1);
			ply_write(oply, mesh->triangles[i].v2);
			ply_write(oply, mesh->triangles[i].v3);
		}

	if (!ply_close(oply))  
		return -1;

	return 0;
}

void transform_mesh(TMesh* Mesh, double Pose[16])
{
	const int N = Mesh->NumVertices;
	int i;

	double R[9], t[3];
	pose_to_rt(Pose, R, t);

#if defined T_OPENMP
#pragma omp parallel for
#endif
	for (i=0; i<N; i++)
	{
		TVertex3D v = (Mesh->vertices[i]);
		TNormal3D nr = (Mesh->normals[i]);
		double pth[4] = {(double)v.vx, (double)v.vy, (double)v.vz, 1};
		double vt[4];
		double n[3] = {(double)nr.nx, (double)nr.ny, (double)nr.nz}, n2[3];

		// transform vertices
		matrix_product441(Pose, pth, vt);

		if (fabs(vt[3])>EPS)
		{
			vt[0] = (float)(vt[0]/vt[3]);
			vt[1] = (float)(vt[1]/vt[3]);
			vt[2] = (float)(vt[2]/vt[3]);
		}

		Mesh->vertices[i].vx = vt[0];
		Mesh->vertices[i].vy = vt[1];
		Mesh->vertices[i].vz = vt[2];

		// rotate and normalize normals		
		matrix_product331(R, n, n2);
		double nNorm = sqrt(n2[0]*n2[0]+n2[1]*n2[1]+n2[2]*n2[2]);

		if (nNorm>EPS)
		{
			Mesh->normals[i].nx=(float)(n2[0]/nNorm);
			Mesh->normals[i].ny=(float)(n2[1]/nNorm);
			Mesh->normals[i].nz=(float)(n2[2]/nNorm);
		}
	}
}

TMesh* clone_mesh(const TMesh* mesh)
{
	TMesh* MeshClone=0;

	create_mesh_3D(&MeshClone, mesh->NumVertices, mesh->NumTriangles);
	
	memcpy(MeshClone->vertices,mesh->vertices,sizeof(TVertex3D)*mesh->NumVertices);
	memcpy(MeshClone->triangles,mesh->triangles,sizeof(TTriangle3D)*mesh->NumTriangles);
	memcpy(MeshClone->normals,mesh->normals,sizeof(TNormal3D)*mesh->NumTriangles);
	memcpy(MeshClone->colors,mesh->colors,sizeof(TRGBColor)*mesh->NumTriangles);
	memcpy(MeshClone->texCoords,mesh->texCoords,sizeof(TTexCoords3D)*mesh->NumVertices);
	MeshClone->isNormalsComputed=mesh->isNormalsComputed;
	MeshClone->isTextureSet=mesh->isTextureSet;
	MeshClone->numConnected=mesh->numConnected;

	return MeshClone;
}

TMesh* transform_mesh_new(TMesh* Mesh, double Pose[16])
{
	TMesh* newMesh = clone_mesh(Mesh);
	transform_mesh(newMesh, Pose);
	return newMesh;
}

void destroy_mesh(TMesh** mesh)
{
	if ((*mesh)->NumVertices)
	{
		memset((*mesh)->vertices,0,sizeof(TVertex3D)*(*mesh)->NumVertices);
		memset((*mesh)->texCoords,0,sizeof(TTexCoords3D)*(*mesh)->NumVertices);
		memset((*mesh)->normals,0,sizeof(TNormal3D)*(*mesh)->NumVertices);

		free((*mesh)->vertices);
		free((*mesh)->texCoords);
	}
	
	if ((*mesh)->NumTriangles)
	{
		memset((*mesh)->triangles,0,sizeof(TTriangle3D)*(*mesh)->NumTriangles);
		memset((*mesh)->colors,0,sizeof(TRGBColor)*(*mesh)->NumTriangles);

		free((*mesh)->triangles);
		free((*mesh)->normals);
		free((*mesh)->colors);
	}
	
	(*mesh)->NumTriangles=0;
	(*mesh)->NumVertices=0;
	(*mesh)->isNormalsComputed=0;
	(*mesh)->isTextureSet=0;
	(*mesh)->numConnected=0;

	free((*mesh));
	(*mesh)=0;
}

int get_mesh_vertices(TMesh* Mesh, unsigned char* verticesPC, int stride, int withNormals, int flipViewpoint)
{
	const int N = Mesh->NumVertices;
	int i;
	withNormals = (Mesh->isNormalsComputed) && withNormals;
		 
	if (withNormals)
	{
#if defined T_OPENMP
#pragma omp parallel for
#endif
		for (i=0; i<N; i++)
		{
			float* vPtr = (float*) (&verticesPC[i*stride]);
			const TVertex3D v = Mesh->vertices[i];
			TNormal3D n = Mesh->normals[i];
			const double norm = sqrt((double)n.nx*n.nx+(double)n.ny*n.ny+(double)n.nz*n.nz);
			if (norm>EPS)
			{
				n.nx/=norm;
				n.ny/=norm;
				n.nz/=norm;
			}

			if (flipViewpoint)
				flip_normal_viewpoint_32f(vPtr, 0, 0, 0, &n.nx, &n.ny, &n.nz);

			vPtr[0] = v.vx;
			vPtr[1] = v.vy;
			vPtr[2] = v.vz;
			vPtr[3] = n.nx;
			vPtr[4] = n.ny;
			vPtr[5] = n.nz;
		}
	}
	else
	{
#if defined T_OPENMP
#pragma omp parallel for
#endif
		for (i=0; i<N; i++)
		{
			float* vPtr = (float*) (&verticesPC[i*stride]);
			const TVertex3D v = Mesh->vertices[i];
			vPtr[0] = v.vx;
			vPtr[1] = v.vy;
			vPtr[2] = v.vz;
		}
	}

	return withNormals;
}
