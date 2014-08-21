//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                          License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2014, OpenCV Foundation, all rights reserved.
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
// Author: Tolga Birdal

#include "tmesh.h"
#include "ppf_match_3d.hpp"
#include "icp.hpp"
#include "visualize_win.h"

using namespace std;
using namespace cv;
using namespace ppf_match_3d;

int main()
{
    int numVert ;
    
    //const char* fn = "../../../data/kinect/model/Frog_ascii2.ply";
    //const char* fn = "../../../data/SpaceTime/Scena1/scene1-model1_0_ascii.ply";
    const char* fn = "data/parasaurolophus_6700.ply";
    
    // Filenames for outputing a transformed mesh (Generally for visualization in Meshlab and etc.)
    const char* outputResultFile = "data/out/PPFICPOutput.ply";
    const char* scenePCFile = "data/out/PPFICPScene.ply";
    
    // Read the model as a mesh
    TMesh* mesh = 0;
    read_mesh_ply(&mesh, fn);
    Mat pc = Mat(mesh->NumVertices, 6, CV_32F);
    get_mesh_vertices(mesh, pc.data, pc.step, 1, 0);
    
    // Now train the model
    printf("Training...");
    int64 tick1 = cv::getTickCount();
    ppf_match_3d::PPF3DDetector detector(0.03, 0.03);
    detector.trainModel(pc);
    int64 tick2 = cv::getTickCount();
    printf("\nTraining complete in %f ms. Loading model...", (double)(tick2-tick1)/ cv::getTickFrequency());
    
    // Read the scene
    
    //numVert = 264310;
    //fn = "../../../data/kinect/scene/frog_scene_5_ascii.ply";
    //numVert = 131834;
    //fn = "../../../data/SpaceTime/Scena1/scene1-scene1_0_ascii.ply";
    numVert = 114373;
    fn = "data/rs1_normals.ply";
    Mat pcTest = loadPLYSimple(fn, 1);
    printf("\nStarting matching...");
    
    // Match the model to the scene and get the pose
    tick1 = cv::getTickCount();
    vector < Pose3D* > results;
    detector.match(pcTest, results, 1.0/10.0, 0.05);
    //t_match_pc_ppf(pcTest, 15, 5, 0.03, ppfModel, results);
    tick2 = cv::getTickCount();
    printf("\nPPF Elapsed Time %f sec\n", (double)(tick2-tick1)/ cv::getTickFrequency());
    
    // Create an instance of ICP
    ICP icp(200, 0.001f, 2.5f, 8);
    
    int64 t1 = cv::getTickCount();
    float residualOutput = 0;
    
    // Get only first N results
    int N=5;
    vector<Pose3D*>::const_iterator first = results.begin();
    vector<Pose3D*>::const_iterator last = results.begin() + N;
    vector<Pose3D*> resultsSub(first, last);
    
    // Register for all selected poses
    icp.registerModelToScene(pc, pcTest, resultsSub);
    int64 t2 = cv::getTickCount();
    
    printf("Elapsed Time: %f\n\n", (double)(t2-t1)/cv::getTickFrequency());
    
    printf("Estimated Poses:\n");
    // debug first five poses
    for (size_t i=0; i<resultsSub.size(); i++)
    {
        Pose3D* pose = resultsSub[i];
        
        printf("Pose Result %d:\n", i);
        pose->printPose();
        
        Mat pct = transformPCPose(pc, pose->Pose);
        
        // write the mesh if desired
        //TMesh* MeshTr = transform_mesh_new(mesh, pose->Pose);
        //writePLY(pct, outputResultFile);
        //t_write_mesh_ply(MeshTr, outputResultFile);
        //writePLY(pcTest, scenePCFile);
        //destroy_mesh(&MeshTr);
        
        // Visualize registration
#if defined (_MSC_VER)
        visualize_registration(pcTest, pct, "Registration");
#endif
    }
    
    return 0;
}

