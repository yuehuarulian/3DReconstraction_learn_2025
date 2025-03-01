//
// Created by caoqi on 2018/10/08.
//

/*
 * Copyright (C) 2015, Simon Fuhrmann
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <iomanip>
#include <iostream>
#include <cstdlib>

#include "mvs/settings.h"
#include "mvs/dmrecon.h"
#include "core/scene.h"
#include "core/view.h"
#include "util/timer.h"
#include "util/arguments.h"
#include "util/system.h"
#include "util/tokenizer.h"
#include "util/file_system.h"

struct AppSettings
{
    std::string scene_path;
    std::string ply_dest = "recon";
    int master_id = -1;
    std::vector<int> view_ids;
    int max_pixels = 1500000;
    bool force_recon = false;
    bool write_ply = false;
    mvs::Settings mvs;
};

int main(int argc, char **argv)
{
    if (argc < 4)
    {
        std::cout << "usage: scendir scale view_id" << std::endl;
        return -1;
    }

    AppSettings conf;

    // 场景文件夹
    conf.scene_path = argv[1];
    // 获取图像尺度
    std::stringstream stream1(argv[2]);
    stream1 >> conf.mvs.scale;

    // 获取重建视角id
    std::stringstream stream2(argv[3]);
    stream2 >> conf.master_id;

    /* Load MVE scene. */
    std::cout << "Loading scene..." << std::endl;
    core::Scene::Ptr scene;
    try
    {
        scene = core::Scene::create(conf.scene_path);
        scene->get_bundle();
    }
    catch (std::exception &e)
    {
        std::cerr << "Error loading scene: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    /* Settings for Multi-view stereo */
    conf.mvs.writePlyFile = conf.write_ply;
    conf.mvs.plyPath = util::fs::join_path(conf.scene_path, conf.ply_dest);
    std::cout << "writing ply file to " << conf.mvs.plyPath << std::endl;

    util::WallTimer timer;
    if (conf.master_id >= 0)
    {

        std::cout << "Reconstructing view ID " << conf.master_id << std::endl;
        conf.mvs.refViewNr = (std::size_t)conf.master_id;
        try
        {
            // start reconstruction
            mvs::DMRecon recon(scene, conf.mvs);
            recon.start();
        }
        catch (std::exception &err)
        {
            std::cerr << err.what() << std::endl;
            return EXIT_FAILURE;
        }
    }

    std::cout << "Reconstruction took "
              << timer.get_elapsed() << "ms." << std::endl;

    /* Save scene */
    std::cout << "Saving views back to disc..." << std::endl;
    scene->save_views();

    return EXIT_SUCCESS;
}
