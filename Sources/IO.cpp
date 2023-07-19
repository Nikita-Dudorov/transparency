// ----------------------------------------------
// Polytechnique - INF584 "Image Synthesis"
//
// Base code for practical assignments.
//
// Copyright (C) 2022 Tamy Boubekeur
// All rights reserved.
// ----------------------------------------------
#include "IO.h" 

#include <iostream>
#include <fstream>
#include <exception>
#include <ios>
#include <algorithm>

#include "Console.h"

using namespace std;

std::shared_ptr<Mesh> IO::loadOFFMesh (const std::string & filename) {
    Console::print ("Start loading mesh <" + filename + ">");
    auto meshPtr = std::make_shared<Mesh>();
    ifstream in (filename.c_str ());
    if (!in) 
        throw std::ios_base::failure ("[IO][loadOFFMesh] Cannot open " + filename);
    string offString;
    unsigned int sizeV, sizeT, tmp;
    in >> offString >> sizeV >> sizeT >> tmp;
    auto & P = meshPtr->vertexPositions ();
    auto & T = meshPtr->triangleIndices ();
    P.resize (sizeV);
    T.resize (sizeT);
    size_t tracker = std::max<size_t> (1, (sizeV + sizeT)/20);
    Console::print (" > [", false);
    for (unsigned int i = 0; i < sizeV; i++) {
        if (i % tracker == 0)
            Console::print ("-",false);
        in >> P[i][0] >> P[i][1] >> P[i][2];
    }
    int s;
    for (unsigned int i = 0; i < sizeT; i++) {
        if ((sizeV + i) % tracker == 0)
            Console::print ("-",false);
        in >> s;
        for (unsigned int j = 0; j < 3; j++) 
            in >> T[i][j];
    }
    Console::print ("]\n", false);
    in.close ();
    meshPtr->recomputePerVertexNormals ();
    Console::print ("Mesh <" + filename + "> loaded");
    return meshPtr;
}