<figure style="text-align: center;">
    <img src="https://github.com/Nikita-Dudorov/transparency/blob/main/images/real1.png" width="400" />
</figure>

<figure style="text-align: center;">
    <img src="https://github.com/Nikita-Dudorov/transparency/blob/main/images/shader_n12.png" width="400" />
</figure>

# MyRenderer
## Base code for INF584 "Image Synthesis" Assignments
## The base code is modified and extended to refraction and transparency effects

"A" solution for the basic ray tracing assignment.

Author: Tamy Boubekeur (tamy.boubekeur@telecom-paris.fr)

Modified by: Nikita Dudorov

### Building

Download the glad library: https://github.com/Dav1dde/glad

This is a standard CMake project. Building it consists in running:

```
cd <path-to-MyRenderer-directory>
mkdir Build
cd Build
cmake ..
cd ..
cmake --build Build
```

The resuling MyRenderer executable is automatically copied to the root MyRenderer directory, so that resources (shaders, meshes) can be loaded easily. By default, the program is compile in Debug mode. For a high performance Release binary, just us:

```
cmake --build Build --config Release
```

### Running

To run the program
```
cd <path-to-MyRenderer-directory>
./MyRenderer [file.off]
```
Note that a collection of example meshes are provided in the Resources/Models directory.

When starting to edit the source code, rerun 

```
cmake --build Build --config Release
```
each time you want to recompile MyRenderer

The resulting binary to use is always the one located in the MyRenderer directory, you can safely ignore whatever is generated in the Build directory. 
