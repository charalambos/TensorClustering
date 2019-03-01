# TensorClustering
### by Charalambos [Charis] Poullis

The source code for TensorClustering as described in the _**IEEE TPAMI 2019 paper: Large-scale Urban Reconstruction with Tensor Clustering and Global Boundary Refinement**_.

####IMPORTANT: To use this software, YOU MUST CITE the following in any resulting publication:

@article{poullis2019large,
  title={Large-scale Urban Reconstruction with Tensor Clustering and Global Boundary Refinement},
  author={Poullis, Charalambos},
  journal={IEEE transactions on pattern analysis and machine intelligence},
  year={2019}
}

####Command-line usage:
`./TensorClustering resampled_xyz_map.pfm`

####Input
The input is a resampled XYZ map in PFM format. To generate such files from XYZ pointclouds (ASCII format) you can use  _**StructurePointcloud**_ available here: https://github.com/charalambos/StructurePointcloud

####Dependencies
To compile and build you will need the following:
1. ImageMagick > 7.x
2. FFTW3
3. OpenCV
4. MPFR  
5. GMP
6. CGAL
7. GLEW
8. OPENGL
9. Boost_SYSTEM
10. EXIV2
11. Qt5::Widgets, Qt5::Core, Qt5::OpenGL

You will also need to download our research lab's "_commons files_" from _**theICTlab**_ repository. 

