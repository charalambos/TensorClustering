cmake_minimum_required(VERSION 3.5)
project(TensorClustering)
add_definitions(-DMAGICKCORE_HDRI_ENABLE=0 -DMAGICKCORE_QUANTUM_DEPTH=16 -DNO_XML -DCMAKE_BUILD_TYPE=RELEASE)

#Set the correct output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")	#include a script to find all Find*.cmake files


## Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -std=c++11")
set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary

set(SOURCE_FILES
        ../theICTlab/Image/bilateral_filtering/fft_3D/convolution_3D.h
        ../theICTlab/Image/bilateral_filtering/fft_3D/fill_3D.h
        ../theICTlab/Image/bilateral_filtering/fft_3D/support_3D.cpp
        ../theICTlab/Image/bilateral_filtering/fft_3D/support_3D.h
        ../theICTlab/Image/bilateral_filtering/include/array.h
        ../theICTlab/Image/bilateral_filtering/include/array_n.h
        ../theICTlab/Image/bilateral_filtering/include/chrono.h
        ../theICTlab/Image/bilateral_filtering/include/fast_color_bf.h
        ../theICTlab/Image/bilateral_filtering/include/fft_3D.h
        ../theICTlab/Image/bilateral_filtering/include/geom.h
        ../theICTlab/Image/bilateral_filtering/include/linear_bf.h
        ../theICTlab/Image/bilateral_filtering/include/math_tools.h
        ../theICTlab/Image/bilateral_filtering/include/mixed_vector.h
        ../theICTlab/Image/bilateral_filtering/include/msg_stream.h
        ../theICTlab/Image/Color.cpp
        ../theICTlab/Image/Color.h
        ../theICTlab/Image/GaborFilter.cpp
        ../theICTlab/Image/GaborFilter.h
        ../theICTlab/Image/GaborJet.cpp
        ../theICTlab/Image/GaborJet.h
        ../theICTlab/Image/Image.cpp
        ../theICTlab/Image/Image.h
        ../theICTlab/Image/ImageManager.cpp
        ../theICTlab/Image/ImageManager.h
        ../theICTlab/Image/ImageProcessing.cpp
        ../theICTlab/Image/ImageProcessing.h

        ../theICTlab/Statistics/alloc_util.c
        ../theICTlab/Statistics/alloc_util.h
        ../theICTlab/Statistics/classify.c
        ../theICTlab/Statistics/classify_util.c
        ../theICTlab/Statistics/classify_util.h
        ../theICTlab/Statistics/clust.h
        ../theICTlab/Statistics/clust_defs.h
        ../theICTlab/Statistics/clust_invert.c
        ../theICTlab/Statistics/clust_io.c
        ../theICTlab/Statistics/clust_io.h
        ../theICTlab/Statistics/clust_util.c
        ../theICTlab/Statistics/clust_util.h
        ../theICTlab/Statistics/eigen.c
        ../theICTlab/Statistics/GaussianDistribution.h
        ../theICTlab/Statistics/GaussianMixtureDistribution.h
        ../theICTlab/Statistics/invert.c
        ../theICTlab/Statistics/mpreal.h
        ../theICTlab/Statistics/Ransac.h
        ../theICTlab/Statistics/subcluster.c
        ../theICTlab/Statistics/subcluster.h
        ../theICTlab/Statistics/WeibullDistribution.h
        ../theICTlab/Statistics/Ransac.h
        ../theICTlab/Statistics/RNG.h

        ../theICTlab/Math/Algebra.h
        ../theICTlab/Math/Tensor.h
        ../theICTlab/Math/Tensor.cpp

        ../theICTlab/Geometry/AdjacencyMap.cpp
        ../theICTlab/Geometry/AdjacencyMap.h
        ../theICTlab/Geometry/BoundingBox.cpp
        ../theICTlab/Geometry/BoundingBox.h
        ../theICTlab/Geometry/Contour.cpp
        ../theICTlab/Geometry/Contour.h
        ../theICTlab/Geometry/ContourExtruder.cpp
        ../theICTlab/Geometry/ContourExtruder.h
        ../theICTlab/Geometry/ComboundPatch.cpp
        ../theICTlab/Geometry/ComboundPatch.h
        ../theICTlab/Geometry/Edge.h
        ../theICTlab/Geometry/Face.h
        ../theICTlab/Geometry/GeometricObject.cpp
        ../theICTlab/Geometry/GeometricObject.h
        ../theICTlab/Geometry/GeometryExporter.cpp
        ../theICTlab/Geometry/GeometryExporter.h
        ../theICTlab/Geometry/GeometryProcessing.cpp
        ../theICTlab/Geometry/GeometryProcessing.h
        ../theICTlab/Geometry/GeometryLoader.cpp
        ../theICTlab/Geometry/GeometryLoader.h
        ../theICTlab/Geometry/GeospatialBoundingBox.cpp
        ../theICTlab/Geometry/GeospatialBoundingBox.h
        ../theICTlab/Geometry/Neighbourhood.cpp
        ../theICTlab/Geometry/Neighbourhood.h
        ../theICTlab/Geometry/ONB.cpp
        ../theICTlab/Geometry/ONB.h
        ../theICTlab/Geometry/Patch.cpp
        ../theICTlab/Geometry/Patch.h
        ../theICTlab/Geometry/PointIndex.cpp
        ../theICTlab/Geometry/PointIndex.h
        ../theICTlab/Geometry/objloader/OBJ.h
        ../theICTlab/Geometry/objloader/OBJ.cpp
        ../theICTlab/Geometry/objloader/OBJParser.h
        ../theICTlab/Geometry/objloader/OBJParser.lex.h
        ../theICTlab/Geometry/objloader/OBJParser.lex.cpp
        ../theICTlab/Geometry/objloader/OBJParser.yacc.hpp
        ../theICTlab/Geometry/objloader/OBJParser.yacc.cpp

        ../theICTlab/Geometry/Ray.h

        ../theICTlab/General/FileUtilities.h
        ../theICTlab/General/Utilities.h

        ../theICTlab/Features/Feature.h
        ../theICTlab/Features/Feature2d.h
        ../theICTlab/Features/FeatureLine.h
        ../theICTlab/Features/FeatureManager.cpp
        ../theICTlab/Features/FeatureManager.h
        ../theICTlab/Features/FeatureMarquee.h
        ../theICTlab/Features/FeaturePoint2d.h
        ../theICTlab/Features/FeaturePolyline.h

        ../theICTlab/Camera/PerspectiveCamera.cpp
        ../theICTlab/Camera/PerspectiveCamera.h

        ../theICTlab/Core/Portal.cpp
        ../theICTlab/Core/Portal.h
        ../theICTlab/Core/Renderable.h
        ../theICTlab/Core/Selectable.h
        ../theICTlab/Core/Singleton.cpp
        ../theICTlab/Core/Singleton.h
        ../theICTlab/Core/Sorter.h
        ../theICTlab/Core/ResourceManager.cpp
        ../theICTlab/Core/ResourceManager.h

        ../theICTlab/OpenGL/GLTextureManager.cpp
        ../theICTlab/OpenGL/GLTextureManager.h
        ../theICTlab/OpenGL/OpenGL.cpp
        ../theICTlab/OpenGL/OpenGL.h
        ../theICTlab/OpenGL/ShaderObject.cpp
        ../theICTlab/OpenGL/ShaderObject.h
        ../theICTlab/OpenGL/ShaderProgram.cpp
        ../theICTlab/OpenGL/ShaderProgram.h

        main.cpp)

include_directories(.)
include_directories(./src)
include_directories(./src/core)
include_directories(./src/opengl)

include_directories(../theICTlab/General)

include_directories(../theICTlab/Statistics)

include_directories(../theICTlab/Math)

include_directories(../theICTlab/Image)
include_directories(../theICTlab/Image/bilateral_filtering)
include_directories(../theICTlab/Image/bilateral_filtering/fft_3D)
include_directories(../theICTlab/Image/bilateral_filtering/include)

include_directories(../theICTlab/Geometry)
include_directories(../theICTlab/Geometry/objloader)

include_directories(../theICTlab/Core)

include_directories(../theICTlab/Camera)

include_directories(../theICTlab/Features)

include_directories(../theICTlab/OpenGL)

include_directories(/usr/include/eigen3)

# Configuring packages
## ImageMagic
find_package(ImageMagick)
find_package(ImageMagick COMPONENTS convert)
find_package(ImageMagick COMPONENTS convert mogrify display)
find_package(ImageMagick COMPONENTS Magick++)
find_package(ImageMagick COMPONENTS Magick++ convert)
#MESSAGE(STATUS "ImageMagick: " ${ImageMagick_INCLUDE_DIRS})
include_directories(${ImageMagick_INCLUDE_DIRS})

## libfftw3
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")	#include a script to find LibFFTW3
set(LIBFFTW3_PRECISION DOUBLE)				#using double precision since some functions are not implemented in float precision
find_package(LibFFTW3 REQUIRED)
include_directories(${LIBFFTW3_INCLUDE_DIR})
#include_directories(/usr/include/fftw3)

## OpenCV
find_package(OpenCV REQUIRED)
#MESSAGE(STATU "OpenCV LIBS:" ${OpenCV_LIBS})
#MESSAGE(STATU "OpenCV Includes:" ${OpenCV_INCLUDE_DIRS})
# Moving all #include "opencv.hpp" to #include "opencv2/opencv.hpp"
include_directories(/usr/include/opencv)

#MPFR
find_package(MPFR REQUIRED)
include_directories(${MPFR_INCLUDES})

#GMP
find_package(GMP REQUIRED)
include_directories(${GMP_INCLUDES})

#CGAL
find_package(CGAL REQUIRED)
include_directories(${CGAL_INCLUDE_DIR})
#MESSAGE(STATUS "CGAL: " ${CGAL_INCLUDE_DIR})

# Find the QtWidgets library
find_package(Qt5 COMPONENTS Core Widgets OpenGL REQUIRED)
#find_package(Qt5Widgets)
#find_package(Qt5Core)

set(CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS_RELEASE} -fprofile-arcs -ftest-coverage")

# set up a mapping so that the Release configuration for the Qt imported target is
# used in the COVERAGE CMake configuration.
set_target_properties(Qt5::Widgets PROPERTIES MAP_IMPORTED_CONFIG_COVERAGE "RELEASE")
set_target_properties(Qt5::Core PROPERTIES MAP_IMPORTED_CONFIG_COVERAGE "RELEASE")
set_target_properties(Qt5::OpenGL PROPERTIES MAP_IMPORTED_CONFIG_COVERAGE "RELEASE")

#Find GLM
find_package(GLM REQUIRED)
include_directories(${GLM_INCLUDE_DIRS})
#MESSAGE(STATUS "GLM: " ${GLM_INCLUDE_DIRS})

#Find GLEW
find_package(GLEW REQUIRED)
include_directories(${GLEW_INCLUDE_DIR})

#Find OPENGL
find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR})

#Find Boost::System
find_package(Boost COMPONENTS system REQUIRED)

#Find EXIV2
find_package(EXIV2 REQUIRED)
include_directories(${EXIV2_INCLUDE_DIR})

add_executable(TensorClustering ${SOURCE_FILES})
target_link_libraries(TensorClustering
        ${ImageMagick_LIBRARIES}
        ${LIBFFTW3_LIBRARY}
        ${OpenCV_LIBS}
        ${MPFR_LIBRARIES}
        ${GMP_LIBRARIES}
        ${CGAL_LIBRARIES}
        ${GLEW_LIBRARY}
        ${OPENGL_LIBRARIES}
        ${Boost_SYSTEM_LIBRARY}
        ${EXIV2_LIBRARIES}
        Qt5::Widgets
        Qt5::Core
        Qt5::OpenGL
        )
