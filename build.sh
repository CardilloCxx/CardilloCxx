cd ~/BA/CardilloMPI/build
rm -rf *
# Enable tinyexr EXR loader by default to support heightmap loading when headers/libs are available.
# Users can override by editing this script or running cmake manually with -DCARDILLO_USE_TINYEXR=OFF
cmake .. -DENABLE_GPROF=OFF -DCARDILLO_USE_TINYEXR=ON -DCMAKE_BUILD_TYPE=RELEASE
make -j
