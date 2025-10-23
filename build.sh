cd ~/BA/CardilloMPI/build
rm -rf *
#cmake .. -DCMAKE_BUILD_TYPE=Release
cmake .. -DENABLE_GPROF=OFF -DCMAKE_BUILD_TYPE=Release
make -j
