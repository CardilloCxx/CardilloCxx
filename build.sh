cd ~/BA/CardilloMPI/build
rm -rf *
#cmake .. -DCMAKE_BUILD_TYPE=Release
cmake .. -DENABLE_GPROF=ON -DCMAKE_BUILD_TYPE=Debug
make -j
