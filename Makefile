

fullbuild:
	@rm -rf ./build/* || true
	@mkdir -p build
	@cmake -B build -S .
	@cmake --build build
	@mpirun -np 8 ./build/bin/main

rebuild:
	@cmake --build build
	@mpirun -np 8 ./build/bin/main


