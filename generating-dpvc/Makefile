boost_build_dir = ./boost_build

gen_no_mods = -DNOREDUCTION_RULES -DNODOMINANCE -DNOADJUSTMENT
gen_no_rr = -DNOREDUCTION_RULES
gen_rr_only = -DNODOMINANCE -DNOADJUSTMENT
gen_rr_dom_only = -DNOADJUSTMENT
gen_no_proof = -DNOPROOF

gcc_nautyc = src/nauty27rc2/nauty.c src/nauty27rc2/nautil.c src/nauty27rc2/naugraph.c src/nauty27rc2/schreier.c src/nauty27rc2/naurng.c
gcc_boosti = -I ${boost_build_dir}/include
gcc_boostl = -L ${boost_build_dir}/lib -static -lboost_iostreams -lz
gcc_ignore_errors = -Wno-write-strings
gcc_flags = -std=c++11 -O4 -march=native -fopenmp -g
gcc_compile = g++ ${gcc_flags} ${gcc_boosti} ${gcc_nautyc} ${gcc_ignore_errors}

build/generating_algorithm: src/generating_algorithm.cpp | build
	${gcc_compile} src/generating_algorithm.cpp ${gcc_boostl} -o build/generating_algorithm

build/generating_algorithm_no_mods: src/generating_algorithm.cpp | build
	${gcc_compile} ${gen_no_mods} src/generating_algorithm.cpp ${gcc_boostl} -o build/generating_algorithm_no_mods

build/generating_algorithm_rr_only: src/generating_algorithm.cpp | build
	${gcc_compile} ${gen_rr_only} src/generating_algorithm.cpp ${gcc_boostl} -o build/generating_algorithm_rr_only

build/generating_algorithm_rr_dom_only: src/generating_algorithm.cpp | build
	${gcc_compile} ${gen_rr_dom_only} src/generating_algorithm.cpp ${gcc_boostl} -o build/generating_algorithm_rr_dom_only

build/generating_algorithm_no_rr: src/generating_algorithm.cpp | build
	${gcc_compile} ${gen_no_rr} src/generating_algorithm.cpp ${gcc_boostl} -o build/generating_algorithm_no_rr

build/generating_algorithm_no_proof: src/generating_algorithm.cpp | build
	${gcc_compile} ${gen_no_proof} src/generating_algorithm.cpp ${gcc_boostl} -o build/generating_algorithm_no_proof

build:
	mkdir -p build

build_all: | build/generating_algorithm build/generating_algorithm_no_mods build/generating_algorithm_rr_only build/generating_algorithm_rr_dom_only build/generating_algorithm_no_rr build/generating_algorithm_no_proof
.PHONY: build_all
