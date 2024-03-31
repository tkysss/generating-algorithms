#pragma once

#include <iostream>
#include <vector>
#include <queue>
#include <set>
#include <bitset>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <map>
#include <fstream>
#include <sstream>
#include <chrono>

#include <omp.h>

using std::array;
using std::vector;
using std::queue;
using std::map;
using std::pair;
using std::bitset;
using std::move;
using std::cout;
using std::cerr;
using std::endl;

#define MAX_GRAPH_N 23

int DPVC_PATH_LEN = -1;
double DPVC_BF = -1;

int GENERATIONS = -1;
int NUM_THREADS = -1;
std::string OUTPUT_FILES_PREFIX = "";
