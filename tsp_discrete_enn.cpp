// -*- C++ -*-
#include "tsp_discrete_enn.hpp"

#include <cstddef>
#include <cstdlib>
#include <chrono>
#include <filesystem>

using TimeMilliS_t = std::chrono::milliseconds;
using TimeMicroS_t = std::chrono::microseconds;
using TimeUnit_t = TimeMilliS_t;
using TimePoint_t = std::chrono::steady_clock::time_point;
const std::string& time_unit{ "ms" };

const std::string& Data_Dir{ "Data/ALL_tsp" };
namespace stdfs = std::filesystem;
stdfs::path Data_Path{ stdfs::current_path() / stdfs::path(Data_Dir) };
const std::string& Data_Filename_berlin{ "berlin52.tsp" };
std::string Data_Filename{ "" };
stdfs::path Data_FilePath;

int main(int argc, char** argv)
{
    std::vector<std::string> args(argc - 1);
    for (int idx{1}; idx < argc; ++idx) {
        args[idx - 1] = std::string(argv[idx]);
#if (DEBUG_PRINT > 1)
        std::cout << "#" << idx << " flag " << args[idx - 1] << '\n';
#endif
    }
    // -------------------------------------------
    // Seed and initialize rng
    // -------------------------------------------
    std::random_device::result_type seed{ std::random_device{}() };
    std::default_random_engine rng{ seed };

    // -------------------------------------------
    // Create data file
    // -------------------------------------------
    if (Data_Filename.empty()) {
        Data_Filename = Data_Filename_berlin;
    }
    Data_FilePath = Data_Path / stdfs::path(Data_Filename);

    // -------------------------------------------
    // Parse cities
    // -------------------------------------------
    Cities_t cities;
    parseCities(cities, Data_FilePath.string());
    std::shuffle(cities.begin(), cities.end(), rng);
    const int num_cities = cities.size();

    // -------------------------------------------
    // Initialize Path
    // -------------------------------------------
    Path_t path;
    initializePath(path, cities);
    if (path.size() == static_cast<std::size_t>(num_cities)) {
        std::printf("[Info]: Algorithm complete. Only %d number of cities provided.", num_cities);
        return EXIT_FAILURE;
    }

    // -------------------------------------------
    // Calculate layer details
    // -------------------------------------------
    int layers{0};
    int layers_val{1};
    while(true) {
        ++layers;
        layers_val *= 4;
        if (layers_val >= num_cities) break;
    }
    std::printf("[Info]: Total number of layers expected %d (power of 4 : %d) for number of cities %d.\n", layers, layers_val, num_cities);

    // -------------------------------------------
    // Construct Stack
    // -------------------------------------------
    Cities_t stack;
    stack.reserve(num_cities);
    createStack(cities, stack, layers);
    std::printf("[Info]: Total number of layers %d.\n", layers);

    // -------------------------------------------
    // Construct Path
    // -------------------------------------------
    constexpr std::size_t num_nodes_initial{3};
    // constructPath(path, cities, num_nodes_initial, rng);
    constructPath(path, stack, num_nodes_initial);
    // for (std::size_t idx{ 0 }; idx != num_nodes_initial; ++idx) {
    //     path[idx]->on_stack = false;
    //     nodeCost(idx, path, true);
    // }
    {

#if (DEBUG_PRINT > 0)
        std::cout << ("\n[Debug] (main): validatePath\n");
#endif
        NodeExp<bool> erased = validatePath(path, Validation_Intersection, Intersection_Recursive);
        if (erased.err()) {
            std::cerr << "[Error] (main): validatePath failed\n";
            return EXIT_FAILURE;
        }

#if (DEBUG_PRINT > 0)
        std::cout << ("\n[Debug] (main): validatePath updateCostAll\n");
#endif
        const auto [success, idx] = updateCostAll(path);
        if (not success) {
            std::cerr
                << "[Error] (main): updateCostAll failed at index "
                << idx << '\n';
            return EXIT_FAILURE;
        }
    }

    // -------------------------------------------
    // Run Discrete ENN
    // -------------------------------------------
    std::cout << ("\n[Info] (main): Run Discrete ENN Algorithm\n");
    TimePoint_t start_time = std::chrono::steady_clock::now();
    const bool success = runDiscreteENN(stack, path, rng);
    TimePoint_t end_time = std::chrono::steady_clock::now();
    if (not success) {
        std::cerr << "[Error] (main): Discrete ENN run failed.\n";
        return EXIT_FAILURE;
    }
    auto delta = std::chrono::duration_cast<TimeUnit_t>(end_time - start_time);
    const auto duration = delta.count();
    std::cout << "\n" + Line_Str + "\n";
    std::cout << "[Info] (main): Algorithm finished in " << duration << time_unit + "\n";
    std::cout << Line_Str + "\n";
    assert("[Error] (main): path size not equal to stack size" &&
           (path.size() == stack.size()));

    const int num_nodes = path.size();
    for (int idx{ 0 }; idx != num_nodes; ++idx) {
        if (not validateNode(idx, path)) {
            std::cerr << "[Error] (main): Algoirthm has not found the optimal path\n";
            return EXIT_FAILURE;
        }
    }

    // -------------------------------------------
    // Show results
    // -------------------------------------------
    std::cout << "[Info]: Print results\n";
    double dist{0.0};
    for (int idx{ 0 }; idx != num_nodes; ++idx) {
        // path[idx]->print();
        dist += distance(*path[idx], *path[properIndex(idx + 1, num_nodes)]);
        if (not validateNode(idx, path)) {
            std::cerr << "[Error]: Final path has invalid node. Exiting\n";
            return EXIT_FAILURE;
        }
    }
    std::cout << "\n" + Line_Str + "\n";
    std::cout << "[Info]: Total distance is : " << dist << '\n';
    std::cout << Line_Str + "\n";

    if (not vectContains(std::string{"--batch"}, args))
    {
        drawPath(path, stack);
    }

    return 0;
}

