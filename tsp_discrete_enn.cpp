// -*- C++ -*-
#include "tsp_discrete_enn.hpp"

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
    for (int idx{ 1 }; idx < argc; ++idx) {
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
    if (utils::vectContains(std::string{ "--input" }, args)) {
        const auto it = std::find_if(args.begin(), args.end(), utils::MatchItem<std::string>{ "--input" });
        Data_Filename = *(it + 1);
    }
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
    // Calculate layer details
    // -------------------------------------------
    int layers{ 0 };
    int layers_val{ 1 };
    while (true) {
        ++layers;
        layers_val *= 4;
        if (layers_val >= num_cities)
            break;
    }
    std::printf(
        "[Info]: Total number of layers expected %d (power of 4 : %d) for number of cities %d.\n",
        layers, layers_val, num_cities);

    // -------------------------------------------
    // Create and setup Discrete ENN Solver
    // -------------------------------------------
    DiscreteENN_TSP enn_tsp;
    enn_tsp.initialSize() = Num_Nodes_Initial;
    enn_tsp.intersection() = Validation_Intersection;
    enn_tsp.recursive() = Intersection_Recursive;
    enn_tsp.iterRandomize() = Iter_Randomize;
    enn_tsp.repeatLength() = Repeat_Check_Length;

    // -------------------------------------------
    // Construct Stack
    // -------------------------------------------
    createStack(cities, enn_tsp.stack(), layers);
    std::printf("[Info]: Total number of layers created %d.\n", layers);

    // -------------------------------------------
    // Initialize Path
    // -------------------------------------------
    enn_tsp.initializePath();
    if (enn_tsp.path().size() == static_cast<std::size_t>(num_cities)) {
        std::printf(
            "[Info]: Algorithm complete. Only %d number of cities provided.",
            num_cities);
        return EXIT_SUCCESS;
    }

    // -------------------------------------------
    // Construct Path
    // -------------------------------------------
    enn_tsp.constructPath();

    {
#if (TSP_DEBUG_PRINT > 0)
        std::cout << ("\n[Debug] (main): validatePath\n");
#endif
        NodeExp_t<bool> erased = enn_tsp.validatePath();
        if (erased.err()) {
            std::cerr << "[Error] (main): validatePath failed\n";
            return EXIT_FAILURE;
        }

#if (TSP_DEBUG_PRINT > 0)
        std::cout << ("\n[Debug] (main): validatePath updateCostAll\n");
#endif
        const auto [idx_fail, err] = enn_tsp.updateCostAll();
        if (err) {
            std::cerr << "[Error] (main): updateCostAll failed at index "
                      << idx_fail << '\n';
            return EXIT_FAILURE;
        }
    }

    // -------------------------------------------
    // Run Discrete ENN
    // -------------------------------------------
    std::cout << ("\n[Info] (main): Run Discrete ENN Algorithm\n");
    TimePoint_t start_time = std::chrono::steady_clock::now();
    const bool success = enn_tsp.run(rng);
    TimePoint_t end_time = std::chrono::steady_clock::now();
    if (not success) {
        std::cerr << "[Error] (main): Discrete ENN run failed.\n";
        return EXIT_FAILURE;
    }
    auto delta = std::chrono::duration_cast<TimeUnit_t>(end_time - start_time);
    const auto duration = delta.count();
    std::cout << "\n" + utils::Line_Str + "\n";
    std::cout << "[Info] (main): Algorithm finished in " << duration
              << time_unit + "\n";
    std::cout << utils::Line_Str + "\n";
    assert("[Error] (main): path size not equal to stack size" &&
           (enn_tsp.path().size() == enn_tsp.stack().size()));
    assert("[Error] (main): path size not equal to number of cities" &&
           (enn_tsp.path().size() == static_cast<std::size_t>(num_cities)));

    const std::size_t num_nodes = enn_tsp.path().size();
    Path_t& path = enn_tsp.path();
    for (std::size_t idx{ 0 }; idx != num_nodes; ++idx) {
        const auto [valid, err] = enn_tsp.validateNode(path[idx]);
        if (err) {
            std::cerr
                << "[Error] (main): Algoirthm has not found the optimal path\n";
            return EXIT_FAILURE;
        }
    }
    NodeExp_t<bool> erased = enn_tsp.validatePath();
    if (erased.err()) {
        std::cerr << "[Error] (main): final validatePath failed\n";
        return EXIT_FAILURE;
    }
    if (erased.has_value()) {
        std::cerr << "[Error] (main): final validatePath removed node(s)\n";
        return EXIT_FAILURE;
    }

    // -------------------------------------------
    // Show results
    // -------------------------------------------
    std::cout << "[Info]: Print results\n";
    Value_t dist{ 0.0 };
    for (std::size_t idx{ 0 }; idx != num_nodes; ++idx) {
        // path[idx]->print();
        dist += distance(*path[idx], *path[enn_tsp.properIndex(idx + 1)]);
    }
    std::cout << "\n" + utils::Line_Str + "\n";
    std::cout << "[Info]: Total distance is : " << dist << '\n';
    std::cout << utils::Line_Str + "\n";

    if (not utils::vectContains(std::string{ "--batch" }, args)) {
        const bool show_coords{ utils::vectContains(std::string{"--show-coords"}, args) };
        drawPath(path, enn_tsp.stack(), show_coords);
    }

    return 0;
}
