// -*- C++ -*-
#include "tsp_discrete_enn.hpp"
#include "utils.hpp"

#include <filesystem>
#include <chrono>

namespace stdfs = std::filesystem;

using TimeMilliS_t = std::chrono::milliseconds;
using TimeMicroS_t = std::chrono::microseconds;
using TimeUnit_t = TimeMilliS_t;
using TimePoint_t = std::chrono::steady_clock::time_point;
const std::string& time_unit{ "ms" };

template <typename F,
          typename std::enable_if<std::is_convertible_v<F, stdfs::path>>::type* =
              nullptr>
stdfs::path getCleanPath(const F& src)
{
    const stdfs::path tmp_src{ src };
    const stdfs::path lexical_src{ tmp_src.lexically_normal() };
    const stdfs::path abs_src{ stdfs::absolute(lexical_src) };
    return stdfs::weakly_canonical(abs_src);
}

const std::string& Data_Optimal_Filename{ "./tsp_optimal_distances.csv" };
const stdfs::path Data_Optimal_Path{ getCleanPath(
    stdfs::current_path() / stdfs::path(Data_Optimal_Filename)) };

const std::string& Data_Dir{ "Data/ALL_tsp" };
stdfs::path Data_Path{ getCleanPath(stdfs::current_path() /
                                          stdfs::path(Data_Dir)) };

const std::string& Data_Filename_berlin{ "berlin52.tsp" };
std::string Data_Filename{ "" };
stdfs::path Data_FilePath;

int runPipelineSingle(const stdfs::path&data_path, std::default_random_engine& rng, bool draw, bool show_coords, DistanceMap_t::mapped_type& distance);
int runPipelineDir(const stdfs::path& data_path, std::default_random_engine& rng, bool draw, bool show_coords, DistanceMap_t& distance_map, const DistanceMap_t& opt_distance_map);

int main(int argc, char** argv)
{
    const auto args_count{ static_cast<std::size_t>(argc) };
    std::vector<std::string> args(args_count - 1);
    for (std::size_t idx{ 1 }; idx < args_count; ++idx) {
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
    // Read optimal distances
    // -------------------------------------------
    DistanceMap_t optimal_distance_map;
    if (not readDistances(Data_Optimal_Path.string(), optimal_distance_map)) {
        utils::printErr("Couldn't read optimal distances from " + Data_Optimal_Path.string(), "main");
    }

    const bool single_input{ utils::vectContains(std::string{ "--single" }, args) };
    const bool draw_path{ not utils::vectContains(std::string{ "--batch" }, args) };
    const bool draw_coords{ utils::vectContains(std::string{ "--show-coords" }, args) };
    if (utils::vectContains(std::string{ "--input" }, args)) {
        const auto it =
            std::find_if(args.begin(), args.end(),
                         utils::MatchItem<std::string>{ "--input" });
        Data_Filename = *(it + 1);
    }
    int runs_failed;
    DistanceMap_t distance_map;
    if(single_input) {
        if (Data_Filename.empty()) {
            Data_Filename = Data_Filename_berlin;
        }
        Data_FilePath = Data_Path / stdfs::path(Data_Filename);
        DistanceMap_t::mapped_type info{};
        runs_failed = runPipelineSingle(Data_FilePath, rng, draw_path, draw_coords, info);
        if (runs_failed == 0) {
            const std::string& filename{Data_FilePath.stem().string()};
            distance_map[filename] = info;
        } else {
            utils::printErr("Single input pipeline failed for path : " + Data_FilePath.string(), "main");
        }
    } else {
        if (not Data_Filename.empty()) {
            Data_Path = getCleanPath(Data_Filename);
        }
        runs_failed = runPipelineDir(Data_Path, rng, draw_path, draw_coords, distance_map, optimal_distance_map);
    }
    if (runs_failed != 0) {
        utils::printErr("Runs failed " + std::to_string(runs_failed) + " Data_Path " + Data_Path.string() + " Data_FilePath " + Data_FilePath.string(), "main");
        return EXIT_FAILURE;
    }
    for (auto it{ distance_map.begin() }; it != distance_map.end(); ++it) {
        const auto [name, info] = *it;
        if (optimal_distance_map.count(name) == 0) {
            utils::printErr("data name : " + name + " not found in optimal_distance_map", "main");
            continue;
        }
        // if (name == "pr2392") {
        //     utils::printInfo("skipping file pr2392 because of long run time.", "runPipelineDir");
        //     continue;;
        // }
        auto& optimal_info{ optimal_distance_map[name] };
        const auto error = std::abs(info[0] - optimal_info[0])/optimal_info[0];
        distance_map[name][2] = error;
    }
    utils::printInfo("name\tdistance\tpoints\terror\ttime(ms)");
    for (auto it{  optimal_distance_map.begin() }; it !=  optimal_distance_map.end(); ++it) {
        const auto [name, discard] = *it;
        const auto info{ distance_map[name] };
        utils::printInfo(name + "\t" + std::to_string(info[0]) + "\t" + std::to_string(info[1]) + "\t" + std::to_string(info[2]) + "\t" + std::to_string(info[3]));
    }
    std::ofstream table_file{"DiscreteENN_TSP_table.txt"};
    table_file << "name\tdistance\tpoints\terror\ttime(ms)\n";
    for (auto it{  optimal_distance_map.begin() }; it !=  optimal_distance_map.end(); ++it) {
        const auto [name, discard] = *it;
        const auto info{ distance_map[name] };
        table_file << (name + "\t" + std::to_string(info[0]) + "\t" + std::to_string(info[1]) + "\t" + std::to_string(info[2]) + "\t" + std::to_string(info[3])) << std::endl;
    }
    table_file.close();
    std::ofstream csv_file{"DiscreteENN_TSP_table.csv"};
    csv_file << "name,distance,points,error,time(ms)\n";
    for (auto it{  optimal_distance_map.begin() }; it !=  optimal_distance_map.end(); ++it) {
        const auto [name, discard] = *it;
        const auto info{ distance_map[name] };
        csv_file << (name + "," + std::to_string(info[0]) + "," + std::to_string(info[1]) + "," + std::to_string(info[2]) + "," + std::to_string(info[3])) << std::endl;
    }
    csv_file.close();

    return 0;
}

int runPipelineDir(const stdfs::path& data_path, std::default_random_engine& rng, bool draw, bool show_coords, DistanceMap_t& distance_map, const DistanceMap_t& opt_distance_map)
{

    if (not stdfs::is_directory(data_path)) {
        utils::printErr("provided path " + data_path.string() + " doesn't exit.", "runPipelineDir");
        return 1;
    }
    distance_map.clear();
    int runs_failed{0};
    for (const auto& entry : stdfs::directory_iterator(data_path)) {
        const auto filepath = entry.path();
        if (entry.is_directory()) {
            utils::printInfo("skipping directory " + filepath.string(), "runPipelineDir");
            continue;;
        }
        const std::string& filename{filepath.stem().string()};
        if (opt_distance_map.count(filename) == 0) {
            utils::printInfo("skipping file " + filepath.string() + " without optimal distance for TSP.", "runPipelineDir");
            continue;;
        }
        // if (filename == "pr2392") {
        //     utils::printInfo("skipping file pr2392 because of long run time.", "runPipelineDir");
        //     continue;;
        // }
        DistanceMap_t::mapped_type info{};
        if (runPipelineSingle(filepath, rng, draw, show_coords, info) != 0) {
            utils::printErr("pipeline failed for the path " + filepath.string(), "runPipelineDir");
            ++runs_failed;
            continue;
        }
        distance_map[filename] = info;
    }

    return runs_failed;
}

int runPipelineSingle(const stdfs::path&data_path, std::default_random_engine& rng, bool draw, bool show_coords, DistanceMap_t::mapped_type& info)
{

    utils::printInfo("Running algorithm for " + data_path.string(), "runPipelineSingle");
    // -------------------------------------------
    // Parse cities
    // -------------------------------------------
    Cities_t cities;
    parseCities(cities, data_path.string());
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
    enn_tsp.validIntersectCK() = Validation_Intersection;
    enn_tsp.rmIntersectRecurse() = Intersection_Recursive;
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
        return 0;
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
            return 1;
        }

#if (TSP_DEBUG_PRINT > 0)
        std::cout << ("\n[Debug] (main): validatePath updateCostAll\n");
#endif
        const auto [idx_fail, err] = enn_tsp.updateCostAll();
        if (err) {
            std::cerr << "[Error] (main): updateCostAll failed at index "
                      << idx_fail << '\n';
            return 1;
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
        return 1;
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
            return 1;
        }
    }
    NodeExp_t<bool> erased = enn_tsp.validatePath();
    if (erased.err()) {
        std::cerr << "[Error] (main): final validatePath failed\n";
        return 1;
    }
    if (erased.has_value()) {
        std::cerr << "[Error] (main): final validatePath removed node(s)\n";
        return 1;
    }

    // -------------------------------------------
    // Show results
    // -------------------------------------------
    std::cout << "[Info]: Print results\n";
    Value_t dist{ 0.0 };
    for (std::size_t idx{ 0 }; idx != num_nodes; ++idx) {
        // path[idx]->print();
        const auto idx_next{ static_cast<std::size_t>(enn_tsp.properIndex(idx + 1)) };
        dist += getDistance(*path[idx], *path[idx_next]);
    }
    info[0] = dist;
    info[1] = num_cities;
    info[3] = duration;
    std::cout << "\n" + utils::Line_Str + "\n";
    std::cout << "[Info]: Total distance is : " << dist << '\n';
    std::cout << utils::Line_Str + "\n";

    if (draw) {
        drawPath(path, enn_tsp.stack(), show_coords);
    }
    utils::printInfo("Finished algorithm for " + data_path.string(), "runPipelineSingle");
    return 0;
}
