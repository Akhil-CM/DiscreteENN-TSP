// -*- C++ -*-
#include "tsp_discrete_enn.hpp"
#include "utils.hpp"

#include <chrono>

namespace stdfs = std::filesystem;

using TimeMilliS_t = std::chrono::milliseconds;
using TimeMicroS_t = std::chrono::microseconds;
using TimeUnit_t = TimeMilliS_t;
using TimePoint_t = std::chrono::steady_clock::time_point;
const std::string& time_unit{ "ms" };

const std::string& Data_Optimal_Filename{ "./tsp_optimal_distances.csv" };
const stdfs::path Data_Optimal_Path{ utils::getCleanPath(
    stdfs::current_path() / stdfs::path(Data_Optimal_Filename)) };

const std::string& Data_Dir{ "Data/ALL_tsp" };

const std::string& Data_Filename_berlin{ "berlin52.tsp" };

int runPipelineSingle(const stdfs::path& data_path,
                      std::default_random_engine& rng, bool draw,
                      bool show_coords, TSPInfo& info);
int runPipelineDir(const stdfs::path& data_path,
                   std::default_random_engine& rng, bool draw, bool show_coords,
                   TSPInfoVect_t& infos, const TSPInfoVect_t& opt_infos);

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
    TSPInfoVect_t optimal_infos;
    if (not readDistances(Data_Optimal_Path.string(), optimal_infos)) {
        utils::printErr("Couldn't read optimal distances from " +
                            Data_Optimal_Path.string(),
                        "main");
        return EXIT_FAILURE;
    }

    const bool single_input{ utils::vectContains(std::string{ "--single" },
                                                 args) };
    const bool draw_path{ not utils::vectContains(std::string{ "--batch" },
                                                  args) };
    const bool draw_coords{ utils::vectContains(std::string{ "--show-coords" },
                                                args) };

    std::string Data_Filename{ Data_Filename_berlin };
    stdfs::path Data_Path{ utils::getCleanPath(stdfs::current_path() /
                                               stdfs::path(Data_Dir)) };
    if (utils::vectContains(std::string{ "--input" }, args)) {
        const auto it =
            std::find_if(args.begin(), args.end(),
                         utils::MatchItem<std::string>{ "--input" });
        const std::string& input{ *(it + 1) };
        if (single_input) {
            Data_Filename = input;
            if (Data_Filename.find(".tsp") == std::string::npos) {
                Data_Filename += ".tsp";
            }
        } else {
            Data_Path = utils::getCleanPath(input);
        }
    }

    int runs_failed;
    TSPInfoVect_t infos;
    infos.reserve(optimal_infos.size());
    if (single_input) {
        stdfs::path Data_FilePath = Data_Path / stdfs::path(Data_Filename);
        TSPInfo info{};
        runs_failed =
            runPipelineSingle(Data_FilePath, rng, draw_path, draw_coords, info);
        if (runs_failed == 0) {
            infos.push_back(info);
        } else {
            utils::printErr("Single input pipeline failed for Data_FilePath: " +
                                Data_FilePath.string(),
                            "main");
        }
    } else {
        runs_failed = runPipelineDir(Data_Path, rng, draw_path, draw_coords,
                                     infos, optimal_infos);
    }
    if (runs_failed != 0) {
        utils::printErr("Runs failed: " + std::to_string(runs_failed) +
                            " Data_Path: " + Data_Path.string(),
                        "main");
        return EXIT_FAILURE;
    }
    for (auto it{ infos.begin() }; it != infos.end(); ++it) {
        TSPInfo& info{ *it };
        int opt_info_pos{ utils::vectFind(info, optimal_infos) };
        if (opt_info_pos == -1) {
            utils::printErr("data name : " + info.m_name +
                                " not found in optimal_distance_map",
                            "main");
            continue;
        }
        // if (name == "pr2392") {
        //     utils::printInfo("skipping file pr2392 because of long run time.", "runPipelineDir");
        //     continue;;
        // }
        TSPInfo& optimal_info{
            optimal_infos[static_cast<std::size_t>(opt_info_pos)]
        };
        const auto error = std::abs(info.m_distance - optimal_info.m_distance) /
                           optimal_info.m_distance;
        info.m_error = error*100;
    }
    utils::printInfo("name\tdistance\tpoints\terror\ttime(ms)");
    for (auto it{ optimal_infos.begin() }; it != optimal_infos.end(); ++it) {
        TSPInfo& opt_info{ *it };
        int info_pos{ utils::vectFind(opt_info, infos) };
        if (info_pos == -1) {
            // utils::printInfo("skipping file " + opt_info.m_name +
            //                      " for which algorithm was not run.",
            //                  "main");
            continue;
        }
        TSPInfo& info{ infos[static_cast<std::size_t>(info_pos)] };
        utils::printInfo(info.m_name + "\t" + std::to_string(info.m_distance) +
                         "\t" + std::to_string(info.m_points) + "\t" +
                         std::to_string(info.m_error) + "\t" +
                         std::to_string(info.m_time));
    }
    std::ofstream table_file{ "DiscreteENN_TSP_table.txt" };
    table_file << "name\tdistance\tpoints\terror\ttime(ms)\n";
    for (auto it{ optimal_infos.begin() }; it != optimal_infos.end(); ++it) {
        TSPInfo& opt_info{ *it };
        int info_pos{ utils::vectFind(opt_info, infos) };
        if (info_pos == -1) {
            // utils::printInfo("skipping file " + opt_info.m_name +
            //                      " for which algorithm was not run.",
            //                  "main");
            continue;
        }
        TSPInfo& info{ infos[static_cast<std::size_t>(info_pos)] };
        table_file << (info.m_name + "\t" + std::to_string(info.m_distance) +
                       "\t" + std::to_string(info.m_points) + "\t" +
                       std::to_string(info.m_error) + "\t" +
                       std::to_string(info.m_time))
                   << std::endl;
    }
    table_file.close();
    std::ofstream csv_file{ "DiscreteENN_TSP_table.csv" };
    csv_file << "name,distance,points,error,time(ms)\n";
    for (auto it{ optimal_infos.begin() }; it != optimal_infos.end(); ++it) {
        TSPInfo& opt_info{ *it };
        int info_pos{ utils::vectFind(opt_info, infos) };
        if (info_pos == -1) {
            // utils::printInfo("skipping file " + opt_info.m_name +
            //                      " for which algorithm was not run.",
            //                  "main");
            continue;
        }
        TSPInfo& info{ infos[static_cast<std::size_t>(info_pos)] };
        csv_file << (info.m_name + "," + std::to_string(info.m_distance) +
                     "," + std::to_string(info.m_points) + "," +
                     std::to_string(info.m_error) + "," +
                     std::to_string(info.m_time))
                 << std::endl;
    }
    csv_file.close();

    return 0;
}

int runPipelineDir(const stdfs::path& data_path,
                   std::default_random_engine& rng, bool draw, bool show_coords,
                   TSPInfoVect_t& infos, const TSPInfoVect_t& opt_infos)
{
    if (not stdfs::is_directory(data_path)) {
        utils::printErr("provided path " + data_path.string() +
                            " doesn't exit.",
                        "runPipelineDir");
        return 1;
    }
    infos.clear();
    int runs_failed{ 0 };
    for (const auto& entry : stdfs::directory_iterator(data_path)) {
        const auto filepath = entry.path();
        if (entry.is_directory()) {
            utils::printInfo("skipping directory " + filepath.string(),
                             "runPipelineDir");
            continue;
        }
        TSPInfo info;
        const std::string& filename{ filepath.stem().string() };
        info.m_name = filename;
        int opt_info_pos{ utils::vectFind(info, opt_infos) };
        if (opt_info_pos == -1) {
            utils::printInfo("skipping file " + filepath.string() +
                                 " without optimal distance for TSP.",
                             "runPipelineDir");
            continue;
        }
        // if (filename == "pr2392") {
        //     utils::printInfo("skipping file pr2392 because of long run time.", "runPipelineDir");
        //     continue;;
        // }
        if (runPipelineSingle(filepath, rng, draw, show_coords, info) != 0) {
            utils::printErr("pipeline failed for the path " + filepath.string(),
                            "runPipelineDir");
            ++runs_failed;
            continue;
        }
        infos.push_back(info);
    }

    return runs_failed;
}

int runPipelineSingle(const stdfs::path& data_path,
                      std::default_random_engine& rng, bool draw,
                      bool show_coords, TSPInfo& info)
{
    utils::printInfo("Running algorithm for " + data_path.string(),
                     "runPipelineSingle");
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

    Path_t& path = enn_tsp.path();
    const std::size_t num_nodes = path.size();
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
        if (draw) {
            drawPath(path, enn_tsp.stack(), show_coords);
        }
        return 1;
    }
    if (erased.has_value()) {
        std::cerr << "[Error] (main): final validatePath removed node(s)\n";
        if (draw) {
            drawPath(path, enn_tsp.stack(), show_coords);
        }
        return 1;
    }
    if (enn_tsp.checkIntersectPath()) {
        std::cerr << "[Error] (main): final checkIntersectPath failed\n";
        if (draw) {
            drawPath(path, enn_tsp.stack(), show_coords);
        }
        return 1;
    }

    // -------------------------------------------
    // Show results
    // -------------------------------------------
    std::cout << "[Info]: Print results\n";
    Value_t dist{ 0.0 };
    for (std::size_t idx{ 0 }; idx != num_nodes; ++idx) {
        // path[idx]->print();
        const auto idx_next{ static_cast<std::size_t>(
            enn_tsp.properIndex(idx + 1)) };
        dist += getDistance(*path[idx], *path[idx_next]);
    }
    const std::string& filename{ data_path.stem().string() };
    info.m_name = filename;
    info.m_distance = dist;
    info.m_points = num_cities;
    info.m_time = duration;
    std::cout << "\n" + utils::Line_Str + "\n";
    std::cout << "[Info]: Total distance is : " << dist << '\n';
    std::cout << utils::Line_Str + "\n";

    if (draw) {
        drawPath(path, enn_tsp.stack(), show_coords);
    }
    utils::printInfo("Finished algorithm for " + data_path.string(),
                     "runPipelineSingle");
    return 0;
}
