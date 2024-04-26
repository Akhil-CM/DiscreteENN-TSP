// -*- C++ -*-
#include "tsp_discrete_enn.hpp"
#include "utils.hpp"
#include <cstddef>

namespace stdfs = std::filesystem;

const std::string& Data_Optimal_Filename{ "./tsp_optimal_distances.csv" };
const stdfs::path Data_Optimal_Path{ utils::getCleanPath(
    stdfs::current_path() / stdfs::path(Data_Optimal_Filename)) };

const std::string& Data_Dir{ "Data/ALL_tsp" };

const std::string& Data_Filename_berlin{ "berlin52.tsp" };

int runPipelineSingle(TSPInfo& info, const stdfs::path& data_path,
                      std::default_random_engine& rng, bool draw,
                      bool draw_failed, bool show_coords, bool randomize);
int runPipelineDir(TSPInfoVect_t& infos, const TSPInfoVect_t& opt_infos,
                   const stdfs::path& data_path,
                   std::default_random_engine& rng, bool draw, bool draw_failed,
                   bool show_coords, bool randomize);

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
    const bool draw_path{ utils::vectContains(std::string{ "--draw" }, args) };
    const bool draw_coords{ utils::vectContains(std::string{ "--show-coords" },
                                                args) };
    const bool draw_failed{ utils::vectContains(std::string{ "--draw-failed" },
                                                args) };
    const bool randomize{ utils::vectContains(std::string{ "--random" },
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
        runs_failed = runPipelineSingle(info, Data_FilePath, rng, draw_path,
                                        draw_failed, draw_coords, randomize);
        if (runs_failed == 0) {
            infos.push_back(info);
        } else {
            utils::printErr("Single input pipeline failed for Data_FilePath: " +
                                Data_FilePath.string(),
                            "main");
        }
    } else {
        runs_failed = runPipelineDir(infos, optimal_infos, Data_Path, rng,
                                     draw_path, draw_failed, draw_coords,
                                     randomize);
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
        info.m_error = utils::getRound2(error * 100);
        info.m_distance = utils::getRound2(info.m_distance);
        optimal_info.m_distance = utils::getRound2(optimal_info.m_distance);
    }
    // const std::string& table_header{ "name\tpoints\terror\ttime(" + time_unit + ")\ttime per city(" + time_unit + ")\ttime per city min(" + time_unit + ")\ttime per city max(" + time_unit + ")\tdistance\toptimal_distance" };
    const std::string& table_header{ "name\tpoints\terror\ttime(" + time_unit +
                                     ")\ttime per city(" + time_unit +
                                     ")\tdistance\toptimal_distance" };
    utils::printInfo(table_header);
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
        utils::printInfo(info.m_name + "\t" + std::to_string(info.m_points) +
                         "\t" + std::to_string(info.m_error) + "\t" +
                         std::to_string(info.m_time) + "\t" +
                         std::to_string(info.m_timePerCity) + "\t" +
                         // std::to_string(info.m_timePerCityMin) + "\t" +
                         // std::to_string(info.m_timePerCityMax) + "\t" +
                         std::to_string(info.m_distance) + "\t" +
                         std::to_string(opt_info.m_distance));
    }
    std::ofstream table_file{ "DiscreteENN_TSP_table.txt" };
    table_file << table_header << std::endl;
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
        table_file << (info.m_name + "\t" + std::to_string(info.m_points) +
                       "\t" + std::to_string(info.m_error) + "\t" +
                       std::to_string(info.m_time) + "\t" +
                       std::to_string(info.m_timePerCity) + "\t" +
                       // std::to_string(info.m_timePerCityMin) + "\t" +
                       // std::to_string(info.m_timePerCityMax) + "\t" +
                       std::to_string(info.m_distance) + "\t" +
                       std::to_string(opt_info.m_distance))
                   << std::endl;
    }
    table_file.close();
    std::ofstream csv_file{ "DiscreteENN_TSP_table.csv" };
    csv_file << table_header << std::endl;
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
        csv_file << (info.m_name + "," + std::to_string(info.m_points) + "," +
                     std::to_string(info.m_error) + "," +
                     std::to_string(info.m_time) + "," +
                     std::to_string(info.m_timePerCity) + "," +
                     // std::to_string(info.m_timePerCityMin) + "," +
                     // std::to_string(info.m_timePerCityMax) + "," +
                     std::to_string(info.m_distance) + "," +
                     std::to_string(opt_info.m_distance))
                 << std::endl;
    }
    csv_file.close();

    return 0;
}

int runPipelineDir(TSPInfoVect_t& infos, const TSPInfoVect_t& opt_infos,
                   const stdfs::path& data_path,
                   std::default_random_engine& rng, bool draw, bool draw_failed,
                   bool show_coords, bool randomize)
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
        if (runPipelineSingle(info, filepath, rng, draw, draw_failed,
                              show_coords, randomize) != 0) {
            utils::printErr("pipeline failed for the path " + filepath.string(),
                            "runPipelineDir");
            ++runs_failed;
            continue;
        }
        infos.push_back(info);
    }

    return runs_failed;
}

int runPipelineSingle(TSPInfo& info, const stdfs::path& data_path,
                      std::default_random_engine& rng, bool draw,
                      bool draw_failed, bool show_coords, bool randomize)
{
    utils::printInfo("Running algorithm for " + data_path.string(),
                     "runPipelineSingle");

    const std::string& filename{ data_path.stem().string() };
    // -------------------------------------------
    // Create and setup Discrete ENN Solver
    // -------------------------------------------
    DiscreteENN_TSP enn_tsp;
    enn_tsp.name() = filename;
    Cities_t& cities{ enn_tsp.cities() };
    Indices_t& path = enn_tsp.path();

    {
        // -------------------------------------------
        // Parse cities
        // -------------------------------------------
        Cities_t cities_tmp;
        parseCities(cities_tmp, data_path.string());
        if (randomize) {
            std::shuffle(cities_tmp.begin(), cities_tmp.end(), rng);
        }
        const int num_cities = cities_tmp.size();

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
        // Construct Stack
        // -------------------------------------------
        createStack(cities_tmp, cities, layers);
        std::printf("[Info]: Total number of layers created %d.\n", layers);
    }
    const std::size_t num_cities = cities.size();

    // -------------------------------------------
    // Initialize Path
    // -------------------------------------------
    enn_tsp.initializePath();
    if (path.size() == static_cast<std::size_t>(num_cities)) {
        std::printf(
            "[Info]: Algorithm complete. Only %zu number of cities provided.",
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
        IndexExp_t<bool> erased = enn_tsp.validatePath();
        if (erased.err()) {
            std::cerr << "[Error] (main): validatePath failed\n";
            return 1;
        }

#if (TSP_DEBUG_PRINT > 0)
        std::cout << ("\n[Debug] (main): validatePath updateCostAll\n");
#endif
        if (erased.has_value() and path.size() > 2) {
            const auto [idx_fail, err] = enn_tsp.updateCostAll();
            if (err) {
                std::cerr << "[Error] (main): updateCostAll failed at index "
                          << idx_fail << '\n';
                return 1;
            }
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
    assert("[Error] (main): path size not equal to number of cities" &&
           (path.size() == num_cities));

    const std::size_t num_nodes = path.size();
    for (std::size_t idx{ 0 }; idx != num_nodes; ++idx) {
        const auto [valid, err] = enn_tsp.validateNode(idx);
        if (err) {
            std::cerr
                << "[Error] (main): Algoirthm has not found the optimal path\n";
            return 1;
        }
    }
    IndexExp_t<bool> erased = enn_tsp.validatePath();
    if (erased.err()) {
        std::cerr << "[Error] (main): final validatePath failed\n";
        if (draw_failed) {
            drawPath(path, cities, show_coords, filename);
        }
        return 1;
    }
    if (erased.has_value()) {
        std::cerr << "[Error] (main): final validatePath removed node(s)\n";
        if (draw_failed) {
            drawPath(path, cities, show_coords, filename);
        }
        return 1;
    }
    if (enn_tsp.checkIntersectPath()) {
        std::cerr << "[Error] (main): final checkIntersectPath failed\n";
        if (draw_failed) {
            drawPath(path, cities, show_coords, filename);
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
        const auto idx_next{ enn_tsp.properIndex(idx + 1) };
        dist += getDistance(cities[path[idx]], cities[path[idx_next]]);
    }
    info.m_name = filename;
    info.m_distance = dist;
    info.m_points = num_cities;
    info.m_time = utils::getRound2(duration);
    info.m_timePerCity =
        utils::getRound2(static_cast<Value_t>(duration) / num_cities);
    // info.m_timePerIter = enn_tsp.timePerCity();
    // std::tie(info.m_timePerCityMin, info.m_timePerCityMax) = enn_tsp.timePerCityMinMax();
    std::cout << "\n" + utils::Line_Str + "\n";
    std::cout << "[Info]: Total distance is : " << dist << '\n';
    std::cout << utils::Line_Str + "\n";

    if (draw) {
        drawPath(path, cities, show_coords, filename);
    }
    utils::printInfo("Finished algorithm for " + data_path.string(),
                     "runPipelineSingle");
    return 0;
}
