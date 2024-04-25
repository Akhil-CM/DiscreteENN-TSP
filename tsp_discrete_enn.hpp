// -*- C++ Header -*-

#pragma once

#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <random>
#include <numeric>
#include <map>
#include <string>
#include <unordered_set>
#include <cmath>
#include <tuple>
#include <utility>
#include <array>
#include <chrono>

#include "utils.hpp"

#define TSP_DRAW

#ifdef TSP_DRAW
#include <SFML/Graphics.hpp>
#endif

#define TSP_DEBUG_PRINT 0
#define TSP_DEBUG_CHECK 0

using TimeMilliS_t = std::chrono::milliseconds;
using TimeMicroS_t = std::chrono::microseconds;
// using TimeUnit_t = TimeMilliS_t;
using TimeUnit_t = TimeMicroS_t;
using TimePoint_t = std::chrono::steady_clock::time_point;
const std::string& time_unit{ "us" };

using namespace std::string_literals;
using utils::ErrorBool, utils::ErrorMsg;

typedef float Value_t;
constexpr Value_t VALUE_ZERO{ Value_t{ 0 } };
constexpr Value_t VALUE_ONE_NEG{ Value_t{ -1 } };

constexpr std::size_t Num_Nodes_Initial{ 3 };
constexpr bool Validation_Intersection{ true };
constexpr bool Intersection_Recursive{ true };
constexpr int Iter_Randomize{ 50 };
constexpr int Repeat_Check_Length{ 10 };

struct City
{
    bool on_stack;
    int id, layer;
    Value_t x, y, cost;

    City(int id, Value_t x, Value_t y)
        : on_stack{ false }
        , id{ id }
        , layer{ -1 }
        , x{ x }
        , y{ y }
        , cost{ 0.0 }
    {
    }

    void print() const
    {
        std::printf("%s\n", utils::Line_Str.c_str());
        std::printf("[Info]: City details\n");
        std::printf("%s\n", utils::Line_Str.c_str());
        std::printf("ID : %d\n", id);
        std::printf("Layer : %d\n", layer);
        std::printf("(x, y) : (%f, %f)\n", x, y);
        std::printf("Cost : %f\n", cost);
        std::printf("on_stack : %s\n", on_stack ? "True" : "False");
    }
};

typedef std::vector<City> Cities_t;
typedef Cities_t::iterator Node_t;
typedef std::vector<Node_t> Path_t;
typedef std::map<int, Cities_t> CityLayers_t;
typedef std::vector<std::size_t> Indices_t;
typedef std::optional<std::size_t> IndexOpt_t;
template <typename T> using IndexExp_t = utils::Expected<std::size_t, T>;

void drawPath(const Indices_t& path, const Cities_t& cities, bool show_coords);

template <> inline bool utils::MatchItem<City>::operator()(const City& city)
{
    return (city.id == m_item.id);
}

template <> inline bool utils::MatchItem<Node_t>::operator()(const Node_t& node)
{
    return (node->id == m_item->id);
}

void parseCities(Cities_t& cities, const std::string& filename)
{
    cities.clear();
    cities.reserve(200);
    std::ifstream in_file{ filename };
    if (not in_file) {
        utils::printErr("Could not open "s + filename + "\nExiting.",
                        "parseCities");
        exit(EXIT_FAILURE);
    }
    std::string line;
    int line_count{ 0 }, line_expected{ -1 }, line_num;
    Value_t x, y;
    bool start_parse;
    while (std::getline(in_file, line)) {
        if (line.empty()) {
            continue;
        }
        if (line == "DIMENSION") {
            std::string::size_type pos = line.find_first_of(":");
            if (pos != std::string::npos) {
                std::string word{ line.substr(pos + 1) };
                utils::Str2Num number{ word };
                if (number.has_value()) {
                    line_expected = number.value();
                }
            }
            continue;
        }
        line = utils::trimLeft(line);
        std::string word{ line.substr(
            0, line.find_first_of(utils::Whitespace_Str)) };
        start_parse = utils::isNumber(word);
        if (start_parse) {
            ++line_count;
            std::istringstream line_stream{ line };
            line_stream >> line_num >> x >> y;
            City city(-1, x, y);
            cities.push_back(city);
        }
    }
    utils::printInfo(
        "Data parsed from " + filename +
            ".\nTotal number of cities : " + std::to_string(line_count),
        "parseCities");
    std::cout << std::endl;
    if ((line_expected != -1) and (line_count != line_expected)) {
        utils::printErr("Number of coords parsed " +
                            std::to_string(line_count) +
                            " not equal to number of coords available " +
                            std::to_string(line_expected),
                        "parseCities");
    }
}

class MinMaxCoords
{
public:
    MinMaxCoords()
        : min_x{ std::numeric_limits<Value_t>::max() }
        , max_x{ VALUE_ONE_NEG }
        , min_y{ std::numeric_limits<Value_t>::max() }
        , max_y{ VALUE_ONE_NEG }
    {
    }
    MinMaxCoords(Value_t minx, Value_t maxx, Value_t miny, Value_t maxy)
        : min_x{ minx }
        , max_x{ maxx }
        , min_y{ miny }
        , max_y{ maxy }
    {
    }

    auto minmax()
    {
        const Value_t min{ std::min(min_x, min_y) };
        const Value_t max{ std::max(max_x, max_y) };
        return std::make_pair(min, max);
    }

    auto value()
    {
        return std::tie(min_x, max_x, min_y, max_y);
    }

    auto value() const
    {
        return std::tie(min_x, max_x, min_y, max_y);
    }

    void reset()
    {
        min_x = std::numeric_limits<Value_t>::max();
        max_x = VALUE_ONE_NEG;
        min_y = std::numeric_limits<Value_t>::max();
        max_y = VALUE_ONE_NEG;
    }

    void update(const Path_t& path)
    {
        for (const Node_t& node : path) {
            if (node->x < min_x)
                min_x = node->x;
            if (node->y < min_y)
                min_y = node->y;
            if (max_x < node->x)
                max_x = node->x;
            if (max_y < node->y)
                max_y = node->y;
        }
    }

    void update(const Cities_t& cities)
    {
        for (const City& city : cities) {
            if (city.x < min_x)
                min_x = city.x;
            if (city.y < min_y)
                min_y = city.y;
            if (max_x < city.x)
                max_x = city.x;
            if (max_y < city.y)
                max_y = city.y;
        }
    }

private:
    Value_t min_x;
    Value_t max_x;
    Value_t min_y;
    Value_t max_y;
};

auto createCityLayers(const Cities_t& cities, CityLayers_t& city_layers,
                      const MinMaxCoords& minmax_coords, int depth)
{
    if (cities.empty()) {
        return std::make_pair(0, depth);
    }
    const auto [min_x, max_x, min_y, max_y] = minmax_coords.value();
    Value_t mid_x = (min_x + max_x) / 2.0;
    Value_t mid_y = (min_y + max_y) / 2.0;

    constexpr int parts{ 4 };
    Cities_t quadrants[parts];
    for (const City& coord : cities) {
        if (coord.x <= mid_x && coord.y <= mid_y)
            quadrants[0].push_back(coord);
        else if (coord.x > mid_x && coord.y <= mid_y)
            quadrants[1].push_back(coord);
        else if (coord.x <= mid_x && coord.y > mid_y)
            quadrants[2].push_back(coord);
        else if (coord.x > mid_x && coord.y > mid_y)
            quadrants[3].push_back(coord);
    }

    int added{ 0 };
    for (int idx{ 0 }; idx != parts; ++idx) {
        Cities_t& cities_quadrant{ quadrants[idx] };
        if (cities_quadrant.empty()) {
            continue;
        }
        City city = cities_quadrant.back();
        cities_quadrant.pop_back();
        city_layers[depth].push_back(city);
        ++added;
    }
    ++depth;
    {
        auto [tmp_add, tmp_depth] = createCityLayers(
            quadrants[0], city_layers, { min_x, mid_x, min_y, mid_y }, depth);
        added += tmp_add;
        depth = std::max(tmp_depth, depth);
    }
    {
        auto [tmp_add, tmp_depth] = createCityLayers(
            quadrants[1], city_layers, { mid_x, max_x, min_y, mid_y }, depth);
        added += tmp_add;
        depth = std::max(tmp_depth, depth);
    }
    {
        auto [tmp_add, tmp_depth] = createCityLayers(
            quadrants[2], city_layers, { min_x, mid_x, mid_y, max_y }, depth);
        added += tmp_add;
        depth = std::max(tmp_depth, depth);
    }
    {
        auto [tmp_add, tmp_depth] = createCityLayers(
            quadrants[3], city_layers, { mid_x, max_x, mid_y, max_y }, depth);
        added += tmp_add;
        depth = std::max(tmp_depth, depth);
    }
    return std::make_pair(added, depth);
}

void createStack(const Cities_t& cities, Cities_t& stack, int& layers)
{
    const std::size_t num_cities = cities.size();
    stack.reserve(static_cast<std::size_t>(num_cities));
    CityLayers_t city_layers;
    MinMaxCoords minmax_coords;
    minmax_coords.update(cities);
    const auto [min_x, max_x, min_y, max_y] = minmax_coords.value();
    const auto [min_coord, max_coord] = minmax_coords.minmax();
    min_x = min_y = static_cast<Value_t>(0);
    max_x = max_y = max_coord;
    const auto [cities_added, depth] =
        createCityLayers(cities, city_layers, minmax_coords, 1);
    if (static_cast<std::size_t>(cities_added) != num_cities) {
        const std::string& error_msg{ "number of cities in city_layers (" +
                                      std::to_string(cities_added) +
                                      ") doesn't match cities size (" +
                                      std::to_string(num_cities) +
                                      ").\nExiting." };
        utils::printErr(error_msg, "createStack");
        exit(EXIT_FAILURE);
    }
    layers = depth;
    int count{ 0 };
    for (int idx_layer{ 1 }; idx_layer != depth; ++idx_layer) {
        const Cities_t& city_layer{ city_layers[idx_layer] };
        const std::size_t layer_size = city_layer.size();
        for (std::size_t idx{ 0 }; idx != layer_size; ++idx) {
            City city = city_layer[idx];
            city.on_stack = true;
            city.id = count;
            city.layer = idx_layer;
            stack.push_back(city);
            ++count;
        }
    }

    const std::size_t stack_size = stack.size();
    if (stack_size != num_cities) {
        const std::string& error_msg{ "stack size (" +
                                      std::to_string(stack_size) +
                                      ") doesn't match number of cities (" +
                                      std::to_string(num_cities) +
                                      ").\nExiting." };
        utils::printErr(error_msg, "createStack");
        exit(EXIT_FAILURE);
    }
}

inline Value_t getDistance(const City& a, const City& b)
{
    const auto x_sqr = (a.x - b.x) * (a.x - b.x);
    const auto y_sqr = (a.y - b.y) * (a.y - b.y);
    return std::sqrt(x_sqr + y_sqr);
}

inline Value_t insertionCost(const City& new_city, const City& cityA,
                             const City& cityB)
{
    return getDistance(new_city, cityA) + getDistance(new_city, cityB) -
           getDistance(cityA, cityB);
}

inline double getArea(const City& a, const City& b, const City& c)
{
    return std::abs(a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y)) /
           2.0;
}
inline bool isBetween(const City& c, const City& a, const City& b)
{
    return std::min(a.x, b.x) <= c.x && c.x <= std::max(a.x, b.x) &&
           std::min(a.y, b.y) <= c.y && c.y <= std::max(a.y, b.y);
}

inline bool isCollinear(const City& cityA, const City& cityB, const City& cityC)
{
    const Value_t area_total = getArea(cityA, cityB, cityC);

    return utils::isEqual(area_total, VALUE_ZERO);
}

inline bool isInside(const City& city, const City& cityA, const City& cityB,
                     const City& cityC)
{
    const Value_t area_total = getArea(cityA, cityB, cityC);

    const Value_t area1 = getArea(city, cityA, cityB);
    const Value_t area2 = getArea(city, cityB, cityC);
    const Value_t area3 = getArea(city, cityC, cityA);

    if (utils::isEqual(area1, VALUE_ZERO)) {
        return isBetween(city, cityA, cityB);
    }
    if (utils::isEqual(area2, VALUE_ZERO)) {
        return isBetween(city, cityB, cityC);
    }
    if (utils::isEqual(area3, VALUE_ZERO)) {
        return isBetween(city, cityC, cityA);
    }

    return utils::isEqual(area_total, (area1 + area2 + area3));
}

// 0 --> collinear
// 1 --> Clockwise
// 2 --> Counterclockwise
inline int getOrientation(const City& a, const City& b, const City& c)
{
    Value_t val = (b.y - a.y) * (c.x - b.x) - (b.x - a.x) * (c.y - b.y);

    if (utils::isEqual(val, VALUE_ZERO))
        return 0;

    return (val > 0) ? 1 : 2;
}

inline bool hasIntersection(const City& cityA1, const City& cityA2,
                            const City& cityB1, const City& cityB2)
{
    int o1 = getOrientation(cityA1, cityA2, cityB1);
    int o2 = getOrientation(cityA1, cityA2, cityB2);
    int o3 = getOrientation(cityB1, cityB2, cityA1);
    int o4 = getOrientation(cityB1, cityB2, cityA2);

    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    // Special Cases
    if (o1 == 0 && isBetween(cityB1, cityA1, cityA2))
        return true;

    if (o2 == 0 && isBetween(cityB2, cityA1, cityA2))
        return true;

    if (o3 == 0 && isBetween(cityA1, cityB1, cityB2))
        return true;

    if (o4 == 0 && isBetween(cityA2, cityB1, cityB2))
        return true;

    return false; // Doesn't fall in any of the above cases
}

class DiscreteENN_TSP
{
public:
    DiscreteENN_TSP() = default;

    Cities_t& cities()
    {
        return m_cities;
    }
    Indices_t& stack()
    {
        return m_stack;
    }
    Indices_t& path()
    {
        return m_path;
    }
    int& initialSize()
    {
        return m_initialSize;
    }
    bool& rmIntersectRecurse()
    {
        return m_rmIntersectRecurse;
    }
    bool& validIntersectCK()
    {
        return m_validIntersectCK;
    }
    int& iterRandomize()
    {
        return m_iterRandomize;
    }
    int& repeatLength()
    {
        return m_repeatLen;
    }
    Value_t timePerCity()
    {
        return m_timePerIter / m_iterNum;
    }
    auto timePerCityMinMax()
    {
        return std::make_pair(m_timePerIterMin, m_timePerIterMax);
    }

    std::size_t stackBack()
    {
        const std::size_t result{ m_stack.back() };
        m_stack.pop_back();
        m_cities[result].on_stack = false;
        return result;
    }

    std::size_t stackAt(std::size_t index)
    {
        const std::size_t result{ m_stack[index] };
        m_stack.erase(m_stack.begin() + index);
        m_cities[result].on_stack = false;
        return result;
    }

    void removeNode(std::size_t index)
    {
        const std::size_t pos{ m_path[index] };
        const std::string& node_idx_str{ "|" + std::to_string(pos) + "|" };
        std::string::size_type idx_erase = m_pattern.find(node_idx_str);
        m_pattern.erase(idx_erase, node_idx_str.length());

        const Indices_t::iterator node_iter{ m_path.begin() + index };
        m_path.erase(node_iter);
        m_cities[pos].on_stack = true;
        m_stack.push_back(pos);
    }

    void addNode(std::size_t index, std::size_t pos)
    {
        const std::string& node_idx_str{ "|" + std::to_string(m_path[index]) +
                                         "|" };
        const std::string& node_idx_added_str{ "|" + std::to_string(pos) +
                                               "|" };
        std::string::size_type idx_insert = m_pattern.find(node_idx_str);
        m_pattern.insert(idx_insert, node_idx_added_str);

        m_path.insert(m_path.begin() + index, pos);
    }

    std::pair<bool, bool> validateNode(std::size_t index)
    {
        const std::size_t num_nodes = m_path.size();
        if (num_nodes < 3) {
            utils::printErr("given path size less than 3", "validateNode");
            return std::make_pair(false, true);
        }
        const std::size_t pos{ m_path[index] };
        const City& city{ m_cities[pos] };
        const Value_t cost_current{ m_cities[pos].cost };
        if (utils::isEqual(cost_current, VALUE_ZERO)) {
            return std::make_pair(true, false);
        }
        if (index != 0 and index != (num_nodes - 1)) {
            const std::size_t pos{ m_path[num_nodes - 1] };
            const std::size_t pos_next{ m_path[0] };
            if (insertionCost(city, m_cities[pos], m_cities[pos_next]) <
                cost_current) {
                return std::make_pair(false, false);
            }
        }
        for (std::size_t idx{ 0 }; idx != (num_nodes - 1); ++idx) {
            if ((idx == index) or (idx + 1 == index)) {
                continue;
            }
            const std::size_t pos{ m_path[idx] };
            const std::size_t pos_next{ m_path[idx + 1] };
            if (insertionCost(city, m_cities[pos], m_cities[pos_next]) <
                cost_current) {
                return std::make_pair(false, false);
            }
        }
        return std::make_pair(true, false);
    }

    std::size_t properIndex(int index)
    {
        const int n = m_path.size();
        index = (index < 0) ? (n + index) : index;
        index = (index > (n - 1)) ? (index - n) : index;
        return static_cast<size_t>(index);
    }

    auto getNeigbhours(int index)
    {
        return std::make_pair(properIndex(index - 1), properIndex(index + 1));
    }

    Value_t updateCost(std::size_t index)
    {
        const int num_nodes = m_path.size();
        if (num_nodes < 3) {
            utils::printErr("given path size less than 3", "nodeCost");
            return VALUE_ONE_NEG;
        }
        const auto [idx_prev, idx_next] = getNeigbhours(index);
        const std::size_t pos_prev{ m_path[idx_prev] },
            pos_curr{ m_path[index] }, pos_next{ m_path[idx_next] };
        const Value_t cost = isCollinear(m_cities[pos_curr], m_cities[pos_prev],
                                         m_cities[pos_next]) ?
                                 VALUE_ZERO :
                                 insertionCost(m_cities[pos_curr],
                                               m_cities[pos_prev],
                                               m_cities[pos_next]);
        m_cities[pos_curr].cost = cost;
        return cost;
    }

    std::pair<int, bool> updateCostNeighbour(std::size_t index, bool self)
    {
        const std::size_t num_nodes = m_path.size();
        if (num_nodes < 3) {
            utils::printErr("given path size less than 3",
                            "updateCostNeighbour");
            return std::make_pair(-1, true);
        }
        const auto [idx_prev, idx_next] = getNeigbhours(index);
        if (self and updateCost(index) < 0.0) {
            utils::printErr("Invalid node cost calculate for node curr at " +
                                std::to_string(index) + " current path size " +
                                std::to_string(m_path.size()),
                            "updateCostNeighbour");
            utils::printErr("Previous node " + std::to_string(idx_prev),
                            "updateCostNeighbour");
            m_cities[m_path[idx_prev]].print();
            utils::printErr("Current node " + std::to_string(index),
                            "updateCostNeighbour");
            m_cities[m_path[index]].print();
            utils::printErr("Next node " + std::to_string(idx_prev),
                            "updateCostNeighbour");
            m_cities[m_path[idx_next]].print();
            utils::printErr(
                "Distances : " +
                    std::to_string(getDistance(m_cities[m_path[index]],
                                               m_cities[m_path[idx_prev]])) +
                    ", " +
                    std::to_string(getDistance(m_cities[m_path[index]],
                                               m_cities[m_path[idx_next]])) +
                    ", " +
                    std::to_string(getDistance(m_cities[m_path[idx_prev]],
                                               m_cities[m_path[idx_next]])),
                "updateCostNeighbour");
            return std::make_pair(index, true);
        }
        if (updateCost(idx_prev) < 0.0) {
            const auto [idx_prev_tmp, idx_next_tmp] = getNeigbhours(idx_prev);
            utils::printErr("Invalid node cost calculate for node prev at " +
                                std::to_string(idx_prev) +
                                " current path size " +
                                std::to_string(m_path.size()),
                            "updateCostNeighbour");
            utils::printErr("Previous node " + std::to_string(idx_prev_tmp),
                            "updateCostNeighbour");
            m_cities[m_path[idx_prev_tmp]].print();
            utils::printErr("Current node " + std::to_string(idx_prev),
                            "updateCostNeighbour");
            m_cities[m_path[idx_prev]].print();
            utils::printErr("Next node " + std::to_string(idx_next_tmp),
                            "updateCostNeighbour");
            m_cities[m_path[idx_next_tmp]].print();
            utils::printErr(
                "Distances : " +
                    std::to_string(
                        getDistance(m_cities[m_path[idx_prev]],
                                    m_cities[m_path[idx_prev_tmp]])) +
                    ", " +
                    std::to_string(
                        getDistance(m_cities[m_path[idx_prev]],
                                    m_cities[m_path[idx_next_tmp]])) +
                    ", " +
                    std::to_string(getDistance(m_cities[m_path[idx_prev_tmp]],
                                               m_cities[m_path[idx_next_tmp]])),
                "updateCostNeighbour");
            return std::make_pair(idx_prev, true);
        }
        if (updateCost(idx_next) < 0.0) {
            const auto [idx_prev_tmp, idx_next_tmp] = getNeigbhours(idx_next);
            utils::printErr("Invalid node cost calculate for node next at " +
                                std::to_string(idx_next) +
                                " current path size " +
                                std::to_string(m_path.size()),
                            "updateCostNeighbour");
            utils::printErr("Previous node " + std::to_string(idx_prev_tmp),
                            "updateCostNeighbour");
            m_cities[m_path[idx_prev_tmp]].print();
            utils::printErr("Current node " + std::to_string(idx_next),
                            "updateCostNeighbour");
            m_cities[m_path[idx_next]].print();
            utils::printErr("Next node " + std::to_string(idx_next_tmp),
                            "updateCostNeighbour");
            m_cities[m_path[idx_next_tmp]].print();
            utils::printErr(
                "Distances : " +
                    std::to_string(
                        getDistance(m_cities[m_path[idx_next]],
                                    m_cities[m_path[idx_prev_tmp]])) +
                    ", " +
                    std::to_string(
                        getDistance(m_cities[m_path[idx_next]],
                                    m_cities[m_path[idx_next_tmp]])) +
                    ", " +
                    std::to_string(getDistance(m_cities[m_path[idx_prev_tmp]],
                                               m_cities[m_path[idx_next_tmp]])),
                "updateCostNeighbour");
            return std::make_pair(idx_next, true);
        }
        return std::make_pair(index, false);
    }

    std::pair<int, bool> updateCostAll()
    {
        const std::size_t num_nodes = m_path.size();
        for (std::size_t idx{ 0 }; idx != num_nodes; ++idx) {
            if (updateCost(idx) < 0.0) {
                utils::printErr("Invalid node cost calculate for node at " +
                                    std::to_string(idx) +
                                    " current path size " +
                                    std::to_string(num_nodes),
                                "updateCostAll");
                return std::make_pair(idx, true);
            }
        }
        return std::make_pair(num_nodes, false);
    }

    int findBestInsertion(std::size_t position)
    {
        const std::size_t num_nodes = m_path.size();
        if (num_nodes < 2) {
            utils::printErr("given path size less than 2", "findBestInsertion");
            return -1;
        }

        const City& new_city{ m_cities[position] };

        Value_t min_cost{ insertionCost(
            new_city, m_cities[m_path[num_nodes - 1]], m_cities[m_path[0]]) };
        int best_index{ 0 };
        if (new_city.id == m_cities[m_path[num_nodes - 1]].id) {
            utils::printErr("trying to add city that already exists in path",
                            "findBestInsertion");
            return -1;
        }

        for (std::size_t idx{ 0 }; idx != (num_nodes - 1); ++idx) {
            const std::size_t pos{ m_path[idx] };
            const std::size_t pos_next{ m_path[idx + 1] };
            if (new_city.id == m_cities[pos].id) {
                utils::printErr(
                    "trying to add city that already exists in path",
                    "findBestInsertion");
                return -1;
            }
            const Value_t cost =
                insertionCost(new_city, m_cities[pos], m_cities[pos_next]);
            if (cost < min_cost) {
                min_cost = cost;
                best_index = idx + 1;
            }
        }

        return best_index;
    }

    IndexOpt_t removeIntersection(std::size_t pos_prev, std::size_t pos_curr,
                                  std::size_t pos_next)
    {
        IndexOpt_t pos_erased{ std::nullopt };
        int num_nodes = m_path.size();
        City& city_prev{ m_cities[pos_prev] };
        City& city_curr{ m_cities[pos_curr] };
        City& city_next{ m_cities[pos_next] };
        for (int idx{ 0 }; idx != num_nodes;) {
            if (num_nodes == 3) {
                break;
            }
            if (city_prev.on_stack or city_next.on_stack) {
                break;
            }
            const std::size_t pos{ m_path[idx] };
            if ((pos == pos_prev) or (pos == pos_curr) or (pos == pos_next)) {
                ++idx;
                continue;
            }
            if (not isInside(m_cities[pos], city_prev, city_curr, city_next)) {
                ++idx;
                continue;
            }
            removeNode(idx);
            --num_nodes;
            if (pos_erased.has_value()) {
                pos_erased = (pos < (*pos_erased)) ? pos : pos_erased;
            } else {
                pos_erased = pos;
            }
            const std::size_t idx_next = properIndex(idx);
            const std::size_t idx_prev =
                properIndex(static_cast<int>(idx_next) - 1);
            updateCost(idx_next);
            updateCost(idx_prev);
            const std::size_t pos_next = m_path[idx_prev];
            const std::size_t pos_prev = m_path[idx_next];
            if (m_rmIntersectRecurse) {
                auto pos_erased_tmp =
                    removeIntersection(pos_prev, pos, pos_next);
                if (pos_erased_tmp.has_value()) {
                    pos_erased = ((*pos_erased_tmp) < (*pos_erased)) ?
                                     pos_erased_tmp :
                                     pos_erased;
                }
                num_nodes = m_path.size();
                idx = 0;
            }
        }
        return pos_erased;
    }

    IndexExp_t<bool> validatePath()
    {
        IndexOpt_t pos_erased{ std::nullopt };
        std::size_t num_nodes = m_path.size();
        [[maybe_unused]] int num_nodes_prev{ 0 };
        for (std::size_t idx{ 0 }; idx != num_nodes;) {
            if (num_nodes == 2) {
                m_fromScratch = true;
                break;
            }
            const std::size_t pos{ m_path[idx] };
            const auto [valid, err1] = validateNode(idx);
            if (err1) {
                utils::printErr(
                    "validateNode failed at " + std::to_string(idx) +
                        " current path size " + std::to_string(m_path.size()),
                    "validatePath");
                return IndexExp_t<bool>{ pos_erased, true };
            }
            if (valid) {
                ++idx;
                continue;
            }
            removeNode(idx);
            num_nodes_prev = num_nodes;
            --num_nodes;
            if (pos_erased.has_value()) {
                pos_erased = (pos < (*pos_erased)) ? pos : pos_erased;
            } else {
                pos_erased = pos;
            }
            const auto idx_next = properIndex(idx);
            const auto idx_prev = properIndex(static_cast<int>(idx_next) - 1);
            updateCost(idx_next);
            updateCost(idx_prev);
            if (m_validIntersectCK) {
                IndexOpt_t pos_erased_tmp =
                    removeIntersection(m_path[idx_prev], pos, m_path[idx_next]);
                if (pos_erased_tmp.has_value()) {
                    pos_erased = ((*pos_erased_tmp) < (*pos_erased)) ?
                                     pos_erased_tmp :
                                     pos_erased;
                }
                num_nodes = m_path.size();
            }
            idx = 0;
        }
        return IndexExp_t<bool>{ pos_erased, false };
    }

    bool sameCoords(const City& cityA, const City& cityB)
    {
        return std::make_pair(cityA.x, cityA.y) ==
               std::make_pair(cityB.x, cityB.y);
    }
    int checkIntersectEdge(std::size_t start, std::size_t end)
    {
        const std::size_t pos_start{ m_path[start] };
        const std::size_t pos_end{ m_path[end] };
        const City& city_start{ m_cities[pos_start] };
        const City& city_end{ m_cities[pos_end] };
        const std::size_t num_nodes{ m_path.size() };
        const bool ok_start = start != 0 and start != (num_nodes - 1);
        const bool ok_end = end != 0 and end != (num_nodes - 1);
        if (ok_start and ok_end) {
            const std::size_t pos{ m_path[num_nodes - 1] };
            const std::size_t pos_next{ m_path[0] };
            const City& city{ m_cities[pos] };
            const City& city_next{ m_cities[pos_next] };
            if (hasIntersection(city_start, city_end, city, city_next)) {
                if (sameCoords(city_start, city) or
                    sameCoords(city_start, city_next) or
                    sameCoords(city_end, city) or
                    sameCoords(city_end, city_next)) {
                } else {
                    return (num_nodes - 1);
                }
            }
        }
        for (std::size_t idx{ 0 }; idx != (num_nodes - 1); ++idx) {
            const bool ok_start = start != idx and start != (idx + 1);
            const bool ok_end = end != idx and end != (idx + 1);
            if (ok_start and ok_end) {
                const std::size_t pos{ m_path[idx] };
                const std::size_t pos_next{ m_path[idx + 1] };
                const City& city{ m_cities[pos] };
                const City& city_next{ m_cities[pos_next] };
                if (hasIntersection(city_start, city_end, city, city_next)) {
                    if (sameCoords(city_start, city) or
                        sameCoords(city_start, city_next) or
                        sameCoords(city_end, city) or
                        sameCoords(city_end, city_next)) {
                    } else {
                        return idx;
                    }
                }
            }
        }
        return -1;
    }

    bool checkIntersectPath()
    {
        const std::size_t num_nodes{ m_path.size() };
        int intersect0{ checkIntersectEdge(num_nodes - 1, 0) };
        if (intersect0 != -1) {
            const std::size_t intersect_next{ properIndex(intersect0 + 1) };
            utils::printErr("Found an intersection for edge (0, " +
                                std::to_string(num_nodes - 1) + ") and edge (" +
                                std::to_string(intersect0) + ", " +
                                std::to_string(intersect_next) + ")",
                            "checkIntersectPath");
            utils::printInfo("node 0");
            m_cities[m_path[0]].print();
            utils::printInfo("node last(" + std::to_string(num_nodes - 1) +
                             ")");
            m_cities[m_path.back()].print();
            utils::printInfo("node " + std::to_string(intersect0));
            m_cities[m_path[static_cast<std::size_t>(intersect0)]].print();
            utils::printInfo("node " + std::to_string(intersect_next));
            m_cities[m_path[intersect_next]].print();
            return true;
        }
        for (std::size_t idx{ 0 }; idx != (num_nodes - 1); ++idx) {
            int intersect{ checkIntersectEdge(idx, idx + 1) };
            if (intersect != -1) {
                const std::size_t intersect_next{ properIndex(intersect + 1) };
                utils::printInfo("Found an intersection for edge (" +
                                     std::to_string(idx) + ", " +
                                     std::to_string(idx + 1) + ") and edge (" +
                                     std::to_string(intersect) + ", " +
                                     std::to_string(intersect_next) + ")",
                                 "checkIntersectPath");
                utils::printInfo("node " + std::to_string(idx));
                m_cities[m_path[idx]].print();
                utils::printInfo("node " + std::to_string(idx + 1));
                m_cities[m_path[idx + 1]].print();
                utils::printInfo("node " + std::to_string(intersect));
                m_cities[m_path[static_cast<std::size_t>(intersect)]].print();
                utils::printInfo("node " + std::to_string(intersect_next));
                m_cities[m_path[intersect_next]].print();
                return true;
            }
        }
        return false;
    }

    void initializePath()
    {
        const std::size_t num_cities = m_cities.size();
        assert("[Error] (initializePath): given number of cities less than 1" &&
               (num_cities > 0));
        m_path.clear();
        m_path.reserve(num_cities);
        if (num_cities < 4) {
            for (std::size_t idx{ 0 }; idx != num_cities; ++idx) {
                m_path.push_back(idx);
            }
            return;
        }
        m_stack.reserve(num_cities);
        [[maybe_unused]] const std::size_t end{ num_cities - 1 };
        for (std::size_t idx{ 0 }; idx != num_cities; ++idx) {
            m_stack.push_back(end - idx);
            // m_stack.push_back(idx);
        }
    }

    void constructPath()
    {
        const std::size_t num_cities = m_cities.size();
        if (m_initialSize < 1) {
            m_initialSize = 0.10f * num_cities;
            m_initialSize = (m_initialSize == 0) ? 3 : m_initialSize;
        }
        for (int idx{ 0 }; idx != m_initialSize; ++idx) {
            const std::size_t pos{ stackBack() };
            m_path.push_back(pos);
            m_pattern += "|" + std::to_string(pos) + "|";
        }
        updateCostAll();
        assert("[Error] (constructPath): constructed path size less than 3" &&
               (m_path.size() >= 3));
    }

    bool run(std::default_random_engine& gen)
    {
        typedef std::uniform_int_distribution<int> distrib_t;
        [[maybe_unused]] const std::size_t num_cities = m_cities.size();
        distrib_t distrib(0, num_cities - 1);
        std::unordered_set<std::string> pattern_hashes;
        [[maybe_unused]] bool flip = false;
        [[maybe_unused]] bool print_pos{ false };
        std::size_t pos{ stackBack() };
        while (true) {
            if (m_fromScratch) {
                m_fromScratch = false;
                addNode(0, pos);
                const auto [discard1, err1] = updateCostNeighbour(0, true);
                if (err1) {
                    utils::printErr(
                        "updateCostNeighbour failed at index 0 for for path size " +
                            std::to_string(m_path.size()),
                        "run");
                    return false;
                }
                continue;
            }
            // TimePoint_t start_time = std::chrono::steady_clock::now();
            const int idx_added = findBestInsertion(pos);
            if (idx_added == -1) {
                utils::printErr("findBestInsertion failed at index " +
                                    std::to_string(idx_added) +
                                    " for for path size " +
                                    std::to_string(m_path.size()),
                                "run");
                return false;
            }
            addNode(idx_added, pos);
            const auto [discard1, err1] = updateCostNeighbour(idx_added, true);
            if (err1) {
                utils::printErr("updateCostNeighbour failed at index " +
                                    std::to_string(idx_added) +
                                    " for for path size " +
                                    std::to_string(m_path.size()),
                                "run");
                return false;
            }
            // if (print_pos) {
            //     utils::printInfo("Adding node for pos " + std::to_string(pos) +
            //                      " at " + std::to_string(idx_added) +
            //                      " path size " +
            //                      std::to_string(m_path.size()));
            // }
            // if (checkIntersectPath()) {
            //     utils::printErr("intersection after adding node at " +
            //                         std::to_string(idx_added) +
            //                         " current path size " +
            //                         std::to_string(m_path.size()),
            //                     "run");
            //     drawPath(m_path, m_cities, false);
            //     // return false;
            // }

            const auto [idx_prev, idx_next] = getNeigbhours(idx_added);
            const auto it_erased1 =
                removeIntersection(m_path[idx_prev], pos, m_path[idx_next]);
            // if (checkIntersectPath()) {
            //     utils::printErr(
            //         "intersection after removeIntersection from adding node at " +
            //             std::to_string(idx_added) + " current path size " +
            //             std::to_string(m_path.size()),
            //         "run");
            //     drawPath(m_path, m_cities, false);
            //     // return false;
            // }
            // if (print_pos and it_erased1.has_value()) {
            //     utils::printInfo("Adding node for pos " + std::to_string(pos) +
            //                      " at " + std::to_string(idx_added) +
            //                      " and removing intersection node " +
            //                      std::to_string(*it_erased1) + " path size " +
            //                      std::to_string(m_path.size()));
            // }

            const auto it_erased2 = validatePath();
            if (it_erased2.err()) {
                utils::printErr("validatePath failed", "run");
                return false;
            }
            // if (print_pos and it_erased2.has_value()) {
            //     utils::printInfo("Adding node for pos " + std::to_string(pos) +
            //                      " at " + std::to_string(idx_added) +
            //                      " and removing validation node " +
            //                      std::to_string(it_erased2.value()) +
            //                      " path size " + std::to_string(m_path.size()));
            // }
            // if (checkIntersectPath()) {
            //     utils::printErr("intersection after validatePath", "run");
            //     drawPath(m_path, m_cities, false);
            //     return false;
            // }
            if (m_fromScratch) {
                // const int idx_rand{ distrib(gen) };
                // pos = stackAt(idx_rand);
                continue;
            }

            // if (it_erased1.has_value() or it_erased2.has_value()) {
            //     distrib.param(distrib_t::param_type(0, m_stack.size() - 1));
            //     const int idx_rand{ distrib(gen) };
            //     pos = stackAt(idx_rand);
            //     pos = stackBack();
            //     continue;
            // }
            // if (it_erased1.has_value()) {
            //     it = *it_erased1;
            // }
            // if (it_erased2.has_value()) {
            //     if (it_erased1.has_value()) {
            //         it = (*it_erased1)->id < (it_erased2.value())->id ?
            //                  *it_erased1 :
            //                  it_erased2.value();
            //     } else {
            //         it = it_erased2.value();
            //     }
            // }

            const std::size_t stack_size{ m_stack.size() };
            if (stack_size == 0) {
                break;
            }
            if (not pattern_hashes.insert(m_pattern).second) {
                print_pos = true;
                utils::printInfo(
                    "Found repeating pattern. Randomize input node", "run");
                utils::printInfo("Path progress " +
                                     std::to_string(m_path.size()) + "/" +
                                     std::to_string(num_cities),
                                 "run");
                std::cout << std::endl;
                // if (flip) {
                //     distrib.param(distrib_t::param_type(0, m_path.size() - 1));
                // } else {
                //     distrib.param(distrib_t::param_type(m_path.size() - 1,
                //                                         num_cities - 1));
                // }
                // flip = not flip;
                distrib.param(distrib_t::param_type(0, stack_size - 1));
                const int idx_rand{ distrib(gen) };
                utils::printInfo(
                    "New starting point " + std::to_string(idx_rand), "run");
                pos = stackAt(idx_rand);
                continue;
#if (TSP_DEBUG_PRINT > 0)
                std::cout << ("\n[Debug] (run): drawPath started\n");
#endif
                // drawPath(path, stack, true);
#if (TSP_DEBUG_PRINT > 0)
                std::cout << ("[Debug] (run): drawPath ended\n");
#endif
            }
            // print_pos = false;
            // TimePoint_t end_time = std::chrono::steady_clock::now();
            // auto delta = std::chrono::duration_cast<TimeUnit_t>(end_time - start_time);
            // const Value_t duration = delta.count();
            // m_timePerIter += duration;
            // m_timePerIterMin = std::min(duration, m_timePerIterMin);
            // m_timePerIterMax = std::max(duration, m_timePerIterMax);
            // ++m_iterNum;

            // #if (TSP_DEBUG_PRINT > 0)
            //             std::cout << ("\n[Debug] (run): validatePath updateCostAll\n");
            // #endif
            //             const auto [success, idx] = updateCostAll(path);
            //             if (not success) {
            //                 std::cerr
            //                     << "[Error] (run): updateCostAll failed at index "
            //                     << idx << '\n';
            //                 return false;
            //             }

            pos = stackBack();
        }
        return true;
    }

private:
    bool m_rmIntersectRecurse{ Intersection_Recursive };
    bool m_validIntersectCK{ Validation_Intersection };
    bool m_fromScratch{ false };
    int m_initialSize{ Num_Nodes_Initial };
    int m_iterRandomize{ Iter_Randomize };
    int m_repeatLen{ Repeat_Check_Length };
    int m_iterNum{ 0 };
    Value_t m_timePerIter{ VALUE_ZERO };
    Value_t m_timePerIterMin{ std::numeric_limits<Value_t>::max() };
    Value_t m_timePerIterMax{ VALUE_ONE_NEG };
    std::string m_pattern;
    Cities_t m_cities;
    Indices_t m_stack;
    Indices_t m_path;
};

// template <typename T>
// void constructPath(Path_t& path, Cities_t& cities, int n, T rng)
// {
//     const int num_cities = cities.size();
//     std::vector<int> indices(num_cities, 0);
//     std::iota(indices.begin(), indices.end(), 0);
//     std::shuffle(indices.begin(), indices.end(), rng);
//     if (n < 1) {
//         n = 0.10f * num_cities;
//         n = (n == 0) ? 3 : n;
//     }
//     const Node_t it_begin{ cities.begin() };
//     for (int idx{ 0 }; idx != n; ++idx) {
//         Node_t it{ it_begin + indices[idx] };
//         it->on_stack = false;
//         path.push_back(it);
//     }
//     updateCostAll(path);
//     assert("[Error] (constructPath): constructed path size less than 3" &&
//            (path.size() >= 3));
// }

typedef std::vector<sf::Vector2f> VectSF_t;
void path2SFVector(const Indices_t& path, const Cities_t& cities,
                   VectSF_t& vect_sf)
{
    vect_sf.clear();
    for (std::size_t pos : path) {
        const City& city{ cities[pos] };
        vect_sf.push_back(
            { static_cast<float>(city.x), static_cast<float>(city.y) });
    }
}

void fitPointsInWindow(VectSF_t& points, const sf::Vector2u& window_size,
                       const MinMaxCoords& minmax_coords, Value_t margin)
{
    if (points.empty())
        return; // Check if the vector is not empty

    const auto [min_x, max_x, min_y, max_y] = minmax_coords.value();

    Value_t width{ max_x - min_x };
    Value_t height{ max_y - min_y };
    // Calculate scale factors for x and y to fit the plot within the window, considering margins
    Value_t scaleX = (window_size.x - 2 * margin) / width;
    Value_t scaleY = (window_size.y - 2 * margin) / height;

    // Apply scale and translation to the points
    for (auto& point : points) {
        point.x = (point.x - min_x) * scaleX + margin;
        point.y = (point.y - min_y) * scaleY + margin;
    }
}

void fitPointsInWindow(Cities_t& cities, const sf::Vector2u& window_size,
                       const MinMaxCoords& minmax_coords, Value_t margin)
{
    if (cities.empty())
        return; // Check if the vector is not empty

    const auto [min_x, max_x, min_y, max_y] = minmax_coords.value();

    Value_t width{ max_x - min_x };
    Value_t height{ max_y - min_y };
    // Calculate scale factors for x and y to fit the plot within the window, considering margins
    Value_t scaleX = (window_size.x - 2 * margin) / width;
    Value_t scaleY = (window_size.y - 2 * margin) / height;

    // Apply scale and translation to the points
    for (auto& city : cities) {
        city.x = (city.x - min_x) * scaleX + margin;
        city.y = (city.y - min_y) * scaleY + margin;
    }
}

void drawPath(const Indices_t& path, const Cities_t& cities, bool show_coords)
{
    // Get the current desktop video mode
    [[maybe_unused]] sf::VideoMode desktopMode =
        sf::VideoMode::getDesktopMode();
    // Create a full-screen window using the current desktop resolution
    MinMaxCoords minmax_coords;
    minmax_coords.update(cities);
    [[maybe_unused]] const auto [min_coord, max_coord] = minmax_coords.minmax();
    // sf::RenderWindow window(sf::VideoMode(minmax_coords.max_x,
    //                                       minmax_coords.max_y),
    //                         "TSP Route. Press Any Key to Close");
    // sf::RenderWindow window(
    // sf::VideoMode(minmax_coords.max_x, minmax_coords.max_y),
    // "TSP Route. Press Any Key to Close", sf::Style::Fullscreen);
    // sf::RenderWindow window(sf::VideoMode(800, 600), "TSP Route. Press Any Key to Close", sf::Style::Fullscreen);
    sf::RenderWindow window(desktopMode, "TSP Route. Press Any Key to Close",
                            sf::Style::Fullscreen);
    // sf::RenderWindow window(sf::VideoMode(max_coord + 50, max_coord + 50),
    //                         "TSP Route. Press Any Key to Close");

    // Main loop
    VectSF_t points;
    points.reserve(path.size());
    path2SFVector(path, cities, points);
    // fitPointsInWindow(points, {static_cast<unsigned int>(max_coord), static_cast<unsigned int>(max_coord)}, minmax_coords, 20);
    fitPointsInWindow(points, window.getSize(), minmax_coords, 50);
    Cities_t cities_copy = cities;
    // fitPointsInWindow(cities, {static_cast<unsigned int>(max_coord), static_cast<unsigned int>(max_coord)}, minmax_coords, 20);
    fitPointsInWindow(cities_copy, window.getSize(), minmax_coords, 50);
    sf::ConvexShape polygon;
    polygon.setPointCount(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        polygon.setPoint(i, points[i]);
    }

    polygon.setFillColor(sf::Color::Transparent); // Transparent fill
    polygon.setOutlineThickness(2); // Line thickness
    polygon.setOutlineColor(sf::Color::Red); // Line color

    sf::Font font;
    font.loadFromFile("/usr/share/fonts/TTF/Roboto-Regular.ttf");
    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            // Check for window close request
            if (event.type == sf::Event::Closed)
                window.close();

            // Check for any key press
            if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::Escape) {
                    window.close(); // Close the window if Escape key is pressed
                }
                if (event.key.code == sf::Keyboard::Q) {
                    window.close(); // Close the window if Escape key is pressed
                }
                // window.close();
            }
        }

        window.clear(); // Clear the window

        for (const auto& point : points) {
            sf::CircleShape marker(5); // Small circle with radius 5 pixels
            marker.setFillColor(sf::Color::Green); // Use a contrasting color
            marker.setPosition(
                point.x,
                point.y); // Adjust position to center the marker on the vertex

            window.draw(marker);
        }
        for (const City& city : cities_copy) {
            if (show_coords) {
                sf::Text label;
                label.setFont(font); // Set the font
                label.setString(std::to_string(
                    city.id)); // Set the point's ID as the label text
                label.setCharacterSize(14); // Set the character size
                label.setFillColor(sf::Color::White); // Set the text color
                label.setPosition(city.x + 10,
                                  city.y -
                                      10); // Position the label near the marker
                window.draw(label); // Draw the label
            }

            if (not city.on_stack)
                continue;
            sf::CircleShape marker(5); // Small circle with radius 5 pixels
            marker.setFillColor(sf::Color::Blue); // Use a contrasting color
            marker.setPosition(
                city.x,
                city.y); // Adjust position to center the marker on the vertex

            window.draw(marker);
        }
        window.draw(polygon); // Draw the polygon

        window.display(); // Display what was drawn
    }
}

struct TSPInfo
{
    int m_points{ -1 };
    Value_t m_distance{ VALUE_ONE_NEG };
    Value_t m_error{ VALUE_ONE_NEG };
    Value_t m_time{ VALUE_ONE_NEG };
    Value_t m_timePerIter{ VALUE_ONE_NEG };
    Value_t m_timePerCity{ VALUE_ONE_NEG };
    Value_t m_timePerCityMin{ VALUE_ONE_NEG };
    Value_t m_timePerCityMax{ VALUE_ONE_NEG };
    std::string m_name{ "" };
};
template <>
inline bool utils::MatchItem<TSPInfo>::operator()(const TSPInfo& info)
{
    return (info.m_name == m_item.m_name);
}
typedef std::vector<TSPInfo> TSPInfoVect_t;

bool readDistances(const std::string& filename, TSPInfoVect_t& infos)
{
    std::ifstream file{ filename };
    if (not file) {
        utils::printErr("couldn't read file : " + filename, "readDistances");
        return false;
    }
    infos.clear();
    std::string line, word;
    int line_count{ 0 };
    std::string name;
    TSPInfo info;
    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }
        std::stringstream line_stream{ line };
        if (not std::getline(line_stream, word, ',')) {
            continue;
        }
        utils::Str2Num index{ word };
        if (not std::getline(line_stream, word, ',')) {
            continue;
        }
        name = word;
        if (not std::getline(line_stream, word, ',')) {
            continue;
        }
        utils::Str2Num distance{ word };
        if (not index.has_value() or not distance.has_value()) {
            continue;
        }
        ++line_count;
        info.m_name = name;
        info.m_distance = static_cast<Value_t>(distance.value());
        infos.push_back(info);
    }
    utils::printInfo(
        "Distance parsed from " + filename +
            ".\nTotal number of cities : " + std::to_string(line_count),
        "readDistances");
    std::cout << std::endl;

    return true;
}
