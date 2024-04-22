// -*- C++ Header -*-

#pragma once

#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <random>
#include <algorithm>
#include <numeric>
#include <map>
#include <unordered_map>
#include <cmath>
#include <tuple>
#include <utility>
#include <array>

#include "utils.hpp"

#define TSP_DRAW

#ifdef TSP_DRAW
#include <SFML/Graphics.hpp>
#endif

#define TSP_DEBUG_PRINT 0
#define TSP_DEBUG_CHECK 0

using namespace std::string_literals;
using utils::ErrorBool, utils::ErrorMsg;

typedef float Value_t;

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
typedef std::optional<Node_t> NodeOpt_t;
template <typename T> using NodeExp_t = utils::Expected<Node_t, T>;

template <> bool utils::MatchItem<City>::operator()(const City& city)
{
    return (city.id == m_item.id);
}

template <> bool utils::MatchItem<Node_t>::operator()(const Node_t& node)
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
    int line_count{ 0 }, line_num{ 0 };
    Value_t x, y;
    bool start_parse{ false };
    while (std::getline(in_file, line)) {
        if (line.empty())
            continue;
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
    assert(
        "[Error] (parseCities): Number of coords parsed not equal to number of coords available" &&
        (line_count == line_num));
}

class MinMaxCoords
{
public:
    MinMaxCoords()
        : min_x{ std::numeric_limits<Value_t>::max() }
        , max_x{ -1.0 }
        , min_y{ std::numeric_limits<Value_t>::max() }
        , max_y{ -1.0 }
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
        max_x = -1.0;
        min_y = std::numeric_limits<Value_t>::max();
        max_y = -1.0;
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
    const int num_cities = cities.size();
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
    if (cities_added != num_cities) {
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

    const int stack_size = stack.size();
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

Value_t getDistance(const City& a, const City& b)
{
    const auto x_sqr = (a.x - b.x) * (a.x - b.x);
    const auto y_sqr = (a.y - b.y) * (a.y - b.y);
    return std::sqrt(x_sqr + y_sqr);
}

Value_t insertionCost(const City& new_city, const City& cityA,
                      const City& cityB)
{
    return getDistance(new_city, cityA) + getDistance(new_city, cityB) -
           getDistance(cityA, cityB);
}

double getArea(const City& a, const City& b, const City& c)
{
    return std::abs(a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y)) /
           2.0;
}
bool isCollinearAndBetween(const City& c, const City& a, const City& b)
{
    return std::min(a.x, b.x) <= c.x && c.x <= std::max(a.x, b.x) &&
           std::min(a.y, b.y) <= c.y && c.y <= std::max(a.y, b.y);
}

bool isCollinear(Node_t nodeA, Node_t nodeB, Node_t nodeC)
{
    const Value_t area_total = getArea(*nodeA, *nodeB, *nodeC);

    return (std::abs(area_total) < utils::tolerance);
}

bool isInside(Node_t node, Node_t nodeA, Node_t nodeB, Node_t nodeC)
{
    const Value_t area_total = getArea(*nodeA, *nodeB, *nodeC);

    const Value_t area1 = getArea(*node, *nodeA, *nodeB);
    const Value_t area2 = getArea(*node, *nodeB, *nodeC);
    const Value_t area3 = getArea(*node, *nodeC, *nodeA);

    if (std::abs(area1) < utils::tolerance) {
        return isCollinearAndBetween(*node, *nodeA, *nodeB);
    }
    if (std::abs(area2) < utils::tolerance) {
        return isCollinearAndBetween(*node, *nodeB, *nodeC);
    }
    if (std::abs(area3) < utils::tolerance) {
        return isCollinearAndBetween(*node, *nodeC, *nodeA);
    }

    const Value_t area_diff = std::abs(area_total - (area1 + area2 + area3));
    return (area_diff < utils::tolerance);
}

class DiscreteENN_TSP
{
public:
    DiscreteENN_TSP() = default;

    Cities_t& stack()
    {
        return m_stack;
    }
    Path_t& path()
    {
        return m_path;
    }
    int& initialSize()
    {
        return m_initialSize;
    }
    bool& recursive()
    {
        return m_recursive;
    }
    bool& intersection()
    {
        return m_intersection;
    }
    int& iterRandomize()
    {
        return m_iterRandomize;
    }
    int& repeatLength()
    {
        return m_repeatLen;
    }

    int properIndex(int index)
    {
        const int n = m_path.size();
        index = (index < 0) ? (n + index) : index;
        index = (index > (n - 1)) ? (index - n) : index;
        return index;
    }

    int findNode(Node_t node)
    {
        const Path_t::const_iterator& it_begin{ m_path.begin() };
        const Path_t::const_iterator& it_end{ m_path.end() };
        const Path_t::const_iterator& it =
            std::find_if(it_begin, it_end, utils::MatchItem<Node_t>{ node });
        return (it == it_end) ? -1 : std::distance(it_begin, it);
    }

    bool removeNode(int index)
    {
#if (TSP_DEBUG_CHECK > 1)
        if (index >= static_cast<int>(m_path.size())) {
            utils::printErr("Request to remove a node not present in path.",
                            "removeNode");
            return false;
        }
        if (index < 0) {
            utils::printErr("given index is negative", "removeNode");
            return false;
        }
#endif
        const Path_t::iterator node_iter{ m_path.begin() + index };
        (*node_iter)->on_stack = true;
#if (TSP_DEBUG_PRINT > 0)
        (*node_iter)->print();
#endif
        m_path.erase(node_iter);
        return true;
    }

    bool removeNode(Node_t node)
    {
        const int index = findNode(node);
        if (index == -1) {
            utils::printErr("Request to remove a node not present in path.",
                            "removeNode");
            return false;
        }
        return removeNode(index);
    }

    void addNode(int index, Node_t node)
    {
        node->on_stack = false;
        m_path.insert(m_path.begin() + index, node);
#if (TSP_DEBUG_PRINT > 1)
        m_path[0]->print();
        std::cout << "FIRST CITY\n";
        node->print();
        std::cout << "Added City\n";
        const std::string& msg{ "path size (" + std::to_string(m_path.size()) +
                                ") position : " + std::to_string(index) };
        std::cout << ("[Debug] (addNode): " + msg + "\n");
#endif
    }

    std::pair<bool, bool> validateNode(Node_t node)
    {
        const int index = findNode(node);
        if (index == -1) {
            utils::printErr("Request to validate a node not present in path.",
                            "validateNode");
            return std::make_pair(false, true);
        }
        const std::size_t num_nodes = m_path.size();
        if (num_nodes < 3) {
            utils::printErr("given path size less than 3", "validateNode");
            return std::make_pair(false, true);
        }
        const Value_t cost_current{ node->cost };
        if (index != 0 and index != static_cast<int>(num_nodes - 1)) {
            if (insertionCost(*node, *(m_path[num_nodes - 1]), *(m_path[0])) <
                cost_current) {
                return std::make_pair(false, false);
            }
        }
        for (std::size_t idx{ 0 }; idx != (num_nodes - 1); ++idx) {
            if ((idx == static_cast<std::size_t>(index)) or
                (idx + 1 == static_cast<std::size_t>(index))) {
                continue;
            }
            if (insertionCost(*node, *(m_path[idx]), *(m_path[idx + 1])) <
                cost_current) {
                return std::make_pair(false, false);
            }
        }
        return std::make_pair(true, false);
    }

    auto getNeigbhours(int index)
    {
        const int curr = properIndex(index);
        return std::make_tuple(properIndex(curr - 1), curr,
                               properIndex(curr + 1));
    }

    auto getNeigbhours(Node_t node)
    {
        const int index = findNode(node);
        if (index == -1) {
            utils::printErr(
                "Request to find neighbours for a node not present in path.",
                "getNeigbhours");
            utils::Expected{ std::make_tuple(Node_t{}, Node_t{}, Node_t{}),
                             true };
        }
        const auto [prev, curr, next] = getNeigbhours(index);
        return utils::Expected{ std::make_tuple(
                                    m_path[static_cast<std::size_t>(prev)],
                                    m_path[static_cast<std::size_t>(curr)],
                                    m_path[static_cast<std::size_t>(next)]),
                                false };
    }

    Value_t nodeCost(int index, bool update)
    {
        const int num_nodes = m_path.size();
        if (num_nodes < 3) {
            utils::printErr("given path size less than 3", "nodeCost");
            return -1.0;
        }
        const auto [idx_prev, idx_curr, idx_next] = getNeigbhours(index);
        const auto node_prev{ m_path[idx_prev] }, node_curr{ m_path[index] },
            node_next{ m_path[idx_next] };
        const Value_t cost =
            isCollinear(node_curr, node_prev, node_next) ?
                Value_t{0.0} :
                insertionCost(*node_curr, *node_prev, *node_next);
        if (update) {
            node_curr->cost = cost;
        }
        return cost;
    }
    Value_t nodeCost(Node_t node, bool update)
    {
        const int index = findNode(node);
        if (index == -1) {
            utils::printErr(
                "Request to calculate cost for a node not present in path. Current path size " +
                    std::to_string(m_path.size()),
                "nodeCost");
            return -1.0;
        }
        return nodeCost(index, update);
    }

    NodeOpt_t removeIntersection(Node_t node_prev, Node_t node_curr,
                                 Node_t node_next)
    {
        NodeOpt_t node_erased{ std::nullopt };
        int num_nodes = m_path.size();
        const int id_prev{ node_prev->id };
        const int id_curr{ node_curr->id };
        const int id_next{ node_next->id };
        for (int idx{ 0 }; idx != num_nodes;) {
#if (TSP_DEBUG_PRINT > 1)
            const std::string& msg_first{
                "idx (" + std::to_string(idx) + ") num_nodes (" +
                std::to_string(num_nodes) + ") path size (" +
                std::to_string(m_path.size()) + ")."
            };
            std::cout << ("[Debug] (removeIntersection): " + msg_first + "\n");
#endif
            const Node_t node{ m_path[static_cast<std::size_t>(idx)] };
            const int id{ node->id };
            if ((id == id_prev) or (id == id_curr) or (id == id_next)) {
                ++idx;
                continue;
            }
            if (not isInside(node, node_prev, node_curr, node_next)) {
                ++idx;
                continue;
            }
#if (TSP_DEBUG_PRINT > 1)
            std::cout << ("[Debug] (removeIntersection): Done 1\n");
#endif
            removeNode(node);
            --num_nodes;
#if (TSP_DEBUG_PRINT > 1)
            std::cout << ("[Debug] (removeIntersection): Done 2\n");
#endif
            if (node_erased.has_value()) {
                node_erased = (node->id < (*node_erased)->id) ? node :
                                                                node_erased;
            } else {
                node_erased = node;
            }
#if (TSP_DEBUG_PRINT > 0)
            path[0]->print();
            std::cout << "FIRST City\n";
            node->print();
            std::cout << "Removed City\n";
            node_prev->print();
            std::cout << "Triangle A City\n";
            node_curr->print();
            std::cout << "Triangle B City\n";
            node_next->print();
            std::cout << "Triangle C City\n";
            const std::string& msg{ "path size (" +
                                    std::to_string(m_path.size()) +
                                    ") num_nodes (" +
                                    std::to_string(num_nodes) + ") idx (" +
                                    std::to_string(idx) + ")." };
            std::cout << ("[Debug] (removeIntersection): " + msg + "\n");
#endif
            const auto idx_next = static_cast<std::size_t>(properIndex(idx));
            const auto idx_prev =
                static_cast<std::size_t>(properIndex(idx_next - 1));
            nodeCost(idx_next, true);
            nodeCost(idx_prev, true);
            if (m_recursive) {
                auto node_erased_tmp = removeIntersection(
                    m_path[idx_prev], node, m_path[idx_next]);
                if (node_erased_tmp.has_value()) {
                    node_erased =
                        ((*node_erased_tmp)->id < (*node_erased)->id) ?
                            node_erased_tmp :
                            node_erased;
                }
                num_nodes = m_path.size();
                idx = 0;
            }
        }
        return node_erased;
    }

    int findBestInsertion(Node_t node)
    {
        const std::size_t num_nodes = m_path.size();
        if (num_nodes < 2) {
            utils::printErr("given path size less than 2", "findBestInsertion");
            return -1;
        }

        const City& new_city{ *node };

        Value_t min_cost{ insertionCost(new_city, *(m_path[num_nodes - 1]),
                                        *(m_path[0])) };
        int best_index{ 0 };

        for (std::size_t idx{ 0 }; idx != num_nodes; ++idx) {
            if (new_city.id == (m_path[idx])->id) {
                utils::printErr(
                    "trying to add city that already exists in path",
                    "findBestInsertion");
                return -1;
            }
            if (idx == (num_nodes - 1)) {
                continue;
            }
            const Value_t cost =
                insertionCost(new_city, *(m_path[idx]), *(m_path[idx + 1]));
            if (cost < min_cost) {
                min_cost = cost;
                best_index = idx + 1;
            }
        }

        return best_index;
    }

    std::pair<int, bool> updateCostNeighbour(Node_t node)
    {
        const std::size_t num_nodes = m_path.size();
        if (num_nodes < 3) {
            utils::printErr("given path size less than 3",
                            "updateCostNeighbour");
            return std::make_pair(-1, true);
        }
        const int index = findNode(node);
        if (index == -1) {
            utils::printErr(
                "Request to update neighbour costs for a node not present in path. Current path size " +
                    std::to_string(m_path.size()),
                "updateCostNeighbour");
            return std::make_pair(-1, true);
        }
        const utils::Expected nodes = getNeigbhours(node);
        if (nodes.err()) {
            utils::printErr("Failed to get the neighbour nodes for index " +
                                std::to_string(index) + " current path size " +
                                std::to_string(m_path.size()),
                            "updateCostNeighbour");
            return std::make_pair(-1, true);
        }
        const auto [node_prev, node_curr, node_next] = nodes.value();
        if (node != node_curr) {
            utils::printErr(
                "The node fetched from getNeigbhours doesn't match original node for index " +
                    std::to_string(index) + " current path size " +
                    std::to_string(m_path.size()),
                "updateCostNeighbour");
            return std::make_pair(-1, true);
        }
        if (nodeCost(node_curr, true) < 0.0) {
            utils::printErr("Invalid node cost calculate for node curr at " +
                                std::to_string(index) + " current path size " +
                                std::to_string(m_path.size()),
                            "updateCostNeighbour");
            utils::printErr("Previous node " +
                                std::to_string(findNode(node_prev)),
                            "updateCostNeighbour");
            node_prev->print();
            utils::printErr("Current node " +
                                std::to_string(findNode(node_curr)),
                            "updateCostNeighbour");
            node_curr->print();
            utils::printErr("Next node " + std::to_string(findNode(node_next)),
                            "updateCostNeighbour");
            node_next->print();
            utils::printErr(
                "Distances : " +
                    std::to_string(getDistance(*node_curr, *node_next)) + ", " +
                    std::to_string(getDistance(*node_curr, *node_prev)) + ", " +
                    std::to_string(getDistance(*node_prev, *node_next)),
                "updateCostNeighbour");
            return std::make_pair(index, true);
        }
        if (nodeCost(node_prev, true) < 0.0) {
            const utils::Expected nodes = getNeigbhours(node_prev);
            if (nodes.err()) {
                utils::printErr("Failed to get the neighbour nodes for index " +
                                    std::to_string(findNode(node_prev)) +
                                    " current path size " +
                                    std::to_string(m_path.size()),
                                "updateCostNeighbour");
            }
            const auto [node_prev, node_curr, node_next] = nodes.value();
            utils::printErr("Invalid node cost calculate for node prev at " +
                                std::to_string(properIndex(index - 1)) +
                                " current path size " +
                                std::to_string(m_path.size()),
                            "updateCostNeighbour");
            utils::printErr("Previous node " +
                                std::to_string(findNode(node_prev)),
                            "updateCostNeighbour");
            node_prev->print();
            utils::printErr("Current node " +
                                std::to_string(findNode(node_curr)),
                            "updateCostNeighbour");
            node_curr->print();
            utils::printErr("Next node " + std::to_string(findNode(node_next)),
                            "updateCostNeighbour");
            node_next->print();
            utils::printErr(
                "Distances : " +
                    std::to_string(getDistance(*node_curr, *node_next)) + ", " +
                    std::to_string(getDistance(*node_curr, *node_prev)) + ", " +
                    std::to_string(getDistance(*node_prev, *node_next)),
                "updateCostNeighbour");
            return std::make_pair(properIndex(index - 1), true);
        }
        if (nodeCost(node_next, true) < 0.0) {
            const utils::Expected nodes = getNeigbhours(node_next);
            if (nodes.err()) {
                utils::printErr("Failed to get the neighbour nodes for index " +
                                    std::to_string(findNode(node_next)) +
                                    " current path size " +
                                    std::to_string(m_path.size()),
                                "updateCostNeighbour");
            }
            const auto [node_prev, node_curr, node_next] = nodes.value();
            utils::printErr("Invalid node cost calculate for node next at " +
                                std::to_string(properIndex(index + 1)) +
                                " current path size " +
                                std::to_string(m_path.size()),
                            "updateCostNeighbour");
            utils::printErr("Previous node " +
                                std::to_string(findNode(node_prev)),
                            "updateCostNeighbour");
            node_prev->print();
            utils::printErr("Current node " +
                                std::to_string(findNode(node_curr)),
                            "updateCostNeighbour");
            node_curr->print();
            utils::printErr("Next node " + std::to_string(findNode(node_next)),
                            "updateCostNeighbour");
            node_next->print();
            utils::printErr(
                "Distances : " +
                    std::to_string(getDistance(*node_curr, *node_next)) + ", " +
                    std::to_string(getDistance(*node_curr, *node_prev)) + ", " +
                    std::to_string(getDistance(*node_prev, *node_next)),
                "updateCostNeighbour");
            return std::make_pair(properIndex(index + 1), true);
        }
        return std::make_pair(index, false);
    }

    std::pair<int, bool> updateCostAll()
    {
        const int num_nodes = m_path.size();
        for (int idx{ 0 }; idx != num_nodes; ++idx) {
            if (nodeCost(m_path[static_cast<std::size_t>(idx)], true) < 0.0) {
                utils::printErr("Invalid node cost calculate for node at " +
                                    std::to_string(idx) +
                                    " current path size " +
                                    std::to_string(m_path.size()),
                                "updateCostAll");
                return std::make_pair(idx, true);
            }
        }
        return std::make_pair(num_nodes, false);
    }

    NodeExp_t<bool> validatePath()
    {
        NodeOpt_t node_erased{ std::nullopt };
        int num_nodes = m_path.size();
        int num_nodes_prev{ 0 };
        for (int idx{ 0 }; idx != num_nodes;) {
#if (TSP_DEBUG_PRINT > 1)
            const std::string& msg{ "idx (" + std::to_string(idx) +
                                    ") num_nodes (" +
                                    std::to_string(num_nodes) +
                                    ") path size (" +
                                    std::to_string(m_path.size()) + ")." };
            std::cout << ("[Debug] (validatePath): " + msg + "\n");
#endif
            const Node_t node{ m_path[static_cast<std::size_t>(idx)] };
            const auto [valid, err1] = validateNode(node);
            if (err1) {
                utils::printErr(
                    "validateNode failed at " + std::to_string(idx) +
                        " current path size " + std::to_string(m_path.size()),
                    "validatePath");
                return NodeExp_t<bool>{ node_erased, true };
            }
            if (valid) {
                ++idx;
                continue;
            }
#if (TSP_DEBUG_PRINT > 1)
            std::cout << ("[Debug] (validatePath): Done 1\n");
#endif
            removeNode(node);
            num_nodes_prev = num_nodes;
            --num_nodes;
            if (node_erased.has_value()) {
                node_erased = (node->id < (*node_erased)->id) ? node :
                                                                node_erased;
            } else {
                node_erased = node;
            }
#if (TSP_DEBUG_PRINT > 1)
            std::cout << ("[Debug] (validatePath): Done 2\n");
#endif
            const auto idx_next = static_cast<std::size_t>(properIndex(idx));
            const auto idx_prev =
                static_cast<std::size_t>(properIndex(idx_next - 1));
            nodeCost(idx_next, true);
            nodeCost(idx_prev, true);
            if (m_intersection) {
                NodeOpt_t node_erased_tmp = removeIntersection(
                    m_path[idx_prev], node, m_path[idx_next]);
                if (node_erased_tmp.has_value()) {
                    node_erased =
                        ((*node_erased_tmp)->id < (*node_erased)->id) ?
                            node_erased_tmp :
                            node_erased;
                }
                num_nodes = m_path.size();
            }
#if (TSP_DEBUG_PRINT > 1)
            std::cout << ("[Debug] (validatePath): Done 3\n");
#endif
            idx = 0;
#if (TSP_DEBUG_PRINT > 1)
            std::cout << ("[Debug] (validatePath): Done 4\n");
#endif
        }
        return NodeExp_t<bool>{ node_erased, false };
    }

    void initializePath()
    {
        const std::size_t num_cities = m_stack.size();
        assert("[Error] (initializePath): given number of cities less than 1" &&
               (num_cities > 0));
        m_path.clear();
        m_path.reserve(num_cities);
        const Node_t it_begin{ m_stack.begin() };
        if (num_cities < 4) {
            for (Node_t it{ it_begin }; it != m_stack.end(); ++it) {
                m_path.push_back(it);
            }
        }
    }

    void constructPath()
    {
        const int num_cities = m_stack.size();
        if (m_initialSize < 1) {
            m_initialSize = 0.10f * num_cities;
            m_initialSize = (m_initialSize == 0) ? 3 : m_initialSize;
        }
        const Node_t it_begin{ m_stack.begin() };
        for (int idx{ 0 }; idx != m_initialSize; ++idx) {
            Node_t it{ it_begin + idx };
            it->on_stack = false;
            m_path.push_back(it);
        }
        updateCostAll();
        assert("[Error] (constructPath): constructed path size less than 3" &&
               (m_path.size() >= 3));
    }

    bool run(std::default_random_engine& gen)
    {
        [[maybe_unused]] const int num_cities = m_stack.size();
        std::uniform_int_distribution<int> distrib(0, num_cities - 1);
        const Node_t& it_begin{ m_stack.begin() };
        const Node_t& it_end{ m_stack.end() };
        [[maybe_unused]] std::size_t repeat_check{ 0 };
        std::vector<int> indices_added(
            static_cast<std::size_t>(m_iterRandomize));
        for (Node_t it{ it_begin }; it != it_end;) {
            if (not it->on_stack) {
                ++it;
            } else {
#if (TSP_DEBUG_PRINT > 0)
                std::cout << ("\n[Debug] (run): addNode\n");
#endif
                const int idx_added = findBestInsertion(it);
                if (idx_added == -1) {
                    utils::printErr("findBestInsertion failed at index " +
                                        std::to_string(idx_added) +
                                        " for for path size " +
                                        std::to_string(m_path.size()),
                                    "run");
                    return false;
                }
                addNode(idx_added, it);
                const auto [discard1, err1] = updateCostNeighbour(it);
                if (err1) {
                    utils::printErr("updateCostNeighbour failed at index " +
                                        std::to_string(idx_added) +
                                        " for for path size " +
                                        std::to_string(m_path.size()),
                                    "run");
                    return false;
                }
                indices_added[repeat_check] = idx_added;
                ++repeat_check;
#if (TSP_DEBUG_PRINT > 0)
                std::cout << ("\n[Debug] (run): addNode removeIntersection\n");
#endif

                const utils::Expected nodes = getNeigbhours(it);
                if (nodes.err()) {
                    utils::printErr(
                        "Failed to get the neighbour nodes for index " +
                            std::to_string(idx_added) + " current path size " +
                            std::to_string(m_path.size()),
                        "run");
                    return false;
                }
                const auto [node_prev, node_curr, node_next] = nodes.value();
                if (it != node_curr) {
                    utils::printErr(
                        "The node fetched from getNeigbhours doesn't match original node for index " +
                            std::to_string(idx_added) + " current path size " +
                            std::to_string(m_path.size()),
                        "run");
                    return false;
                }
                const auto it_erased1 =
                    removeIntersection(node_prev, node_curr, node_next);
                if (it_erased1.has_value()) {
#if (TSP_DEBUG_PRINT > 0)
                    std::cout << ("\n[Debug] (run): addNode updateCostAll\n");
#endif
                    // const auto [idx_fail, err] = updateCostAll();
                    // if (err) {
                    //     utils::printErr("updateCostAll failed at index " +
                    //                         std::to_string(idx_fail) +
                    //                         " path size " +
                    //                         std::to_string(m_path.size()),
                    //                     "run");
                    //     return false;
                    // }
                }
#if (TSP_DEBUG_PRINT > 0)
                std::cout << ("\n[Debug] (run): validatePath\n");
#endif
                const auto it_erased2 = validatePath();
                if (it_erased2.err()) {
                    utils::printErr("validatePath failed", "run");
                    return false;
                }

                // if (it_erased1.has_value() or it_erased2.has_value()) {
                //     it = it_begin;
                // }
                if (it_erased1.has_value()) {
                    it = *it_erased1;
                }
                if (it_erased2.has_value()) {
                    if (it_erased1.has_value()) {
                        it = (*it_erased1)->id < (it_erased2.value())->id ?
                                 *it_erased1 :
                                 it_erased2.value();
                    } else {
                        it = it_erased2.value();
                    }
                }

                if (repeat_check >= static_cast<std::size_t>(m_iterRandomize)) {
                    // utils::printInfo("Many indices added : " +
                    //                      std::to_string(indices_added.size()),
                    //                  "run");
                    repeat_check = 0;
                    const bool repeating{ utils::hasRepeatingPattern(
                        indices_added, m_repeatLen) };
                    if (repeating) {
                        utils::printInfo(
                            "Found repeating pattern. Randomize input node",
                            "run");
                        utils::printInfo("Path progress " +
                                             std::to_string(m_path.size()) +
                                             "/" + std::to_string(num_cities),
                                         "run");
                        utils::printInfo("Repeating pattern:",
                                         "run");
                        const auto idx_it_end{ indices_added.end() };
                        for (auto idx_it{idx_it_end - m_repeatLen}; idx_it != idx_it_end; ++idx_it) {
                            std::cout << *idx_it << ",";

                        }
                        std::cout << std::endl;
                        const int idx_rand{ distrib(gen) };
                        utils::printInfo("New starting point " +
                                             std::to_string(idx_rand),
                                         "run");
                        it = it_begin + idx_rand;
                        continue;
#if (TSP_DEBUG_PRINT > 0)
                        std::cout << ("\n[Debug] (run): drawPath started\n");
#endif
                        // drawPath(path, stack, true);
#if (TSP_DEBUG_PRINT > 0)
                        std::cout << ("[Debug] (run): drawPath ended\n");
#endif
                    }
                }
            }

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

            if (it == it_end) {
                const bool finished = m_path.size() ==
                                      static_cast<std::size_t>(num_cities);
                if (not finished) {
                    it = it_begin;
                }
            }
        }
        return true;
    }

private:
    bool m_recursive;
    bool m_intersection;
    int m_initialSize;
    int m_iterRandomize;
    int m_repeatLen;
    Cities_t m_stack;
    Path_t m_path;
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
void path2SFVector(const Path_t& path, VectSF_t& vect_sf)
{
    vect_sf.clear();
    for (const auto& node : path) {
        vect_sf.push_back(
            { static_cast<float>(node->x), static_cast<float>(node->y) });
    }
}

void fitPointsInWindow(VectSF_t& points, const sf::Vector2u& windowSize,
                       const MinMaxCoords& minmax_coords, Value_t margin)
{
    if (points.empty())
        return; // Check if the vector is not empty

    const auto [min_x, max_x, min_y, max_y] = minmax_coords.value();

    // Calculate scale factors for x and y to fit the plot within the window, considering margins
    Value_t scaleX = (windowSize.x - 2 * margin) / (max_x - min_x);
    Value_t scaleY = (windowSize.y - 2 * margin) / (max_y - min_y);
    Value_t scale = std::min(
        scaleX,
        scaleY); // Use the smaller scale factor to maintain aspect ratio

    // Calculate translation to center the plot
    Value_t translateX = (windowSize.x - (max_x + min_x) * scale) / 2.0f;
    Value_t translateY = (windowSize.y - (max_y + min_y) * scale) / 2.0f;

    // Apply scale and translation to the points
    for (auto& point : points) {
        point.x = (point.x - min_x) * scale + margin + translateX;
        point.y = (point.y - min_y) * scale + margin + translateY;
    }
}

void fitPointsInWindow(Cities_t& cities, const sf::Vector2u& windowSize,
                       const MinMaxCoords& minmax_coords, Value_t margin)
{
    if (cities.empty())
        return; // Check if the vector is not empty

    const auto [min_x, max_x, min_y, max_y] = minmax_coords.value();

    // Calculate scale factors for x and y to fit the plot within the window, considering margins
    Value_t scaleX = (windowSize.x - 2 * margin) / (max_x - min_x);
    Value_t scaleY = (windowSize.y - 2 * margin) / (max_y - min_y);
    Value_t scale = std::min(
        scaleX,
        scaleY); // Use the smaller scale factor to maintain aspect ratio

    // Calculate translation to center the plot
    Value_t translateX = (windowSize.x - (max_x + min_x) * scale) / 2.0f;
    Value_t translateY = (windowSize.y - (max_y + min_y) * scale) / 2.0f;

    // Apply scale and translation to the points
    for (auto& city : cities) {
        city.x = (city.x - min_x) * scale + margin + translateX;
        city.y = (city.y - min_y) * scale + margin + translateY;
    }
}

void drawPath(const Path_t& path, const Cities_t& stack, bool show_coords)
{
    // Get the current desktop video mode
    [[maybe_unused]] sf::VideoMode desktopMode =
        sf::VideoMode::getDesktopMode();
    // Create a full-screen window using the current desktop resolution
    MinMaxCoords minmax_coords;
    minmax_coords.update(stack);
    const auto [min_coord, max_coord] = minmax_coords.minmax();
    // sf::RenderWindow window(sf::VideoMode(minmax_coords.max_x,
    //                                       minmax_coords.max_y),
    //                         "Polygon Plot. Press Any Key to Close");
    // sf::RenderWindow window(
    // sf::VideoMode(minmax_coords.max_x, minmax_coords.max_y),
    // "Polygon Plot. Press Any Key to Close", sf::Style::Fullscreen);
    // sf::RenderWindow window(sf::VideoMode(800, 600), "Polygon Plot. Press Any Key to Close", sf::Style::Fullscreen);
    // sf::RenderWindow window(desktopMode, "Polygon Plot. Press Any Key to Close",
    //                         sf::Style::Fullscreen);
    sf::RenderWindow window(sf::VideoMode(max_coord, max_coord),
                            "Polygon Plot. Press Any Key to Close");

    // Main loop
    VectSF_t points;
    points.reserve(path.size());
    path2SFVector(path, points);
    fitPointsInWindow(points, window.getSize(), minmax_coords, 20);
    Cities_t cities = stack;
    fitPointsInWindow(cities, window.getSize(), minmax_coords, 20);
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
        for (const City& city : cities) {
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

typedef std::unordered_map<std::string, std::array<Value_t, 4>> DistanceMap_t;
bool readDistances(const std::string& filename, DistanceMap_t& distance_map)
{
    std::ifstream file{ filename };
    if (not file) {
        utils::printErr("couldn't read file : " + filename, "readDistances");
        return false;
    }
    distance_map.clear();
    std::string line, word;
    int line_count{ 0 }, line_num{ 0 };
    std::string name;
    while (std::getline(file, line)) {
        if (line.empty())
            continue;
        std::stringstream line_stream{ line };
        if (not std::getline(line_stream, word, ',')) {
            continue;
            ;
        }
        utils::Str2Num index{ word };
        if (not std::getline(line_stream, word, ',')) {
            continue;
            ;
        }
        name = word;
        if (not std::getline(line_stream, word, ',')) {
            continue;
            ;
        }
        utils::Str2Num distance{ word };
        if (not index.has_value() or not distance.has_value()) {
            continue;
            ;
        }
        ++line_count;
        line_num = index.value();
        distance_map[name] =
            DistanceMap_t::mapped_type{ static_cast<float>(distance.value()) };
    }
    utils::printInfo(
        "Distance parsed from " + filename +
            ".\nTotal number of cities : " + std::to_string(line_count),
        "readDistances");
    std::cout << std::endl;
    assert(
        "[Error] (readDistances): Number of distances parsed not equal to number of distances available" &&
        (line_count == line_num));

    return true;
}
