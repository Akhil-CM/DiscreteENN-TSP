// -*- C++ Header -*-

#pragma once

#define TSP_DRAW

#define TSP_DEBUG_PRINT 2
#define TSP_DRAW_DISPLAY 0
#define TSP_DRAW_SAVE 0
#define TSP_DRAW_ROOT 1
#define TSP_ERROR 0

#define TSP_TIME 2

#include <cstddef>
#include <limits>
#include <cstdlib>
#include <set>
#include <random>
#include <unordered_map>
#include <string>
#include <unordered_set>
#include <cmath>
#include <tuple>
#include <utility>
#include <chrono>

#if TSP_DRAW_DISPLAY > 0
#include <thread>
#endif

#if TSP_DRAW_ROOT > 0
#include <TApplication.h>
#endif

#if TSP_DEBUG_PRINT > 0
#include <thread>
#endif

#if TSP_ERROR > 1
#include <stdexcept>
#endif

#include "utils.hpp"

using TimePoint_t = std::chrono::steady_clock::time_point;
using TimeMicroS_t = std::chrono::microseconds;
using TimeMilliS_t = std::chrono::milliseconds;
using TimeS_t = std::chrono::seconds;
#if TSP_TIME > 1
using TimeUnit_t = TimeMicroS_t;
inline const std::string& time_unit{ "us" };
#elif TSP_TIME > 0
using TimeUnit_t = TimeMilliS_t;
inline const std::string& time_unit{ "ms" };
#else
using TimeUnit_t = TimeS_t;
inline const std::string& time_unit{ "s" };
#endif

inline const std::string& pos_sep{ "sep" };

typedef float Value_t;
typedef std::size_t Index_t;
constexpr Value_t VALUE_ZERO{ Value_t{ 0 } };
constexpr Value_t VALUE_ONE{ Value_t{ 1 } };
constexpr Value_t VALUE_ONE_NEG{ Value_t{ -1 } };

constexpr Index_t Num_Nodes_Initial{ 3 };
constexpr bool Validation_Intersection{ true };
constexpr bool Intersection_Recursive{ true };

#if TSP_DEBUG_PRINT > 0
inline bool global_print{ false };
#endif

typedef std::vector<Index_t> Indices_t;
typedef std::optional<Index_t> IndexOpt_t;
template <typename T> using IndexExp_t = utils::Expected<Index_t, T>;

struct City
{
    bool on_stack{false};
    int id, layer{-1};
    int id_prev{-1};
    int id_next{-1};
    int x_cell{-1}, y_cell{-1};
    Value_t x, y, cost{VALUE_ZERO};

    City() = default;

    City(int id, Value_t x, Value_t y)
        : id{ id }
        , x{ x }
        , y{ y }
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
typedef std::unordered_map<int, Indices_t> CityLayers_t;
typedef CityLayers_t Stack_t;
// typedef Indices_t Stack_t;
// typedef std::set<Index_t> Stack_t;
// typedef std::set<Index_t, std::greater<Index_t>> Stack_t;
typedef std::vector<std::vector<Indices_t>> Grid_t;
typedef std::pair<Index_t, Index_t> Edge_t;
typedef std::set<Edge_t> Edges_t;
typedef std::vector<std::vector<Edges_t>> GridEdge_t;
typedef std::pair<int, Value_t> IndexCost_t;
typedef std::array<IndexCost_t, 8> Faces_t;
typedef std::vector<Faces_t> FacesVect_t;
typedef std::pair<Edge_t, Value_t> EdgeScore_t;
typedef std::vector<EdgeScore_t> EdgeScores_t;
typedef std::vector<EdgeScores_t> EdgeScoresVect_t;

template <> inline bool utils::MatchItem<City>::operator()(const City& city)
{
    return (city.id == m_item.id);
}

template <> inline bool utils::MatchItem<Node_t>::operator()(const Node_t& node)
{
    return (node->id == m_item->id);
}

class DiscreteENN_TSP;
#if TSP_DRAW_ROOT > 0
template<typename ListType = Indices_t>
void drawState(const std::string& state_file, const DiscreteENN_TSP& enn_tsp, const ListType& highlight_list = ListType{}, std::size_t timeout=-1, bool wait=false);
#endif
void drawPath(const DiscreteENN_TSP& enn_tsp, bool draw_coords,
              const std::string& title, float close_time = 3600.0,
              int highlight = -1, int pos_start = -1, int pos_end = -1);
void savePath(const DiscreteENN_TSP& enn_tsp, bool draw_coords,
              const std::string& title, int highlight = -1, int pos_start = -1, int pos_end = -1);

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

    auto minmaxConst() const
    {
        const Value_t min{ std::min(min_x, min_y) };
        const Value_t max{ std::max(max_x, max_y) };
        return std::make_pair(min, max);
    }

    auto value()
    {
        return std::tie(min_x, max_x, min_y, max_y);
    }

    auto valueConst() const
    {
        return std::make_tuple(min_x, max_x, min_y, max_y);
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

inline Value_t getDistance(const City& a, const City& b)
{
    const auto x_sqr = (a.x - b.x) * (a.x - b.x);
    const auto y_sqr = (a.y - b.y) * (a.y - b.y);
    return std::sqrt(x_sqr + y_sqr);
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

// 0 --> collinear
// 1 --> Clockwise
// 2 --> Counterclockwise
inline int getOrientation(const City& a, const City& b, const City& c)
{
    const Value_t val = (b.y - a.y) * (c.x - b.x) - (b.x - a.x) * (c.y - b.y);

    if (utils::isEqual(val, VALUE_ZERO))
        return 0;

    return (val > 0) ? 1 : 2;
}

inline bool sameCoords(const City& cityA, const City& cityB)
{
    return std::make_pair(cityA.x, cityA.y) == std::make_pair(cityB.x, cityB.y);
}

bool isInside(const City& city, const City& cityA, const City& cityB,
              const City& cityC);

bool hasIntersection(const City& cityA1, const City& cityA2, const City& cityB1,
                     const City& cityB2);

Value_t insertionCost(const City& new_city, const City& cityA,
                      const City& cityB);

class DiscreteENN_TSP
{
public:
    DiscreteENN_TSP() = default;

    Cities_t& cities()
    {
        return m_cities;
    }
    Stack_t& stack()
    {
        return m_stack;
    }
    Indices_t& path()
    {
        return m_path;
    }
    const Cities_t& cities() const
    {
        return m_cities;
    }
    const Stack_t& stack() const
    {
        return m_stack;
    }
    const Indices_t& path() const
    {
        return m_path;
    }
    Grid_t& grid()
    {
        return m_grid;
    }
    GridEdge_t& gridEdges()
    {
        return m_gridEdges;
    }
    EdgeScores_t& edgeScores(std::size_t index)
    {
        return m_edgeScores[index];
    }
    EdgeScoresVect_t& edgeScoresAll()
    {
        return m_edgeScores;
    }
    std::string& name()
    {
        return m_name;
    }
    const std::string& name() const
    {
        return m_name;
    }
    int& initialSize()
    {
        return m_initialSize;
    }
    int& layers()
    {
        return m_layers;
    }
    int& gridSize()
    {
        return m_gridSize;
    }
    auto gridSteps()
    {
        return std::tie(m_gridStepX, m_gridStepY);
    }
    auto gridSteps() const
    {
        return std::tuple(m_gridStepX, m_gridStepY);
    }
    auto gridStepMinMax()
    {
        return std::tie(m_gridStepMin, m_gridStepMax);
    }
    auto gridStepMinMax() const
    {
        return std::tuple(m_gridStepMin, m_gridStepMax);
    }
    auto gridStart()
    {
        return std::tie(m_gridStartX, m_gridStartY);
    }
    bool& rmIntersectRecurse()
    {
        return m_rmIntersectRecurse;
    }
    bool& validIntersectCK()
    {
        return m_validIntersectCK;
    }
    bool& shuffleCities()
    {
        return m_shuffleCities;
    }
    bool& draw()
    {
        return m_draw;
    }
    bool& drawFailed()
    {
        return m_drawFailed;
    }
    bool& drawCoords()
    {
        return m_drawCoords;
    }

    Index_t properIndex(int index) const
    {
        const int n = m_path.size();
        index = (index < 0) ? (n + index) : index;
        index = (index > (n - 1)) ? (index - n) : index;
        return static_cast<size_t>(index);
    }

    auto getNeigbhours(int index) const
    {
        return std::make_pair(properIndex(index - 1), properIndex(index + 1));
    }

    Index_t stackPopBack()
    {
        const Index_t result{ m_stack[m_layer].back() };
        m_stack[m_layer].pop_back();
        // const Stack_t::const_iterator back{ --(m_stack.end()) };
        // const Index_t result{ *back };
        // m_stack.erase(back);
        m_cities[result].on_stack = false;
        return result;
    }

    Index_t stackPopAt(Index_t index)
    {
        Indices_t& stack_slice{ m_stack[m_layer] };
        const Indices_t::const_iterator it{ stack_slice.begin() + index };
        const Index_t result{ *it };
        stack_slice.erase(it);
        // const Stack_t::const_iterator it{ std::next(m_stack.begin(), index) };
        // const Index_t result{ *it };
        // m_stack.erase(it);
        m_cities[result].on_stack = false;
        return result;
    }

    void removeNode(Index_t index)
    {
        auto tmp = m_pattern;
        const Index_t pos{ m_path[index] };
        const std::string& node_idx_str{ pos_sep + std::to_string(pos) +
                                         pos_sep };
        std::string::size_type idx_erase = m_pattern.find(node_idx_str);
        m_pattern.erase(idx_erase, node_idx_str.length());

        const Indices_t::iterator node_iter{ m_path.begin() + index };
        m_path.erase(node_iter);
        m_cities[pos].on_stack = true;
        // m_stack.push_back(pos);
        // m_stack.insert(m_stack.begin(), pos);
        // m_stack.insert(pos);
        m_stack[m_cities[pos].layer].push_back(pos);

        const Index_t pos_next = m_path[properIndex(index)];
        const Index_t pos_prev = m_path[properIndex(int(index) - 1)];
        m_cities[pos_next].id_prev = pos_prev;
        m_cities[pos_prev].id_next = pos_next;
    }

    void addNode(Index_t index, Index_t pos)
    {
        const std::string& node_idx_str{
            pos_sep + std::to_string(m_path[index]) + pos_sep
        };
        const std::string& node_idx_added_str{ pos_sep + std::to_string(pos) +
                                               pos_sep };
        // utils::printInfoFmt("node_idx_str: %s and node_idx_added_str: %s", "addNode", node_idx_str.c_str(), node_idx_added_str.c_str());
        std::string::size_type idx_insert = m_pattern.find(node_idx_str);
        m_pattern.insert(idx_insert, node_idx_added_str);

        if (index == m_path.size()) {
            m_path.push_back(pos);
        } else {
            m_path.insert(m_path.begin() + index, pos);
        }

        const Index_t pos_next = m_path[properIndex(index + 1)];
        const Index_t pos_prev = m_path[properIndex(int(index) - 1)];
        m_cities[pos].id_next = pos_next;
        m_cities[pos].id_prev = pos_prev;
        m_cities[pos_next].id_prev = pos;
        m_cities[pos_prev].id_next = pos;

        m_cities[pos].on_stack = false;
    }

    std::pair<bool, bool> validateNode(Index_t index) const
    {
        const Index_t num_nodes = m_path.size();
        if (num_nodes == 3) {
            return std::make_pair(true, false);
        }
#if TSP_ERROR > 0
        if (num_nodes < 3) {
            utils::printErr("given path size less than 3", "validateNode");
            std::cout.flush();
            std::cerr.flush();
#if TSP_ERROR > 1
            throw std::runtime_error{ "Invalid path size" };
#endif
            return std::make_pair(false, true);
        }
#endif
        const Index_t position{ m_path[index] };
        const City& city{ m_cities[position] };
        const Value_t cost{ city.cost };
        for (Index_t idx{ 0 }; idx != num_nodes; ++idx) {
            const Index_t idx_next{ properIndex(idx + 1) };
            if ((idx == index) or (idx_next == index)) {
                continue;
            }
            const Index_t pos{ m_path[idx] };
            const Index_t pos_next{ m_path[idx_next] };
            const Value_t cost_current{ insertionCost(city, m_cities[pos],
                                                      m_cities[pos_next]) };
            if (cost_current < cost) {
#if TSP_DEBUG_PRINT > 0
                if (global_print) {
                    utils::printInfo("cost_current " +
                                         std::to_string(cost_current) +
                                         " and cost " + std::to_string(cost),
                                     "validateNode");
                    utils::printInfo(
                        "Removing " + std::to_string(position) +
                            " from index " + std::to_string(index) +
                            " due to edge (" + std::to_string(idx) + ", " +
                            std::to_string(idx_next) + ")",
                        "validateNode");
                }
#endif
                return std::make_pair(false, false);
            }
        }
        return std::make_pair(true, false);
    }

    Value_t updateCost(Index_t index)
    {
#if TSP_ERROR > 0
        const int num_nodes = m_path.size();
        if (num_nodes < 3) {
            utils::printErr("given path size (" + std::to_string(num_nodes) +
                                ") less than 3",
                            "updateCost");
            utils::printErr("cities size (" + std::to_string(m_cities.size()) +
                                ")",
                            "updateCost");
            utils::printErr("index passed (" + std::to_string(index) + ")",
                            "updateCost");
            std::cout.flush();
            std::cerr.flush();
#if TSP_ERROR > 1
            throw std::runtime_error{ "Invalid path size" };
#endif
            return VALUE_ONE_NEG;
        }
#endif
        const auto [idx_prev, idx_next] = getNeigbhours(index);
        const Index_t pos_prev{ m_path[idx_prev] }, pos_curr{ m_path[index] },
            pos_next{ m_path[idx_next] };
        const Value_t cost = insertionCost(
            m_cities[pos_curr], m_cities[pos_prev], m_cities[pos_next]);
        m_cities[pos_curr].cost = cost;
        return cost;
    }

    std::pair<int, bool> updateCostNeighbour(Index_t index, bool self)
    {
#if TSP_ERROR > 0
        const Index_t num_nodes = m_path.size();
        if (num_nodes < 3) {
            utils::printErr("given path size less than 3",
                            "updateCostNeighbour");
            std::cout.flush();
            std::cerr.flush();
#if TSP_ERROR > 1
            throw std::runtime_error{ "Invalid path size" };
#endif
            return std::make_pair(-1, true);
        }
#endif
        const auto [idx_prev, idx_next] = getNeigbhours(index);
        if (self and updateCost(index) < VALUE_ZERO) {
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
        if (updateCost(idx_prev) < VALUE_ZERO) {
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
        if (updateCost(idx_next) < VALUE_ZERO) {
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
        const Index_t num_nodes = m_path.size();
        for (Index_t idx{ 0 }; idx != num_nodes; ++idx) {
            if (updateCost(idx) < VALUE_ZERO) {
                const auto [idx_prev, idx_next] = getNeigbhours(idx);
                utils::printErr("Invalid node cost calculate for node at " +
                                    std::to_string(idx) +
                                    " current path size " +
                                    std::to_string(m_path.size()),
                                "updateCostAll");
                utils::printErr("Previous node " + std::to_string(idx_prev),
                                "updateCostAll");
                m_cities[m_path[idx_prev]].print();
                utils::printErr("Current node " + std::to_string(idx),
                                "updateCostAll");
                m_cities[m_path[idx]].print();
                utils::printErr("Next node " + std::to_string(idx_next),
                                "updateCostAll");
                m_cities[m_path[idx_next]].print();
                return std::make_pair(idx, true);
            }
        }
        return std::make_pair(num_nodes, false);
    }

    IndexExp_t<bool> validatePath()
    {
        IndexOpt_t pos_erased{ std::nullopt };
        Index_t num_nodes = m_path.size();
#if TSP_ERROR > 0
        if (num_nodes < 3) {
            utils::printErr("given path size less than 3",
                            "updateCostNeighbour");
            std::cout.flush();
            std::cerr.flush();
#if TSP_ERROR > 1
            throw std::runtime_error{ "Invalid path size" };
#endif
            return IndexExp_t<bool>{ pos_erased, true };
        }
#endif
        for (Index_t idx{ 0 }; idx != num_nodes;) {
            if (m_fromScratch or (num_nodes == 3)) {
                break;
            }
            const auto [valid, err] = validateNode(idx);
            if (err) {
                utils::printErr(
                    "first validateNode failed at " + std::to_string(idx) +
                        " current path size " + std::to_string(m_path.size()),
                    "validatePath");
                return IndexExp_t<bool>{ pos_erased, true };
            }
            if (valid) {
                ++idx;
                continue;
            }
            const Index_t pos = m_path[idx];
            removeNode(idx);
            --num_nodes;
            if (pos_erased.has_value()) {
                pos_erased = std::min(*pos_erased, pos);
            } else {
                pos_erased = pos;
            }
            const Index_t end = properIndex(idx);
            const Index_t start = properIndex(properIndex(int(end) - 1));
            updateCost(start);
            updateCost(end);
            if (num_nodes == 3) {
                break;
            }
            const auto erased = validateEdge(start, end);
            if (erased.err()) {
                utils::printErr(
                    "validateEdge failed with (" + std::to_string(start) +
                        ", " + std::to_string(end) + ") current path size " +
                        std::to_string(m_path.size()),
                    "validatePath");
                return IndexExp_t<bool>{ pos_erased, true };
            }
            if (erased.has_value()) {
                pos_erased = pos_erased.has_value() ?
                                 std::min(*pos_erased, erased.value()) :
                                 erased.opt();
                num_nodes = m_path.size();
            }
            idx = 0;
        }
        return IndexExp_t<bool>{ pos_erased, false };
    }

    int findBestInsertion(Index_t position) const
    {
#if TSP_ERROR > 0
        if (m_path.size() < 2) {
            utils::printErr("given path size less than 2", "findBestInsertion");
#if TSP_ERROR > 1
            throw std::runtime_error{ "Invalid path size" };
#endif
            return -1;
        }
#endif
        const City& new_city{ m_cities[position] };
        const int gridX = new_city.x_cell;
        const int gridY = new_city.y_cell;

        int search_len{0}, best_position{-1};
        Index_t num_nodes{0};
        Value_t min_cost{ std::numeric_limits<Value_t>::max() };
        const Indices_t& cell_positions = m_grid[gridX][gridY];
        utils::printInfoFmt("cell size %u, current position %i, (gridX, gridY) (%i, %i), best_position %i, num_nodes %u, grids %i", "findBestInsertion", cell_positions.size(), position, gridX, gridY, best_position, num_nodes, m_gridSize);
        for (const Index_t& cell_position : cell_positions) {
            // utils::printInfoFmt("cell position %u, current position %i, (gridX, gridY) (%i, %i), best_position %i, num_nodes %u", "findBestInsertion", cell_position, position, gridX, gridY, best_position, num_nodes);
            if (cell_position == position) {
                continue;
            }
            const City& city{ m_cities[cell_position] };
            if (city.on_stack) {
                continue;
            }
            const Index_t pos_prev = city.id_prev;
            const Index_t pos_next = city.id_next;
            // utils::printInfoFmt("cell position %u, (pos_prev, pos_next) (%i, %i)", "findBestInsertion", cell_position, pos_prev, pos_next);
            Value_t cost{ insertionCost(new_city, city, m_cities[pos_next]) };
            if (cost < min_cost) {
                min_cost = cost;
                best_position = pos_next;
            }
            // utils::printInfoFmt("cell position %u, current position %i, (gridX, gridY) (%i, %i), best_position %i, num_nodes %u, cost %f", "findBestInsertion", cell_position, position, gridX, gridY, best_position, num_nodes, cost);
            cost = insertionCost(new_city, m_cities[pos_prev], city);
            if (cost < min_cost) {
                min_cost = cost;
                best_position = cell_position;
            }
            ++num_nodes;
            // utils::printInfoFmt("cell position %u, current position %i, (gridX, gridY) (%i, %i), best_position %i, num_nodes %u, cost %f", "findBestInsertion", cell_position, position, gridX, gridY, best_position, num_nodes, cost);
        }
        while (num_nodes == 0) {
            ++search_len;
            for (int x_offset{-search_len}; x_offset != (search_len + 1); ++x_offset) {
                for (int y_offset{-search_len}; y_offset != (search_len + 1); ++y_offset) {
                    if ( std::abs(x_offset) != search_len and std::abs(y_offset) != search_len) {
                        continue;
                    }
                    const int x_cell = gridX + x_offset;
                    const int y_cell = gridY + y_offset;
                    if ((x_cell < 0) or (x_cell >= m_gridSize)) {
                        continue;
                    }
                    if ((y_cell < 0) or (y_cell >= m_gridSize)) {
                        continue;
                    }
                    const Indices_t& cell_positions = m_grid[x_cell][y_cell];
                    // utils::printInfoFmt("cell size %u, current position %i, (x_cell, y_cell) (%i, %i), best_position %i, num_nodes %u", "findBestInsertion", cell_positions.size(), position, x_cell, y_cell, best_position, num_nodes);
                    for (const Index_t& cell_position : cell_positions) {
                        // utils::printInfoFmt("cell position %u, current position %i, (x_cell, y_cell) (%i, %i), best_position %i, num_nodes %u", "findBestInsertion", cell_position, position, x_cell, y_cell, best_position, num_nodes);
                        if (cell_position == position) {
                            continue;
                        }
                        const City& city{ m_cities[cell_position] };
                        if (city.on_stack) {
                            continue;
                        }
                        const int pos_prev = city.id_prev;
                        const int pos_next = city.id_next;
                        // utils::printInfoFmt("cell position %u, (pos_prev, pos_next) (%i, %i)", "findBestInsertion", cell_position, pos_prev, pos_next);
                        Value_t cost{ insertionCost(new_city, city, m_cities[pos_next]) };
                        if (cost < min_cost) {
                            min_cost = cost;
                            best_position = pos_next;
                        }
                        // utils::printInfoFmt("cell position %u, current position %i, (x_cell, y_cell) (%i, %i), best_position %i, num_nodes %u, cost %f", "findBestInsertion", cell_position, position, x_cell, y_cell, best_position, num_nodes, cost);
                        cost = insertionCost(new_city, m_cities[pos_prev], city);
                        if (cost < min_cost) {
                            min_cost = cost;
                            best_position = cell_position;
                        }
                        ++num_nodes;
                        // utils::printInfoFmt("cell position %u, current position %i, (x_cell, y_cell) (%i, %i), best_position %i, num_nodes %u, cost %f", "findBestInsertion", cell_position, position, x_cell, y_cell, best_position, num_nodes, cost);
                    }
                }
            }
        }

        for (Index_t idx{0}; idx != m_path.size(); ++idx) {
            if (static_cast<int>(m_path[idx]) == best_position) {
                return idx;
            }
        }

        return -1;
    }

    auto findEdge(Index_t pos_start, Index_t pos_end)
    {
        const Index_t num_nodes = m_path.size();
        if (pos_start == m_path.back() and pos_end == m_path[0]) {
            return std::make_pair(int(num_nodes - 1), 0);
        }
        for (Index_t idx{ 0 }; idx != (num_nodes - 1); ++idx) {
            if (pos_start == m_path[idx] and pos_end == m_path[idx + 1]) {
                return std::make_pair(int(idx), int(idx + 1));
            }
        }
        return std::make_pair(-1, -1);
    }

    IndexExp_t<bool> removeIntersectionEdge(Index_t start, Index_t end)
    {
        IndexOpt_t pos_erased{ std::nullopt };
        const Index_t pos_start{ m_path[start] };
        const Index_t pos_end{ m_path[end] };
        const City& city_start{ m_cities[pos_start] };
        const City& city_end{ m_cities[pos_end] };
        Index_t num_nodes = m_path.size();
        for (Index_t idx{ 0 }; idx != num_nodes;) {
            if (m_fromScratch or city_start.on_stack or city_end.on_stack) {
                break;
            }
            const Index_t idx_next{ properIndex(idx + 1) };
            const Index_t pos{ m_path[idx] };
            const Index_t pos_next{ m_path[idx_next] };
            const bool ok_start = pos_start != pos and pos_start != pos_next;
            const bool ok_end = pos_end != pos and pos_end != pos_next;
            if ((not ok_start) or (not ok_end)) {
                ++idx;
                continue;
            }
            const City& city{ m_cities[pos] };
            const City& city_next{ m_cities[pos_next] };
            if (not hasIntersection(city_start, city_end, city, city_next)) {
                ++idx;
                continue;
            }
#if TSP_DEBUG_PRINT > 0
            if (global_print) {
                const auto [start_tmp, end_tmp] = findEdge(pos_start, pos_end);
                utils::printInfoFmt("hasIntersection true for the edges (%u, %u) and (%u, %u)", "removeIntersectionEdge", start_tmp, end_tmp, idx, idx_next);
                utils::printInfoFmt("with the associated cities (%i, %i) and (%i, %i)", "removeIntersectionEdge", pos_start, pos_end, pos, pos_next);
                city_start.print();
                city_end.print();
                city.print();
                city_next.print();
            }
#endif
            Index_t idx_remove{ idx };
            Index_t idx_remove_next{ idx_next };
            // const Value_t dist1{ getDistance(city_start, city_end) };
            // const Value_t dist2{ getDistance(city, city_next) };
            // if ((not utils::isEqual(dist1, dist2)) and (dist1 > dist2)) {
            //     auto [start_new, end_new] = findEdge(pos_start, pos_end);
            //     if (start_new == -1) {
            //         utils::printErr("Tried to find a non-existing edge");
            //         throw std::runtime_error{
            //             "Tried to find a non-existing edge"
            //         };
            //     }
            //     idx_remove = start_new;
            //     idx_remove_next = end_new;
            // }
            if (pos_erased.has_value()) {
                pos_erased =
                    std::min(*pos_erased, std::min(m_path[idx_remove],
                                                   m_path[idx_remove_next]));
            } else {
                pos_erased =
                    std::min(m_path[idx_remove], m_path[idx_remove_next]);
            }
            removeNode(idx_remove);
            removeNode(properIndex(idx_remove));
            num_nodes -= 2;
            if (num_nodes == 2) {
                m_fromScratch = true;
                break;
            }
            end = idx_remove_next == 0 ? 0 : properIndex(idx_remove);
            start = properIndex(int(end) - 1);
            updateCost(start);
            updateCost(end);
            const auto erased = validateEdge(start, end);
            if (erased.err()) {
                utils::printErr(
                    "validateEdge failed with (" + std::to_string(start) +
                        ", " + std::to_string(end) + ") current path size " +
                        std::to_string(m_path.size()),
                    "removeIntersectionEdge");
                return IndexExp_t<bool>{ pos_erased, true };
            }
            if (erased.has_value()) {
                pos_erased = pos_erased.has_value() ?
                                 std::min(*pos_erased, erased.value()) :
                                 erased.opt();
                num_nodes = m_path.size();
            }
            idx = 0;
        }
        return IndexExp_t<bool>{ pos_erased, false };
    }

    IndexExp_t<bool> validateEdgeNodes(Index_t start, Index_t end)
    {
        IndexOpt_t pos_erased{ std::nullopt };
        const auto [valid1, err1] = validateNode(start);
        if (err1) {
            utils::printErr("first validateNode failed at " +
                                std::to_string(start) + " current path size " +
                                std::to_string(m_path.size()),
                            "validateEdgeNodes");
            return IndexExp_t<bool>{ pos_erased, true };
        }
        if (valid1) {
            const auto [valid2, err2] = validateNode(end);
            if (err2) {
                utils::printErr(
                    "second validateNode failed at " + std::to_string(end) +
                        " current path size " + std::to_string(m_path.size()),
                    "validateEdgeNodes");
                return IndexExp_t<bool>{ pos_erased, true };
            }
            if (not valid2) {
                const Index_t pos = m_path[end];
                removeNode(end);
                if (pos_erased.has_value()) {
                    pos_erased = std::min(*pos_erased, pos);
                } else {
                    pos_erased = pos;
                }
                end = properIndex(end);
                start = properIndex(properIndex(int(end) - 1));
                updateCost(start);
                updateCost(end);
                if (m_path.size() == 3) {
                    return IndexExp_t<bool>{ pos_erased, false };
                }
                const auto erased1 = validateEdge(start, end);
                if (erased1.err()) {
                    utils::printErr("validateEdge failed with (" +
                                        std::to_string(start) + ", " +
                                        std::to_string(end) +
                                        ") current path size " +
                                        std::to_string(m_path.size()),
                                    "validateEdgeNodes");
                    return IndexExp_t<bool>{ pos_erased, true };
                }
                if (erased1.has_value()) {
                    pos_erased = pos_erased.has_value() ?
                                     std::min(*pos_erased, erased1.value()) :
                                     erased1.opt();
                }
            }
        } else {
            const Index_t pos = m_path[start];
            removeNode(start);
            if (pos_erased.has_value()) {
                pos_erased = std::min(*pos_erased, pos);
            } else {
                pos_erased = pos;
            }
            end = properIndex(start);
            start = properIndex(properIndex(int(end) - 1));
            updateCost(start);
            updateCost(end);
            if (m_path.size() == 3) {
                return IndexExp_t<bool>{ pos_erased, false };
            }
            const auto erased2 = validateEdge(start, end);
            if (erased2.err()) {
                utils::printErr(
                    "validateEdge failed with (" + std::to_string(start) +
                        ", " + std::to_string(end) + ") current path size " +
                        std::to_string(m_path.size()),
                    "validateEdgeNodes");
                return IndexExp_t<bool>{ pos_erased, true };
            }
            if (erased2.has_value()) {
                pos_erased = pos_erased.has_value() ?
                                 std::min(*pos_erased, erased2.value()) :
                                 erased2.opt();
            }
        }
        return IndexExp_t<bool>{ pos_erased, false };
    }
    IndexExp_t<bool> validateWithEdge(Index_t start, Index_t end)
    {
        IndexOpt_t pos_erased{ std::nullopt };
        const Index_t pos_start{ m_path[start] };
        const Index_t pos_end{ m_path[end] };
        const City& city_start{ m_cities[pos_start] };
        const City& city_end{ m_cities[pos_end] };
        for (Index_t idx{ 0 }; idx != m_path.size();) {
            if (m_fromScratch or (m_path.size() == 3) or city_start.on_stack or city_end.on_stack) {
                break;
            }
            const Index_t pos{ m_path[idx] };
            if (pos == pos_start or pos == pos_end) {
                ++idx;
                continue;
            }
            const City& city{ m_cities[pos] };
            const Value_t cost{ insertionCost(city, city_start, city_end) };
            if (cost < city.cost) {
#if TSP_DEBUG_PRINT > 0
                if (global_print) {
                    const auto [start_tmp, end_tmp] = findEdge(pos_start, pos_end);
                    utils::printInfoFmt("validation failed for the node %u with edge (%u, %u)", "validateWithEdge", idx, start_tmp, end_tmp);
                    utils::printInfoFmt("with the associated cities %i and (%i, %i)", "validateWithEdge", pos, pos_start, pos_end);
                    utils::printInfoFmt("with the associated cost %f and city.cost %f", "validateWithEdge", cost, city.cost);
                    city_start.print();
                    city_end.print();
                    city.print();
                }
#endif
                removeNode(idx);
                if (pos_erased.has_value()) {
                    pos_erased = std::min(*pos_erased, pos);
                } else {
                    pos_erased = pos;
                }
                end = properIndex(idx);
                start = properIndex(properIndex(int(end) - 1));
                updateCost(start);
                updateCost(end);
                if (m_path.size() == 3) {
                    return IndexExp_t<bool>{ pos_erased, false };
                }
                const auto erased = validateEdge(start, end);
                if (erased.err()) {
                    utils::printErr("validateEdge failed with (" +
                                        std::to_string(start) + ", " +
                                        std::to_string(end) +
                                        ") current path size " +
                                        std::to_string(m_path.size()),
                                    "validateWithEdge");
                    return IndexExp_t<bool>{ pos_erased, true };
                }
                if (erased.has_value()) {
                    pos_erased = pos_erased.has_value() ?
                                     std::min(*pos_erased, erased.value()) :
                                     erased.opt();
                }
                idx = 0;
            } else {
                ++idx;
            }
        }
        return IndexExp_t<bool>{ pos_erased, false };
    }
    IndexExp_t<bool> validateEdge(Index_t start, Index_t end)
    {
        IndexOpt_t pos_erased{ std::nullopt };
        if (m_path.size() == 3) {
            return IndexExp_t<bool>{ pos_erased, false };
        }
        const Index_t pos_start{ m_path[start] };
        const Index_t pos_end{ m_path[end] };
        const auto erased1 = removeIntersectionEdge(start, end);
        if (erased1.err()) {
            utils::printErr("removeIntersectionEdge failed with (" +
                                std::to_string(start) + ", " +
                                std::to_string(end) + ") current path size " +
                                std::to_string(m_path.size()),
                            "validateEdge");
            return IndexExp_t<bool>{ pos_erased, true };
        }
        if (erased1.has_value()) {
            pos_erased = pos_erased.has_value() ?
                             std::min(*pos_erased, erased1.value()) :
                             erased1.opt();
        }
        if ((not m_fromScratch) and (not m_cities[pos_start].on_stack) and
            (not m_cities[pos_end].on_stack)) {
            auto [start_new, end_new] = findEdge(pos_start, pos_end);
            const auto erased2 = validateEdgeNodes(start_new, end_new);
            if (erased2.err()) {
                utils::printErr("validateEdgeNodes failed with (" +
                                    std::to_string(start_new) + ", " +
                                    std::to_string(end_new) +
                                    ") current path size " +
                                    std::to_string(m_path.size()),
                                "validateEdge");
                return IndexExp_t<bool>{ pos_erased, true };
            }
            if (erased2.has_value()) {
                pos_erased = pos_erased.has_value() ?
                                 std::min(*pos_erased, erased2.value()) :
                                 erased2.opt();
            }
        }
        if ((not m_fromScratch) and (not m_cities[pos_start].on_stack) and
            (not m_cities[pos_end].on_stack)) {
            auto [start_new, end_new] = findEdge(pos_start, pos_end);
            const auto erased3 = validateWithEdge(start_new, end_new);
            if (erased3.err()) {
                utils::printErr("validateEdgeNodes failed with (" +
                                    std::to_string(start_new) + ", " +
                                    std::to_string(end_new) +
                                    ") current path size " +
                                    std::to_string(m_path.size()),
                                "validateEdge");
                return IndexExp_t<bool>{ pos_erased, true };
            }
            if (erased3.has_value()) {
                pos_erased = pos_erased.has_value() ?
                                 std::min(*pos_erased, erased3.value()) :
                                 erased3.opt();
            }
        }
        return IndexExp_t<bool>{ pos_erased, false };
    }
    auto toGridCell(Value_t x, Value_t y)
    {
        int x_cell = std::abs(x - m_gridStartX)/m_gridStepX;
        int y_cell = std::abs(y - m_gridStartY)/m_gridStepY;
        x_cell = x_cell == m_gridSize ? m_gridSize - 1 : x_cell;
        y_cell = y_cell == m_gridSize ? m_gridSize - 1 : y_cell;
        return std::make_pair(x_cell, y_cell);
    }
    std::vector<int> validateEdge(Index_t pos_start, Index_t pos_end, const Edges_t& edges)
    {
        const City& city_start{ m_cities[pos_start] };
        const City& city_end{ m_cities[pos_end] };
        const std::size_t num_edges = edges.size();
        std::vector<int> result(num_edges, 0);
        for (int idx{0}; idx < static_cast<int>(num_edges); ++idx) {
            const auto [pos, pos_next] = *std::next(edges.begin(),idx);
            if (pos == pos_start or pos == pos_end) {
                continue;
            }
            if (pos_next == pos_start or pos_next == pos_end) {
                continue;
            }
            if (hasIntersection(city_start, city_end, m_cities[pos], m_cities[pos_next])) {
                result[idx] = 1;
            }
        }
        return result;
    }
    IndexExp_t<bool> markEdgeCells(Index_t pos_start, Index_t pos_end)
    {
        IndexOpt_t pos_erased{ std::nullopt };
        // if (m_path.size() == 3) {
        //     return IndexExp_t<bool>{ pos_erased, false };
        // }
        const City& city_start{ m_cities[pos_start] };
        const City& city_end{ m_cities[pos_end] };
        const Value_t p1_x = city_start.x;
        const Value_t p1_y = city_start.y;
        const Value_t p2_x = city_end.x;
        const Value_t p2_y = city_end.y;
        const Value_t width_x = p2_x - p1_x;
        const Value_t width_y = p2_y - p1_y;
        const Value_t length = std::sqrt(width_x*width_x + width_y*width_y) ;
        Value_t step_startX = p1_x;
        Value_t step_startY = p1_y;
        const int multiplierX = (p1_x < p2_x) ? 1 : -1;
        const int multiplierY = (p1_y < p2_y) ? 1 : -1;
        while ((multiplierX*step_startX < multiplierX*p2_x or utils::isEqual(step_startX, p2_x)) and (multiplierY*step_startY < multiplierY*p2_y or utils::isEqual(step_startY, p2_y))) {
            if (m_fromScratch) {
                break;
            }
            if (city_start.on_stack or city_end.on_stack) {
                unmarkEdgeCells(pos_start, pos_end);
                break;
            }
            const auto [x_cell, y_cell] = toGridCell(step_startX, step_startY);
            if (x_cell >= m_gridSize or y_cell >= m_gridSize) {
                utils::printInfoFmt("Invaid cell indices : (%i, %i)", "markEdgeCells", x_cell, y_cell);
                throw std::runtime_error{utils::stringFmt("Invaid cell indices : (%i, %i)", x_cell, y_cell)};
            }
            Edges_t& edges = m_gridEdges[x_cell][y_cell];
            auto inserted = edges.insert(std::make_pair(pos_start, pos_end));
            if (not inserted.second) {
                step_startX += utils::getRound2((width_x/length)*m_gridStepMin/2.);
                step_startY += utils::getRound2((width_y/length)*m_gridStepMin/2.);
                continue;
            }
            for (auto it{ edges.begin() }; it != edges.end();) {
                if (city_start.on_stack or city_end.on_stack) {
                    unmarkEdgeCells(pos_start, pos_end);
                    break;
                }
                const auto [pos, pos_next] = *it;
                if (pos == pos_start or pos == pos_end) {
                    ++it;
                    continue;
                }
                if (pos_next == pos_start or pos_next == pos_end) {
                    ++it;
                    continue;
                }
                if (hasIntersection(city_start, city_end, m_cities[pos], m_cities[pos_next])) {
                    auto [idx, idx_next] = findEdge(pos, pos_next);
                    if (idx == -1 or idx_next == -1) {
                        // utils::printInfoFmt("Current cities for path size %u", "markEdgeCells", m_path.size());
                        // for (const City& tmp: m_cities) {
                        //     tmp.print();
                        // }
                        utils::printInfoFmt("Found intersection for (%i, %i) and (%i, %i)", "markEdgeCells", pos_start, pos_end, pos, pos_next);
                        utils::printInfoFmt("Current path size %u", "markEdgeCells", m_path.size());
                        for (const Index_t& tmp: m_path) {
                            m_cities[tmp].print();
                        }
                        utils::printInfoFmt("Edges size : %u\nEdges elements", "markEdgeCells", edges.size());
                        for (const auto& edge : edges) {
                            utils::printInfoFmt("Edge : (%i, %i)", "markEdgeCells", edge.first, edge.second);
                        }
                        throw std::runtime_error{utils::stringFmt("[Error] (markEdgeCells): invalid removal indices: %i, %i from positions %i and %i", idx, idx_next, pos, pos_next)};
                    }
                    {
                        const Index_t pos_prev = m_path[properIndex(int(idx) - 1)];
                        unmarkEdgeCells(pos_prev, pos);
                    }
                    {
                        const Index_t pos_next_next = m_path[properIndex(idx_next + 1)];
                        unmarkEdgeCells(pos_next, pos_next_next);
                    }
                    unmarkEdgeCells(pos, pos_next);
                    removeNode(idx);
                    removeNode(properIndex(idx));
                    // edges.erase(it);
                    if (pos_erased.has_value()) {
                        pos_erased =
                            std::min(*pos_erased, std::min(pos,
                                                           pos_next));
                    } else {
                        pos_erased =
                            std::min(pos, pos_next);
                    }
                    if (m_path.size() == 2) {
                        m_fromScratch = true;
                        break;
                    }
                    idx_next = idx_next == 0 ? 0 : properIndex(idx);
                    idx = properIndex(int(idx_next) - 1);
                    updateCost(idx);
                    updateCost(idx_next);
                    auto erased = markEdgeCells(m_path[idx], m_path[idx_next]);
                    if (erased.has_value()) {
                        pos_erased = pos_erased.has_value() ?
                                         std::min(*pos_erased, erased.value()) :
                                         erased.opt();
                    }
                    if (city_start.on_stack or city_end.on_stack) {
                        break;
                    }
                    it = edges.begin();
                } else {
                    ++it;
                }
            }
            step_startX += utils::getRound2((width_x/length)*m_gridStepMin/2.);
            step_startY += utils::getRound2((width_y/length)*m_gridStepMin/2.);
        }

        return IndexExp_t<bool>{ pos_erased, false };
    }
    void unmarkEdgeCells(Index_t pos_start, Index_t pos_end)
    {
        const City& city_start{ m_cities[pos_start] };
        const City& city_end{ m_cities[pos_end] };
        m_gridEdges[city_start.x_cell][city_start.y_cell].erase(std::make_pair(pos_start, pos_end));
        m_gridEdges[city_end.x_cell][city_end.y_cell].erase(std::make_pair(pos_start, pos_end));
        const Value_t p1_x = city_start.x;
        const Value_t p1_y = city_start.y;
        const Value_t p2_x = city_end.x;
        const Value_t p2_y = city_end.y;
        const Value_t width_x = p2_x - p1_x;
        const Value_t width_y = p2_y - p1_y;
        const Value_t length = std::sqrt(width_x*width_x + width_y*width_y) ;
        Value_t step_startX = p1_x;
        Value_t step_startY = p1_y;
        const int multiplierX = (p1_x < p2_x) ? 1 : -1;
        const int multiplierY = (p1_y < p2_y) ? 1 : -1;
        while ((multiplierX*step_startX < multiplierX*p2_x or utils::isEqual(step_startX, p2_x)) and (multiplierY*step_startY < multiplierY*p2_y or utils::isEqual(step_startY, p2_y))) {
            const auto [x_cell, y_cell] = toGridCell(step_startX, step_startY);
            if (x_cell >= m_gridSize or y_cell >= m_gridSize) {
                utils::printInfoFmt("Invaid cell indices : (%i, %i)", "unmarkEdgeCells", x_cell, y_cell);
                throw std::runtime_error{utils::stringFmt("Invaid cell indices : (%i, %i)", x_cell, y_cell)};
            }
            m_gridEdges[x_cell][y_cell].erase(std::make_pair(pos_start, pos_end));
            step_startX += utils::getRound2((width_x/length)*m_gridStepMin/2.);
            step_startY += utils::getRound2((width_y/length)*m_gridStepMin/2.);
        }
        step_startX += utils::getRound2((width_x/length)*m_gridStepMin/2.);
        step_startY += utils::getRound2((width_y/length)*m_gridStepMin/2.);
        const auto [x_cell, y_cell] = toGridCell(step_startX, step_startY);
        if (x_cell >= m_gridSize or y_cell >= m_gridSize) {
            utils::printInfoFmt("Invaid cell indices : (%i, %i)", "unmarkEdgeCells", x_cell, y_cell);
            throw std::runtime_error{utils::stringFmt("Invaid cell indices : (%i, %i)", x_cell, y_cell)};
        }
        m_gridEdges[x_cell][y_cell].erase(std::make_pair(pos_start, pos_end));

    }

    int checkIntersectEdge(Index_t start, Index_t end)
    {
        const Index_t pos_start{ m_path[start] };
        const Index_t pos_end{ m_path[end] };
        const City& city_start{ m_cities[pos_start] };
        const City& city_end{ m_cities[pos_end] };
        const Index_t num_nodes = m_path.size();
        const bool ok_start = start != 0 and start != (num_nodes - 1);
        const bool ok_end = end != 0 and end != (num_nodes - 1);
        if (ok_start and ok_end) {
            const Index_t pos{ m_path[num_nodes - 1] };
            const Index_t pos_next{ m_path[0] };
            const City& city{ m_cities[pos] };
            const City& city_next{ m_cities[pos_next] };
            if (hasIntersection(city_start, city_end, city, city_next)) {
                return (num_nodes - 1);
            }
        }
        for (Index_t idx{ 0 }; idx != (num_nodes - 1); ++idx) {
            const bool ok_start = start != idx and start != (idx + 1);
            const bool ok_end = end != idx and end != (idx + 1);
            if (ok_start and ok_end) {
                const Index_t pos{ m_path[idx] };
                const Index_t pos_next{ m_path[idx + 1] };
                const City& city{ m_cities[pos] };
                const City& city_next{ m_cities[pos_next] };
                if (hasIntersection(city_start, city_end, city, city_next)) {
                    return idx;
                }
            }
        }
        return -1;
    }

    bool checkIntersectPath()
    {
        const Index_t num_nodes = m_path.size();
        int intersect0{ checkIntersectEdge(num_nodes - 1, 0) };
        if (intersect0 != -1) {
            const Index_t intersect_next{ properIndex(intersect0 + 1) };
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
            m_cities[m_path[static_cast<Index_t>(intersect0)]].print();
            utils::printInfo("node " + std::to_string(intersect_next));
            m_cities[m_path[intersect_next]].print();
            return true;
        }
        for (Index_t idx{ 0 }; idx != (num_nodes - 1); ++idx) {
            int intersect{ checkIntersectEdge(idx, idx + 1) };
            if (intersect != -1) {
                const Index_t intersect_next{ properIndex(intersect + 1) };
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
                m_cities[m_path[static_cast<Index_t>(intersect)]].print();
                utils::printInfo("node " + std::to_string(intersect_next));
                m_cities[m_path[intersect_next]].print();
                return true;
            }
        }
        return false;
    }

    void initialize()
    {
        const Index_t num_cities = m_cities.size();
        assert("[Error] (initialize): given number of cities less than 1" &&
               (num_cities > 0));
        if (num_cities < 4) {
            for (Index_t idx{ 0 }; idx != num_cities; ++idx) {
                m_path.push_back(idx);
            }
            return;
        }
        m_path.clear();
        m_path.reserve(num_cities);
        m_pattern = "";
    }

    void constructInitialPath()
    {
        int id_min{-1};
        auto coords_min{ std::pair<Value_t, Value_t>(std::numeric_limits<Value_t>::max(), std::numeric_limits<Value_t>::max()) };
        for (const City& city : m_cities) {
            const auto coords_curr{ std::pair<Value_t, Value_t>{city.x, city.y} };
            if (coords_curr < coords_min) {
                coords_min = coords_curr;
                id_min = city.id;
            }
        }
        // constructFirstConnection(id_min);
    }

    template<typename PatternsType>
    void constructFirstConnection(Index_t id, PatternsType& patterns)
    {
        const EdgeScores_t& edge_scores{ m_edgeScores[id] };
        if (edge_scores.size() == 0) {
            const std::string& error_msg{ utils::stringFmt("constructFirstConnection id %i has edge_scores 0.\nnum_nodes/num_cities %u/%u", id, m_path.size(), m_cities.size()) };
            throw std::runtime_error{error_msg};
        }
        patterns.clear();
        m_path.clear();
        m_pattern = "";
        m_path.push_back(id);
        m_pattern = pos_sep + std::to_string(id) + pos_sep;
        const auto [id_prev, id_next] = change_child ? edge_scores[1].first : edge_scores[0].first;
        m_path.push_back(id_prev);
        m_pattern += pos_sep + std::to_string(id_prev) + pos_sep;
        m_path.push_back(id_next);
        m_pattern += pos_sep + std::to_string(id_next) + pos_sep;
        patterns.insert(m_pattern);

        m_cities[id].on_stack = false;
        m_cities[id_prev].on_stack = false;
        m_cities[id_next].on_stack = false;
        m_indexLeft = 1;
        m_indexRight = 2;

        m_cities[id].id_prev = id_next;
        m_cities[id].id_next = id_prev;

        m_cities[id_prev].id_prev = id;
        m_cities[id_prev].id_next = id_next;

        m_cities[id_next].id_prev = id_prev;
        m_cities[id_next].id_next = id;
    }

    int checkNeighbours()
    {
        for (std::size_t idx{0}; idx != m_path.size(); ++idx) {
            const Index_t id{ m_path[idx] };
            const int id_prev = m_path[properIndex(int(idx) - 1)];
            const int id_next = m_path[properIndex(idx+1)];
            if (m_cities[id].id_prev != id_prev or m_cities[id].id_next != id_next) {
                return idx;
            }
        }
        return -1;
    }

    void updateChildsRemove()
    {
        if (m_indexLeft < m_indexRight) {
            if (m_left) {
                m_indexLeft = properIndex(int(m_indexLeft) - 1);
                m_indexRight = properIndex(int(m_indexRight) - 1);
            } else {
                m_indexRight = properIndex(m_indexRight);
            }
        } else {
            if (m_left) {
                m_indexLeft = properIndex(int(m_indexLeft) - 1);
            } else {
                m_indexLeft = properIndex(int(m_indexLeft) - 1);
                m_indexRight = properIndex(m_indexRight);
            }
        }
    }

    void updateChildsAdd()
    {
        if (m_indexLeft < m_indexRight) {
            if (m_left) {
                m_indexLeft = properIndex(m_indexLeft + 1);
                m_indexRight = properIndex(m_indexRight + 1);
            } else {
                m_indexRight = properIndex(m_indexRight);
            }
        } else {
            if (m_left) {
                m_indexLeft = properIndex(m_indexLeft + 1);
            } else {
                m_indexLeft = properIndex(m_indexLeft + 1);
                m_indexRight = properIndex(m_indexRight);
            }
        }
    }

    void updateParentAndChild(int index_now)
    {

        if (m_left) {
            m_indexChild = properIndex(index_now + 1);
            m_indexParent = properIndex(index_now - 1);
            // id_parent = m_cities[id_now].id_prev;
        } else {
            m_indexChild = properIndex(index_now);
            m_indexParent = properIndex(index_now + 1);
            // id_parent = m_cities[id_now].id_next;
        }
    }

    bool run(std::default_random_engine& gen)
    {
        [[maybe_unused]] const Index_t num_cities = m_cities.size();
        std::unordered_set<std::string> pattern_hashes;
        {

            int id_min{-1};
            auto coords_min{ std::pair<Value_t, Value_t>(std::numeric_limits<Value_t>::max(), std::numeric_limits<Value_t>::max()) };
            for (const City& city : m_cities) {
                const auto coords_curr{ std::pair<Value_t, Value_t>{city.x, city.y} };
                if (coords_curr < coords_min) {
                    coords_min = coords_curr;
                    id_min = city.id;
                }
            }
            constructFirstConnection(id_min, pattern_hashes);
        }
        std::size_t count{0};
        const std::size_t draw_timeout{10};
        while (true) {
            if (m_path.size() == 1) {
                constructFirstConnection(m_path[0], pattern_hashes);
                continue;
            }

            if (m_path.size() == 2) {
                if (m_cities[m_path[0]].id_next != m_path[1] or m_cities[m_path[1]].id_next != m_path[0] or m_cities[m_path[0]].id_prev != m_path[1]  or m_cities[m_path[1]].id_prev != m_path[0] ) {
                    const std::string& error_msg{ utils::stringFmt("(id_prev, id_next) (%u, %u) mismatch with value from m_path[1] %u (%u, %u) for m_path[0] %u.\nnum_nodes/num_cities %u/%u", m_cities[m_path[0]].id_prev, m_cities[m_path[0]].id_next, m_path[1], m_cities[m_path[1]].id_prev, m_cities[m_path[1]].id_next, m_path[0], m_path.size(), m_cities.size()) };
                    utils::printErr(error_msg, "run");
                    drawState("", *this, {}, -1, true);
                    throw std::runtime_error{""};
                }
                // if (child_left != 0 or child_right != 1) {
                //     const std::string& error_msg1{ utils::stringFmt("(child_left, child_right) (%u, %u) mismatch with values (0, 1).\nnum_nodes/num_cities %u/%u", child_left, child_right, m_path.size(), m_cities.size()) };
                //     utils::printErr(error_msg1, "run");
                //     const std::string& error_msg2{ utils::stringFmt("(id_prev, id_next) (%u, %u) mismatch with value from m_path[1] %u (%u, %u) for m_path[0] %u.\nnum_nodes/num_cities %u/%u", m_cities[m_path[0]].id_prev, m_cities[m_path[0]].id_next, m_path[1], m_cities[m_path[1]].id_prev, m_cities[m_path[1]].id_next, m_path[0], m_path.size(), m_cities.size()) };
                //     utils::printErr(error_msg2, "run");
                //     drawState("", *this, {id_now}, -1, true);
                //     throw std::runtime_error{""};
                // }
                // index_child = 1;
            }

            utils::printInfoFmt("Running iteration #%u with left %i", "run", ++count, static_cast<int>(m_left));
            utils::printInfoFmt("m_indexLeft, m_indexRight %u, %u", "run", m_indexLeft, m_indexRight);
            if (int tmp = checkNeighbours(); tmp != -1) {
                utils::printErrFmt("checkNeighbours first Index %i with city %u doesn't have correct (id_prev, id_next) (%u, %u) which should be (%u, %u)", "run", tmp, m_path[tmp], m_cities[m_path[tmp]].id_prev, m_cities[m_path[tmp]].id_next, m_path[properIndex(int(tmp) - 1)], m_path[properIndex(tmp + 1)]);
                throw std::runtime_error{""};
            }
            int& index_now = m_left ? m_indexLeft : m_indexRight;
            updateParentAndChild(index_now);
            Index_t id_now = m_path[index_now];
            Index_t id_parent = m_path[m_indexParent];

            if (id_parent != m_path[m_indexParent]) {
                const std::string& error_msg{ utils::stringFmt("id_parent %u mismatch with value from m_indexParent %i which is %u.\nnum_nodes/num_cities %u/%u", id_parent, m_indexParent, m_path[m_indexParent], m_path.size(), m_cities.size()) };
                utils::printErr(error_msg, "run");
                drawState("", *this, {id_now}, -1, true);
                throw std::runtime_error{""};
            }

            const std::string& info_msg1{ utils::stringFmt("index_now %i, id_now %u, m_indexParent %i, id_parent %u, m_indexChild %i.\nnum_nodes/num_cities %u/%u", index_now, id_now, m_indexParent, id_parent, m_indexChild, m_path.size(), m_cities.size()) };
            utils::printInfo(info_msg1, "run");
            int id_child = -1;
            const EdgeScores_t& edge_scores{ m_edgeScores[id_now] };
            for (int idx{0}; idx != static_cast<int>(edge_scores.size()); ++idx) {
                const auto [id, id_next] = edge_scores[idx].first;
                // const std::string& info_msg{ utils::stringFmt("Searching for edge.\nidx %i, (id, id_next) (%u, %u), index_now %i, id_now %u, index_child %i, id_child %i, index_parent %i, id_parent %u.\nnum_nodes/num_cities %u/%u", idx, id, id_next, index_now, id_now, index_child, id_now, index_parent, id_parent, m_path.size(), m_cities.size()) };
                // utils::printInfo(info_msg, "run");
                if ((id == id_parent and (not m_cities[id_next].on_stack)) or (id_next == id_parent and (not m_cities[id].on_stack))) {
                    continue;
                }
                if (id == id_parent) {
                    id_child = id_next;
                } else if (id_next == id_parent) {
                    id_child = id;
                }
                if (id_child != -1) {
                    const std::string& info_msg1{ utils::stringFmt("Selected edge.\nidx %i, (id, id_next) (%u, %u), index_now %i, id_now %u, m_indexParent %i, id_parent %u, m_indexChild %i, id_child %i.\nnum_nodes/num_cities %u/%u", idx, id, id_next, index_now, id_now, m_indexParent, id_parent, m_indexChild, id_child, m_path.size(), m_cities.size()) };
                    utils::printInfo(info_msg1, "run");
                    const std::string& info_msg2{ utils::stringFmt("Adding node.\nindex_now %i, id_now %u, m_indexParent %i, id_parent %u, m_indexChild %i, id_child %i.\nnum_nodes/num_cities %u/%u", index_now, id_now, m_indexParent, id_parent, m_indexChild, id_child, m_path.size(), m_cities.size()) };
                    utils::printInfo(info_msg2, "run");
                    addNode(m_indexChild, id_child);
                    updateChildsAdd();
                    updateParentAndChild(index_now);
                    Index_t id_now = m_path[index_now];
                    Index_t id_parent = m_path[m_indexParent];
#if TSP_DRAW_ROOT > 0
                    drawState("", *this, {id_now, static_cast<Index_t>(id_child)}, draw_timeout);
#endif
                    const std::string& info_msg3{ utils::stringFmt("Added node.\nindex_now %i, id_now %u, m_indexParent %i, id_parent %u, m_indexChild %i, id_child %i.\nnum_nodes/num_cities %u/%u", index_now, id_now, m_indexParent, id_parent, m_indexChild, id_child, m_path.size(), m_cities.size()) };
                    utils::printInfo(info_msg3, "run");
                    utils::printInfoFmt("m_indexLeft, m_indexRight %u, %u", "run", m_indexLeft, m_indexRight);
                    if (not pattern_hashes.insert(m_pattern).second) {
                        const std::string& info_msg1{ utils::stringFmt("Removing node.\nindex_now %i, id_now %u, m_indexParent %i, id_parent %u, m_indexChild %i, id_child %i.\nnum_nodes/num_cities %u/%u", index_now, id_now, m_indexParent, id_parent, m_indexChild, id_child, m_path.size(), m_cities.size()) };
                        utils::printInfo(info_msg1, "run");
                        id_child = -1;
                        removeNode(m_indexChild);
                        updateChildsRemove();
                        updateParentAndChild(index_now);
                        Index_t id_now = m_path[index_now];
                        Index_t id_parent = m_path[m_indexParent];
                        drawState("", *this, {id_now}, draw_timeout);
                        const std::string& info_msg2{ utils::stringFmt("Removed node.\nindex_now %i, id_now %u, m_indexParent %i, id_parent %u, m_indexChild %i, id_child %i.\nnum_nodes/num_cities %u/%u", index_now, id_now, m_indexParent, id_parent, m_indexChild, id_child, m_path.size(), m_cities.size()) };
                        utils::printInfo(info_msg2, "run");
                        utils::printInfoFmt("m_indexLeft, m_indexRight %u, %u", "run", m_indexLeft, m_indexRight);
                        continue;
                    }
                    break;
                }
            }

            if (int tmp = checkNeighbours(); tmp != -1) {
                utils::printErrFmt("checkNeighbours second Index %i with city %u doesn't have correct (id_prev, id_next) (%u, %u) which should be (%u, %u)", "run", tmp, m_path[tmp], m_cities[m_path[tmp]].id_prev, m_cities[m_path[tmp]].id_next, m_path[properIndex(int(tmp) - 1)], m_path[properIndex(tmp + 1)]);
                throw std::runtime_error{""};
            }

            if (id_child == -1) {
                removeNode(m_indexChild);
                const std::string& info_msg1{ utils::stringFmt("Removing node.\nindex_now %i, id_now %u, m_indexParent %i, id_parent %u, m_indexChild %i, id_child %i.\nnum_nodes/num_cities %u/%u", index_now, id_now, m_indexParent, id_parent, m_indexChild, id_child, m_path.size(), m_cities.size()) };
                utils::printInfo(info_msg1, "run");
                drawState("", *this, {id_now}, draw_timeout);
                updateChildsRemove();
                updateParentAndChild(index_now);
                Index_t id_now = m_path[index_now];
                Index_t id_parent = m_path[m_indexParent];
                const std::string& info_msg2{ utils::stringFmt("Removed node.\nindex_now %i, id_now %u, m_indexParent %i, id_parent %u, m_indexChild %i, id_child %i.\nnum_nodes/num_cities %u/%u", index_now, id_now, m_indexParent, id_parent, m_indexChild, id_child, m_path.size(), m_cities.size()) };
                utils::printInfo(info_msg2, "run");
                utils::printInfoFmt("m_indexLeft, m_indexRight %u, %u", "run", m_indexLeft, m_indexRight);
                continue;
                // const std::string& error_msg{ utils::stringFmt("id_child is -1 for index_now %i, edge_scores size %u, index_parent %i, index_child %i, id_now %u and id_parent %u, change_child %i.\nnum_nodes/num_cities %u/%u", index_now, edge_scores.size(), index_parent, index_child, id_now, id_parent, static_cast<int>(change_child), m_path.size(), m_cities.size()) };
                // throw std::runtime_error{error_msg};
            }

            if (int tmp = checkNeighbours(); tmp != -1) {
                utils::printErrFmt("checkNeighbours second Index %i with city %u doesn't have correct (id_prev, id_next) (%u, %u) which should be (%u, %u)", "run", tmp, m_path[tmp], m_cities[m_path[tmp]].id_prev, m_cities[m_path[tmp]].id_next, m_path[properIndex(int(tmp) - 1)], m_path[properIndex(tmp + 1)]);
                throw std::runtime_error{""};
            }

            m_left = not m_left;

#if TSP_DEBUG_PRINT > 1
            if (checkIntersectPath()) {
                utils::printErrFmt(
                    "intersection from adding node at %i for city %i.\nnum_nodes/num_cities %u/%u", "run", m_indexChild, id_child, m_path.size(), m_cities.size());
                // std::chrono::seconds sleep_time{ 5 };
                // std::this_thread::sleep_for(sleep_time);
                // drawState("", *this);
            }
#endif

            const Index_t num_nodes = m_path.size();
            if (num_nodes == num_cities) {
                break;
            }
        }

        // for (std::size_t idx{ 0 }; idx != num_cities; ++idx) {
        //     const auto [valid, err] = validateNode(idx);
        //     if (err) {
        //         std::cerr
        //             << "[Error] (run): Algoirthm has not found the optimal path\n";
        //         return false;
        //     }
        // }
        // IndexExp_t<bool> erased = validatePath();
        // if (erased.err()) {
        //     std::cerr << "[Error] (run): final validatePath failed\n";
        //     return false;
        // }
        // if (erased.has_value()) {
        //     std::cerr << "[Error] (run): final validatePath removed node(s)\n";
        //     return false;
        // }
        if (checkIntersectPath()) {
            std::cerr << "[Error] (run): final checkIntersectPath failed\n";
            return false;
        }

        return true;
    }

private:
    bool m_rmIntersectRecurse{ Intersection_Recursive };
    bool m_validIntersectCK{ Validation_Intersection };
    bool m_fromScratch{ false };
    bool m_shuffleCities{ false };
    bool m_draw{ false };
    bool m_drawFailed{ false };
    bool m_drawCoords{ false };
    bool change_child{false}, change_parent{false}, m_left;
    int m_initialSize{ Num_Nodes_Initial };
    int m_layer{ -1 };
    int m_layers{ -1 };
    int m_gridSize{ -1 };
    int m_indexLeft{-1}, m_indexRight{-1};
    int m_indexParent{-1}, m_indexChild{-1};
    Value_t m_gridStepX{ -1 };
    Value_t m_gridStepY{ -1 };
    Value_t m_gridStepMin{ -1 };
    Value_t m_gridStepMax{ -1 };
    Value_t m_gridStartX{ -1 };
    Value_t m_gridStartY{ -1 };
    Indices_t* m_stackSlice;
    std::string m_name;
    std::string m_pattern;
    Cities_t m_cities;
    Stack_t m_stack;
    Indices_t m_path;
    Grid_t m_grid;
    GridEdge_t m_gridEdges;
    EdgeScoresVect_t m_edgeScores;
};

void parseCities(Cities_t& cities, const std::string& filename);

int createStack(DiscreteENN_TSP& enn_tsp, const MinMaxCoords& minmax_coords);

void makeGridInfo(DiscreteENN_TSP& enn_tsp, int depth, const MinMaxCoords& minmax_coords);

void makeConnections(DiscreteENN_TSP& enn_tsp);

struct TSPInfo
{
    int m_points{ -1 };
    Value_t m_distance{ VALUE_ONE_NEG };
    Value_t m_error{ VALUE_ONE_NEG };
    Value_t m_time{ VALUE_ONE_NEG };
    Value_t m_timePerIter{ VALUE_ONE_NEG };
    Value_t m_timePerCity{ VALUE_ONE_NEG };
    Value_t m_timePerIterMin{ VALUE_ONE_NEG };
    Value_t m_timePerIterMax{ VALUE_ONE_NEG };
    std::string m_name{ "" };
};
template <>
inline bool utils::MatchItem<TSPInfo>::operator()(const TSPInfo& info)
{
    return (info.m_name == m_item.m_name);
}
typedef std::vector<TSPInfo> TSPInfoVect_t;

bool readDistances(const std::string& filename, TSPInfoVect_t& infos);
