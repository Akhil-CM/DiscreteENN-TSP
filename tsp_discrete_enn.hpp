// -*- C++ Header -*-
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <optional>
#include <vector>
#include <map>
#include <cmath>
#include <cassert>

const std::string line_str = std::string{}.assign(30, '-');

struct City
{
    bool on_stack;
    int id, layer;
    double x, y, cost;

    City(int id, double x, double y)
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
        std::printf("%s\n", line_str.c_str());
        std::printf("[Info]: City details\n");
        std::printf("%s\n", line_str.c_str());
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

double distance(const City& a, const City& b)
{
    const auto x_sqr = (a.x - b.x) * (a.x - b.x);
    const auto y_sqr = (a.y - b.y) * (a.y - b.y);
    return std::sqrt(x_sqr + y_sqr);
}

double insertionCost(const City& new_city, const City& cityA, const City& cityB)
{
    return distance(new_city, cityA) + distance(new_city, cityB) -
           distance(cityA, cityB);
}

int properIndex(int index, int n)
{
    index = (index == -1) ? (n - 1) : index;
    index = (index == n) ? 0 : index;
    return index;
}

auto getNeigbhours(int index, Path_t& path)
{
    typedef std::tuple<Node_t, Node_t, Node_t> result_t;
    const int n = path.size();
    index = properIndex(index, n);
    const int prev{ properIndex(index - 1, n) },
        next{ properIndex(index + 1, n) };
    return result_t{ path[prev], path[index], path[next] };
}

auto getNeigbhours(int index, int n)
{
    typedef std::tuple<int, int, int> result_t;
    index = properIndex(index, n);
    return result_t{ properIndex(index - 1, n), index,
                     properIndex(index + 1, n) };
}

double neighbourCost(int index, Path_t& path, bool update)
{
    assert("[Error] (neighbourCost): given path size less than 3" &&
           (path.size() >= 3));
    auto [node_prev, node, node_next] = getNeigbhours(index, path);
    const double cost{ insertionCost(*node, *node_prev, *node_next) };
    if (update) {
        node->cost = cost;
    }
    return cost;
}

void removeNode(int index, Path_t& path)
{
    path[index]->on_stack = true;
    path[index]->print();
    path.erase(path.begin() + index);
    neighbourCost(index, path, true);
    neighbourCost(index - 1, path, true);
}

double getArea(const City& a, const City& b, const City& c)
{
    return std::abs(a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y)) /
           2.0;
}
bool isCollinearAndBetween(const City& c, const City& a, const City& b) {
    return std::min(a.x, b.x) <= c.x && c.x <= std::max(a.x, b.x) &&
           std::min(a.y, b.y) <= c.y && c.y <= std::max(a.y, b.y);
}
bool isInside(const Node_t& node, const Node_t& nodeA, const Node_t& nodeB,
              const Node_t& nodeC)
{
    constexpr double tolerance = 1e-8;
    const double area_total = getArea(*nodeA, *nodeB, *nodeC);

    const double area1 = getArea(*node, *nodeA, *nodeB);
    const double area2 = getArea(*node, *nodeB, *nodeC);
    const double area3 = getArea(*node, *nodeC, *nodeA);

    if (std::abs(area1) < tolerance) {
        return isCollinearAndBetween(*node, *nodeA, *nodeB);
    }
    if (std::abs(area2) < tolerance) {
        return isCollinearAndBetween(*node, *nodeB, *nodeC);
    }
    if (std::abs(area3) < tolerance) {
        return isCollinearAndBetween(*node, *nodeC, *nodeA);
    }

    const double area_diff = std::abs(area_total - (area1 + area2 + area3));
    return (area_diff < tolerance);
}

void removeIntersection(int index, Path_t& path)
{
    int num_nodes = path.size();
    auto [node_prev, node_curr, node_next] = getNeigbhours(index, path);
    auto [idx_prev, idx_curr, idx_next] = getNeigbhours(index, num_nodes);
    for (int idx{ 0 }; idx != num_nodes;) {
        if ((idx == idx_prev) or (idx == idx_curr) or (idx == idx_next)) {
            ++idx;
            continue;
        }
        const Node_t& node{path[idx]};
        if (not isInside(node, node_prev, node_curr, node_next)) {
            ++idx;
            continue;
        }
        path[idx]->print();
        removeNode(idx, path);
        node_prev->print();
        node_curr->print();
        node_next->print();
        --num_nodes;
        const std::string& msg{ "path size (" + std::to_string(path.size()) +
                                ") num_nodes (" + std::to_string(num_nodes) +
                                ") index (" + std::to_string(index) +
                                ") idx (" + std::to_string(idx) + ")." };
        std::cout << ("[Info] (removeIntersection): " + msg + "\n");
        if (idx < index) {
            --index;
            std::tie(node_prev, node_curr, node_next) = getNeigbhours(index, path);
            std::tie(idx_prev, idx_curr, idx_next) = getNeigbhours(index, num_nodes);
        }
    }
    // const std::string& msg{ "path size (" + std::to_string(path.size()) + ")" };
    // std::cout << ("[Info] (removeIntersection): " + msg + "\n");
}

int findBestInsertion(const City& new_city, const Path_t& path)
{
    const int num_nodes = path.size();
    assert("[Error] (findBestInsertion): given path size less than 2" &&
           (num_nodes >= 2));
    double min_cost{ insertionCost(new_city, *path.back(), *path[0]) };
    int best_index{ 0 };

    for (int idx{ 0 }; idx != num_nodes; ++idx) {
        if (new_city.id == path[idx]->id) {
            std::printf(
                "[Error] (findBestInsertion): trying to add city that already exists in path.");
            new_city.print();
            exit(EXIT_FAILURE);
        }
        if (idx == (num_nodes - 1)) {
            continue;
        }
        const double cost = insertionCost(new_city, *path[idx], *path[idx + 1]);
        if (cost < min_cost) {
            min_cost = cost;
            best_index = idx + 1;
        }
    }

    assert("[Error] (findBestInsertion): invalid index position found" &&
           (best_index >= 0));

    return best_index;
}

void updatePath(Path_t& path, const Node_t& node)
{
    int position = findBestInsertion(*node, path);
    path.insert(path.begin() + position, node);
    neighbourCost(position, path, true);
    neighbourCost(position - 1, path, true);
    neighbourCost(position + 1, path, true);

    path[position]->print();
    node->print();
    const std::string& msg{ "path size (" + std::to_string(path.size()) +
                            ") position : " + std::to_string(position) };
    std::cout << ("[Info] (updatePath): " + msg + "\n");
    removeIntersection(position, path);
}

bool validateNode(int index, const Path_t& path)
{
    const int num_nodes = path.size();
    const Node_t& node{ path[index] };
    const double cost_current{ node->cost };
    if (index != 0 and index != (num_nodes - 1)) {
        if (insertionCost(*node, *path.back(), *path[0]) < cost_current) {
            return false;
        }
    }
    for (int idx{ 0 }; idx != (num_nodes - 1); ++idx) {
        if ((idx == index) or (idx + 1 == index)) {
            continue;
        }
        if (insertionCost(*node, *path[idx], *path[idx + 1]) < cost_current) {
            return false;
        }
    }
    return true;
}

std::optional<Node_t> validatePath(Path_t& path)
{
    std::optional<Node_t> node_erased{ std::nullopt };
    int num_nodes = path.size();
    for (int idx{ 0 }; idx != num_nodes;) {
        if (validateNode(idx, path)) {
            ++idx;
            // const std::string& msg{ "idx (" + std::to_string(idx) + ")." };
            // std::cout << ("[Info] (validatePath): " + msg + "\n");
            continue;
        }
        if (node_erased.has_value()) {
            node_erased = (path[idx]->id < (*node_erased)->id) ? path[idx] :
                                                                 node_erased;
        } else {
            node_erased = path[idx];
        }
        path[idx]->print();
        removeNode(idx, path);
        --num_nodes;
        idx = 0;
        const std::string& msg{ "path size (" + std::to_string(path.size()) +
                                ") num_nodes (" + std::to_string(num_nodes) +
                                ") idx (" + std::to_string(idx) + ")." };
        std::cout << ("[Info] (validatePath): " + msg + "\n");
    }
    return node_erased;
}

void initializePath(Path_t& path, Cities_t& cities)
{
    const int num_cities = cities.size();
    assert("[Error] (initializePath): given number of cities less than 1" &&
           (num_cities > 0));
    path.clear();
    path.reserve(num_cities);
    const Node_t it_begin{ cities.begin() };
    if (num_cities < 4) {
        for (Node_t it{ it_begin }; it != cities.end(); ++it) {
            path.push_back(it);
        }
    }
}

void constructPath(Path_t& path, Cities_t& cities, int n)
{
    const int num_cities = cities.size();
    const Node_t it_begin{ cities.begin() };
    if (n < 1) {
        n = 0.10f * num_cities;
        n = (n == 0) ? 3 : n;
    }
    for (int idx{ 0 }; idx != n; ++idx) {
        Node_t it{ it_begin + idx };
        path.push_back(it);
        it->on_stack = false;
    }
    assert("[Error] (constructPath): constructed path size less than 3" &&
           (path.size() >= 3));
}

template <typename T>
void constructPath(Path_t& path, Cities_t& cities, int n, T rng)
{
    const int num_cities = cities.size();
    const Node_t it_begin{ cities.begin() };
    std::vector<int> indices(num_cities, 0);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);
    if (n < 1) {
        n = 0.10f * num_cities;
        n = (n == 0) ? 3 : n;
    }
    for (int idx{ 0 }; idx != n; ++idx) {
        Node_t it{ it_begin + indices[idx] };
        path.push_back(it);
        it->on_stack = false;
    }
    assert("[Error] (constructPath): constructed path size less than 3" &&
           (path.size() >= 3));
}

void runDiscreteENN(Cities_t& stack, Path_t& path)
{
    const int num_cities = stack.size();
    // int count_cities{0};
    for (Node_t it{ stack.begin() }; it != stack.end();) {
        // if ((it + 1) == stack.end() and (count_cities + 1 != num_cities)) {
        //     it = stack.begin();
        // }
        bool& on_stack{ it->on_stack };
        if (not on_stack) {
            ++it;
            continue;
        }
        on_stack = false;

        std::cout << ("[Info] (runDiscreteENN): updatePath\n");
        updatePath(path, it);

        std::cout << ("[Info] (runDiscreteENN): validatePath\n");
        validatePath(path);
        it = stack.begin();

        // std::cout << ("[Info] (runDiscreteENN): updatePath\n");
        // updatePath(path, it);

        // std::cout << ("[Info] (runDiscreteENN): validatePath\n");
        // std::optional<Node_t> node_erased = validatePath(path);
        // if (node_erased.has_value()) {
        //     it = *node_erased;
        //     continue;
        // }
        // ++it;
    }
}

bool isNumber(const std::string& str)
{
    std::istringstream is(str);
    double number{ 0.00 };
    is >> std::noskipws >> number;
    return is.eof() and !is.fail();
}

void parseCities(Cities_t& cities, const std::string& filename)
{
    cities.clear();
    cities.reserve(200);
    std::ifstream in_file{ filename };
    if (not in_file) {
        std::cout << "[Error] (parseCities): Could not open " << filename
                  << ". Exiting\n";
        exit(EXIT_FAILURE);
    }
    std::string line;
    int line_count{ 0 }, line_num{ 0 };
    double x, y;
    bool start_parse{ false };
    while (std::getline(in_file, line)) {
        std::string word;
        std::istringstream tmp_stream(line);
        tmp_stream >> word;
        start_parse = isNumber(word);

        if (start_parse) {
            ++line_count;
            std::istringstream line_stream(line);
            line_stream >> line_num >> x >> y;
            City city(-1, x, y);
            cities.push_back(city);
        }
    }
    assert(
        "[Error] (parseCities): Number of coords parsed not equal to number of coords" &&
        (line_count == line_num));
    std::cout << "[Info] (parseCities): " << filename
              << " parsed. Total number of coords : " << line_count << '\n';
}

struct MinMaxCoords
{
    double min_x{ std::numeric_limits<double>::max() };
    double max_x{ -1.0 };
    double min_y{ std::numeric_limits<double>::max() };
    double max_y{ -1.0 };
};

void getMinMaxCoords(const Cities_t& cities, MinMaxCoords& a)
{
    double& min_x{ a.min_x };
    double& max_x{ a.max_x };
    double& min_y{ a.min_y };
    double& max_y{ a.max_y };
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

std::pair<int, int> createCityLayers(const Cities_t& cities,
                                     CityLayers_t& city_layers,
                                     MinMaxCoords minmax_coords, int depth = 1)
{
    if (cities.empty()) {
        return std::make_pair(0, depth);
    }
    const double& min_x{ minmax_coords.min_x };
    const double& max_x{ minmax_coords.max_x };
    const double& min_y{ minmax_coords.min_y };
    const double& max_y{ minmax_coords.max_y };
    double mid_x = (min_x + max_x) / 2.0;
    double mid_y = (min_y + max_y) / 2.0;

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

    int added{ 0 }, sum{ 0 };
    for (int idx{ 0 }; idx != parts; ++idx) {
        Cities_t& cities_quadrant{ quadrants[idx] };
        // std::printf("[Info]: Layer : %d Quadrant %d size : %zu.\n", depth, idx, cities_quadrant.size());
        if (cities_quadrant.empty())
            continue;
        City city = cities_quadrant.back();
        cities_quadrant.pop_back();
        city_layers[depth].push_back(city);
        ++added;
        sum += cities_quadrant.size();
        // std::printf("[Info]: Current layer : %d, quadrant : %d, and quadrant size %zu.\n", depth, idx,
        //             quadrants[idx].size());
    }
    // std::printf("[Info]: Current layer : %d and layer size %d.\n", depth,
    //             added);
    // std::printf("[Info]: Coords in all quadrants : %d.\n", sum);
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
    // std::printf("[Info]: Current layer : %d and layer size %zu.\n", depth,
    //             city_layers[depth].size());
    return std::make_pair(added, depth);
}

void createStack(const Cities_t& cities, Cities_t& stack, int& layers)
{
    CityLayers_t city_layers;
    MinMaxCoords minmax_coords;
    getMinMaxCoords(cities, minmax_coords);
    const auto [num_cities, depth] =
        createCityLayers(cities, city_layers, minmax_coords);
    if (cities.size() != static_cast<std::size_t>(num_cities)) {
        const std::string& error_msg{ "number of cities in city_layers (" +
                                      std::to_string(num_cities) +
                                      ") doesn't match cities size (" +
                                      std::to_string(cities.size()) + ")." };
        std::cerr << ("[Error] (createStack): " + error_msg + " Exiting.\n");
        exit(EXIT_FAILURE);
    }
    layers = depth;
    int count{ 0 };
    for (int idx{ 1 }; idx <= depth; ++idx) {
        const Cities_t& city_layer{ city_layers[idx] };
        const int layer_size = city_layer.size();
        for (int jdx{ 0 }; jdx != layer_size; ++jdx) {
            City city = city_layer[jdx];
            city.on_stack = true;
            city.id = count;
            city.layer = idx;
            ++count;
            stack.push_back(city);
        }
    }
    if (stack.size() != cities.size()) {
        const std::string& error_msg{ "stack size (" +
                                      std::to_string(stack.size()) +
                                      ") doesn't match number of cities (" +
                                      std::to_string(cities.size()) + ")." };
        std::cerr << ("[Error] (createStack): " + error_msg + " Exiting.\n");
        exit(EXIT_FAILURE);
    }
}
