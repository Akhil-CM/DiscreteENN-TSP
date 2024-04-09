// -*- C++ Header -*-
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <random>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <optional>
#include <vector>
#include <map>
#include <cmath>
#include <cassert>

#include <SFML/Graphics.hpp>

#define DEBUG_PRINT 0

const std::string line_str = std::string{}.assign(30, '-');

bool isNumber(const std::string& str)
{
    std::istringstream is(str);
    double number{ 0.00 };
    is >> std::noskipws >> number;
    return is.eof() and !is.fail();
}

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
typedef std::optional<Node_t> NodeOpt_t;

template <typename T> class MatchItem
{
public:
    MatchItem(T item)
        : m_item{ item }
    {
    }
    bool operator()(const T& item)
    {
        return (item == m_item);
    }
    T& item()
    {
        return m_item;
    }

    const T& item() const
    {
        return m_item;
    }

private:
    T m_item;
};

template <> bool MatchItem<City>::operator()(const City& city)
{
    return (city.id == m_item.id);
}

template <> bool MatchItem<Node_t>::operator()(const Node_t& node)
{
    return (node->id == m_item->id);
}

int findNode(Node_t node, const Path_t& path)
{
    const Path_t::const_iterator& it_begin{ path.begin() };
    const Path_t::const_iterator& it_end{ path.end() };
    const Path_t::const_iterator& it =
        std::find_if(it_begin, it_end, MatchItem{ node });
    return (it == it_end) ? -1 : std::distance(it_begin, it);
}

void removeNode(std::size_t index, Path_t& path)
{
    path[index]->on_stack = true;
#if (DEBUG_PRINT > 0)
    path[index]->print();
#endif
    path.erase(path.begin() + static_cast<int>(index));
}

void removeNode(Node_t node, Path_t& path)
{
    const int index = findNode(node, path);
    if (index == -1) {
        std::cerr
            << "[Error] (removeNode): Request to remove a node not present in path. Exiting\n";
        exit(EXIT_FAILURE);
    }
    removeNode(static_cast<std::size_t>(index), path);
}

int properIndex(int index, int n)
{
    index = (index < 0) ? (n + index) : index;
    index = (index > (n - 1)) ? (index - n) : index;
    return index;
}

auto getNeigbhours(int index, int n)
{
    typedef std::tuple<int, int, int> result_t;
    const int curr = properIndex(index, n);
    return result_t{ properIndex(curr - 1, n), curr, properIndex(curr + 1, n) };
}

auto getNeigbhours(int index, const Path_t& path)
{
    typedef std::tuple<Node_t, Node_t, Node_t> result_t;
    const int n = path.size();
    const auto [prev, curr, next] = getNeigbhours(index, n);
    return result_t{ path[prev], path[curr], path[next] };
}

struct MinMaxCoords
{
    double min_x{ std::numeric_limits<double>::max() };
    double max_x{ -1.0 };
    double min_y{ std::numeric_limits<double>::max() };
    double max_y{ -1.0 };
};

void getMinMaxCoords(const Path_t& path, MinMaxCoords& a)
{
    double& min_x{ a.min_x };
    double& max_x{ a.max_x };
    double& min_y{ a.min_y };
    double& max_y{ a.max_y };
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

double nodeCost(int index, const Path_t& path, bool update)
{
    const int num_nodes = path.size();
    if (num_nodes < 3) {
        std::cerr
            << "[Error] (nodeCost): given path size less than 3. Exiting\n";
        exit(EXIT_FAILURE);
    }
    if (index >= num_nodes) {
        std::cerr
            << "[Error] (nodeCost): given index size greater than number of nodes in path. Exiting\n";
        exit(EXIT_FAILURE);
    }
    const auto [node_prev, node, node_next] = getNeigbhours(index, path);
    const double cost{ insertionCost(*node, *node_prev, *node_next) };
    if (update) {
        node->cost = cost;
    }
    return cost;
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
bool isInside(Node_t node, Node_t nodeA, Node_t nodeB, Node_t nodeC)
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

NodeOpt_t removeIntersection(Node_t node_prev, Node_t node_curr,
                             Node_t node_next, Path_t& path, bool recursive)
{
    NodeOpt_t node_erased{ std::nullopt };
    int num_nodes = path.size();
    const int id_prev{ node_prev->id };
    const int id_curr{ node_curr->id };
    const int id_next{ node_next->id };
    for (int idx{ 0 }; idx != num_nodes;) {
#if (DEBUG_PRINT > 1)
        const std::string& msg_first{ "idx (" + std::to_string(idx) +
                                      ") num_nodes (" +
                                      std::to_string(num_nodes) +
                                      ") path size (" +
                                      std::to_string(path.size()) + ")." };
        std::cout << ("[Debug] (removeIntersection): " + msg_first + "\n");
#endif
        const Node_t node{ path[static_cast<std::size_t>(idx)] };
        const int id{ node->id };
        if ((id == id_prev) or (id == id_curr) or (id == id_next)) {
            ++idx;
            continue;
        }
        if (not isInside(node, node_prev, node_curr, node_next)) {
            ++idx;
            continue;
        }
#if (DEBUG_PRINT > 1)
        std::cout << ("[Debug] (removeIntersection): Done 1\n");
#endif
        removeNode(node, path);
        --num_nodes;
#if (DEBUG_PRINT > 1)
        std::cout << ("[Debug] (removeIntersection): Done 2\n");
#endif
        if (node_erased.has_value()) {
            node_erased = (node->id < (*node_erased)->id) ? node : node_erased;
        } else {
            node_erased = node;
        }
#if (DEBUG_PRINT > 0)
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
        const std::string& msg{ "path size (" + std::to_string(path.size()) +
                                ") num_nodes (" + std::to_string(num_nodes) +
                                ") idx (" + std::to_string(idx) + ")." };
        std::cout << ("[Debug] (removeIntersection): " + msg + "\n");
#endif
        if (recursive) {
            const int idx_next{ properIndex(idx, num_nodes) };
            const int idx_prev{ properIndex(idx_next - 1, num_nodes) };
            NodeOpt_t node_erased_tmp = removeIntersection(
                path[idx_prev], node, path[idx_next], path, recursive);
            if (node_erased_tmp.has_value()) {
                node_erased = ((*node_erased_tmp)->id < (*node_erased)->id) ?
                                  node_erased_tmp :
                                  node_erased;
            }
            num_nodes = path.size();
            idx = 0;
        }
    }
    return node_erased;
}

template <typename T>
bool vectContains(const T& item, const std::vector<T>& vect)
{
    for (const T& elmnt : vect) {
        if (elmnt == item) {
            return true;
        }
    }
    return false;
}

int findBestInsertion(Node_t node, const Path_t& path)
{
    const int num_nodes = path.size();
    assert("[Error] (findBestInsertion): given path size less than 2" &&
           (num_nodes >= 2));

    const City& new_city{ *node };

    double min_cost{ insertionCost(new_city, *(path[num_nodes - 1]),
                                   *(path[0])) };
    int best_index{ 0 };

    for (int idx{ 0 }; idx != num_nodes; ++idx) {
        if (new_city.id == path[idx]->id) {
            std::cerr
                << ("[Error] (findBestInsertion): trying to add city that already exists in path.\n");
            new_city.print();
            exit(EXIT_FAILURE);
        }
        if (idx == (num_nodes - 1)) {
            continue;
        }
        const double cost =
            insertionCost(new_city, *(path[idx]), *(path[idx + 1]));
        if (cost < min_cost) {
            min_cost = cost;
            best_index = idx + 1;
        }
    }

    assert("[Error] (findBestInsertion): invalid index position found" &&
           (best_index != -1));

    return best_index;
}

int addNode(Node_t node, Path_t& path)
{
    int position = findBestInsertion(node, path);
    node->on_stack = false;
    path.insert(path.begin() + position, node);
#if (DEBUG_PRINT > 0)
    path[0]->print();
    std::cout << "FIRST CITY\n";
    node->print();
    std::cout << "Added City\n";
    const std::string& msg{ "path size (" + std::to_string(path.size()) +
                            ") position : " + std::to_string(position) };
    std::cout << ("[Debug] (addNode): " + msg + "\n");
#endif

    return position;
}

void updateCostNeighbour(int index, const Path_t& path)
{
    const int idx_curr = properIndex(index, path.size());
    nodeCost(idx_curr, path, true);
    nodeCost(idx_curr - 1, path, true);
    nodeCost(idx_curr + 1, path, true);
}

void updateCostAll(const Path_t& path)
{
    const int num_nodes = path.size();
    for (int idx{ 0 }; idx != num_nodes; ++idx) {
        nodeCost(idx, path, true);
    }
}

bool validateNode(int index, const Path_t& path)
{
    const int num_nodes = path.size();
    const Node_t& node{ path[index] };
    const double cost_current{ node->cost };
    if (index != 0 and index != (num_nodes - 1)) {
        if (insertionCost(*node, *(path[num_nodes - 1]), *(path[0])) <
            cost_current) {
            return false;
        }
    }
    for (int idx{ 0 }; idx != (num_nodes - 1); ++idx) {
        if ((idx == index) or (idx + 1 == index)) {
            continue;
        }
        if (insertionCost(*node, *(path[idx]), *(path[idx + 1])) <
            cost_current) {
            return false;
        }
    }
    return true;
}

NodeOpt_t validatePath(Path_t& path, bool intersection, bool recursive)
{
    NodeOpt_t node_erased{ std::nullopt };
    int num_nodes = path.size();
    for (int idx{ 0 }; idx != num_nodes;) {
#if (DEBUG_PRINT > 1)
        const std::string& msg{ "idx (" + std::to_string(idx) +
                                ") num_nodes (" + std::to_string(num_nodes) +
                                ") path size (" + std::to_string(path.size()) +
                                ")." };
        std::cout << ("[Debug] (validatePath): " + msg + "\n");
#endif
        if (validateNode(idx, path)) {
            ++idx;
            continue;
        }
#if (DEBUG_PRINT > 1)
        std::cout << ("[Debug] (validatePath): Done 1\n");
#endif
        const Node_t node{ path[static_cast<std::size_t>(idx)] };
        removeNode(node, path);
        --num_nodes;
        if (node_erased.has_value()) {
            node_erased = (node->id < (*node_erased)->id) ? node : node_erased;
        } else {
            node_erased = node;
        }
#if (DEBUG_PRINT > 1)
        std::cout << ("[Debug] (validatePath): Done 2\n");
#endif
        if (intersection) {
            const int idx_next{ properIndex(idx, num_nodes) };
            const int idx_prev{ properIndex(idx_next - 1, num_nodes) };
            NodeOpt_t node_erased_tmp = removeIntersection(
                path[idx_prev], node, path[idx_next], path, recursive);
            if (node_erased_tmp.has_value()) {
                node_erased = ((*node_erased_tmp)->id < (*node_erased)->id) ?
                                  node_erased_tmp :
                                  node_erased;
            }
            num_nodes = path.size();
        }
#if (DEBUG_PRINT > 1)
        std::cout << ("[Debug] (validatePath): Done 3\n");
#endif
        idx = 0;
        updateCostAll(path);
#if (DEBUG_PRINT > 1)
        std::cout << ("[Debug] (validatePath): Done 4\n");
#endif
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
    if (n < 1) {
        n = 0.10f * num_cities;
        n = (n == 0) ? 3 : n;
    }
    const Node_t it_begin{ cities.begin() };
    for (int idx{ 0 }; idx != n; ++idx) {
        Node_t it{ it_begin + idx };
        it->on_stack = false;
        path.push_back(it);
    }
    updateCostAll(path);
    assert("[Error] (constructPath): constructed path size less than 3" &&
           (path.size() >= 3));
}

template <typename T>
void constructPath(Path_t& path, Cities_t& cities, int n, T rng)
{
    const int num_cities = cities.size();
    std::vector<int> indices(num_cities, 0);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);
    if (n < 1) {
        n = 0.10f * num_cities;
        n = (n == 0) ? 3 : n;
    }
    const Node_t it_begin{ cities.begin() };
    for (int idx{ 0 }; idx != n; ++idx) {
        Node_t it{ it_begin + indices[idx] };
        it->on_stack = false;
        path.push_back(it);
    }
    updateCostAll(path);
    assert("[Error] (constructPath): constructed path size less than 3" &&
           (path.size() >= 3));
}

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
                       MinMaxCoords minmax_coords, float margin)
{
    if (points.empty())
        return; // Check if the vector is not empty

    float minX = minmax_coords.min_x, maxX = minmax_coords.max_x;
    float minY = minmax_coords.min_y, maxY = minmax_coords.max_y;

    // Calculate scale factors for x and y to fit the plot within the window, considering margins
    float scaleX = (windowSize.x - 2 * margin) / (maxX - minX);
    float scaleY = (windowSize.y - 2 * margin) / (maxY - minY);
    float scale = std::min(
        scaleX,
        scaleY); // Use the smaller scale factor to maintain aspect ratio

    // Calculate translation to center the plot
    float translateX = (windowSize.x - (maxX + minX) * scale) / 2.0f;
    float translateY = (windowSize.y - (maxY + minY) * scale) / 2.0f;

    // Apply scale and translation to the points
    for (auto& point : points) {
        point.x = (point.x - minX) * scale + margin + translateX;
        point.y = (point.y - minY) * scale + margin + translateY;
    }
}

void fitPointsInWindow(Cities_t& cities, const sf::Vector2u& windowSize,
                       MinMaxCoords minmax_coords, float margin)
{
    if (cities.empty())
        return; // Check if the vector is not empty

    double minX = minmax_coords.min_x, maxX = minmax_coords.max_x;
    double minY = minmax_coords.min_y, maxY = minmax_coords.max_y;

    // Calculate scale factors for x and y to fit the plot within the window, considering margins
    float scaleX = (windowSize.x - 2 * margin) / (maxX - minX);
    float scaleY = (windowSize.y - 2 * margin) / (maxY - minY);
    float scale = std::min(
        scaleX,
        scaleY); // Use the smaller scale factor to maintain aspect ratio

    // Calculate translation to center the plot
    float translateX = (windowSize.x - (maxX + minX) * scale) / 2.0f;
    float translateY = (windowSize.y - (maxY + minY) * scale) / 2.0f;

    // Apply scale and translation to the points
    for (auto& city : cities) {
        city.x = (city.x - minX) * scale + margin + translateX;
        city.y = (city.y - minY) * scale + margin + translateY;
    }
}

void drawPath(const Path_t& path, const Cities_t& stack)
{
    // Get the current desktop video mode
    [[maybe_unused]] sf::VideoMode desktopMode =
        sf::VideoMode::getDesktopMode();
    // Create a full-screen window using the current desktop resolution
    MinMaxCoords minmax_coords;
    getMinMaxCoords(stack, minmax_coords);
    // sf::RenderWindow window(sf::VideoMode(minmax_coords.max_x,
    //                                       minmax_coords.max_y),
    //                         "Polygon Plot. Press Any Key to Close");
    // sf::RenderWindow window(
    // sf::VideoMode(minmax_coords.max_x, minmax_coords.max_y),
    // "Polygon Plot. Press Any Key to Close", sf::Style::Fullscreen);
    // sf::RenderWindow window(sf::VideoMode(800, 600), "Polygon Plot. Press Any Key to Close", sf::Style::Fullscreen);
    // sf::RenderWindow window(sf::VideoMode(800, 600), "Polygon Plot. Press Any Key to Close");
    sf::RenderWindow window(desktopMode, "Polygon Plot. Press Any Key to Close",
                            sf::Style::Fullscreen);

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
            sf::Text label;
            label.setFont(font); // Set the font
            label.setString(std::to_string(
                city.id)); // Set the point's ID as the label text
            label.setCharacterSize(14); // Set the character size
            label.setFillColor(sf::Color::White); // Set the text color
            label.setPosition(
                city.x + 10, city.y - 10); // Position the label near the marker
            window.draw(label); // Draw the label

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

void runDiscreteENN(Cities_t& stack, Path_t& path,
                    [[maybe_unused]] std::default_random_engine& gen)
{
    int iter_count{ 0 };
    [[maybe_unused]] const int num_cities = stack.size();
    const Node_t& it_begin{ stack.begin() };
    const Node_t& it_end{ stack.end() };
    Node_t it{ it_begin };
    std::vector<int> indices;
    bool randomize_node{ true };
    constexpr bool recursive{ false };
    while (true) {
        if (it == it_end) {
            const bool finished = path.size() ==
                                  static_cast<std::size_t>(num_cities);
            if (finished) {
                break;
            } else {
                it = it_begin;
            }
        }
        if (not it->on_stack) {
            ++it;
            continue;
        }

#if (DEBUG_PRINT > 0)
        std::cout << ("\n[Debug] (runDiscreteENN): addNode\n");
#endif
        const int idx_added = addNode(it, path);
        updateCostNeighbour(idx_added, path);

#if (DEBUG_PRINT > 0)
        std::cout
            << ("\n[Debug] (runDiscreteENN): addNode removeIntersection\n");
#endif
        const auto [node_prev, node_curr, node_next] =
            getNeigbhours(idx_added, path);
        const auto& it_erased1 = removeIntersection(node_prev, node_curr,
                                                    node_next, path, recursive);

        if (it_erased1.has_value()) {
#if (DEBUG_PRINT > 0)
            std::cout << ("\n[Debug] (runDiscreteENN): addNode updateCostAll\n");
#endif
            updateCostAll(path);
        }

#if (DEBUG_PRINT > 0)
        std::cout << ("\n[Debug] (runDiscreteENN): validatePath\n");
#endif
        const auto& it_erased2 = validatePath(path, true, recursive);

#if (DEBUG_PRINT > 0)
        std::cout
            << ("\n[Debug] (runDiscreteENN): validatePath updateCostAll\n");
#endif
        updateCostAll(path);

        if (it_erased1.has_value() or it_erased2.has_value()) {
            it = it_begin;
            ++iter_count;
        }

        if (iter_count > 10000) {
            if (randomize_node) {
                // randomize_node = false;
                std::uniform_int_distribution<int> distrib(1, num_cities);
                it = it_begin + distrib(gen);
                continue;
            }
            iter_count = 0;
            // iter_count = 10001;
#if (DEBUG_PRINT > 0)
            std::cout << ("\n[Debug] (runDiscreteENN): drawPath started\n");
#endif
            // drawPath(path, stack);
#if (DEBUG_PRINT > 0)
            std::cout << ("[Debug] (runDiscreteENN): drawPath ended\n");
#endif
        }
    }
}
