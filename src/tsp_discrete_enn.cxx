// -*- C++ -*-

#include "tsp_discrete_enn.hpp"

#include <fstream>
#include <sstream>

#include "utils.hpp"

#ifdef TSP_DRAW
#include <SFML/Graphics.hpp>
#endif

void parseCities(Cities_t& cities, const std::string& filename)
{
    cities.clear();
    cities.reserve(200);
    std::ifstream in_file{ filename };
    if (not in_file) {
        utils::printErr("Could not open " + filename + "\nExiting.",
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
    int depth_current{depth};
    {
        auto [tmp_add, tmp_depth] = createCityLayers(
            quadrants[0], city_layers, { min_x, mid_x, min_y, mid_y }, depth_current);
        added += tmp_add;
        depth = std::max(tmp_depth, depth);
    }
    {
        auto [tmp_add, tmp_depth] = createCityLayers(
            quadrants[1], city_layers, { mid_x, max_x, min_y, mid_y }, depth_current);
        added += tmp_add;
        depth = std::max(tmp_depth, depth);
    }
    {
        auto [tmp_add, tmp_depth] = createCityLayers(
            quadrants[2], city_layers, { min_x, mid_x, mid_y, max_y }, depth_current);
        added += tmp_add;
        depth = std::max(tmp_depth, depth);
    }
    {
        auto [tmp_add, tmp_depth] = createCityLayers(
            quadrants[3], city_layers, { mid_x, max_x, mid_y, max_y }, depth_current);
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

void createStack(const Cities_t& cities, Cities_t& stack)
{
    const Index_t num_cities = cities.size();
    stack.reserve(static_cast<std::size_t>(num_cities));
    auto sorter = [](const City& a, const City& b) { return (std::pair{a.x, a.y} < std::pair{b.x, b.y}); };
    auto cities_copy = cities;
    std::sort(cities_copy.begin(), cities_copy.end(), sorter);
    bool change{ false };
    Index_t start{0}, end{num_cities - 1};
    for (Index_t idx{0}; idx != num_cities; ++idx) {
        City city;
        if (change) {
            city = cities_copy[end];
            --end;
        } else {
            city = cities_copy[start];
            ++start;
        }
        city.on_stack = true;
        city.id = idx;
        stack.push_back(city);
        change = not change;
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

void createStack(const Cities_t& cities, Cities_t& stack, int parts)
{
    const Index_t num_cities = cities.size();
    stack.reserve(static_cast<std::size_t>(num_cities));
    auto sorter = [](const City& a, const City& b) { return (std::pair{a.x, a.y} < std::pair{b.x, b.y}); };
    auto cities_copy = cities;
    std::sort(cities_copy.begin(), cities_copy.end(), sorter);
    MinMaxCoords minmax_coords;
    minmax_coords.update(cities);
    const auto [min_x, max_x, min_y, max_y] = minmax_coords.value();
    const auto [min_coord, max_coord] = minmax_coords.minmax();
    min_x = min_y = static_cast<Value_t>(0);
    max_x = max_y = max_coord;
    Value_t mid_x = (min_x + max_x) / 2.0;
    Value_t mid_y = (min_y + max_y) / 2.0;
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

    int count{ 0 };
    for (int idx{ 0 }; idx != parts; ++idx) {
        Cities_t& cities_quadrant{ quadrants[idx] };
        if (cities_quadrant.empty()) {
            continue;
        }
        City city = cities_quadrant.back();
        cities_quadrant.pop_back();
        city.id = count;
        city.on_stack = true;
        stack.push_back(city);
        ++count;
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


Value_t insertionCost(const City& new_city, const City& cityA,
                             const City& cityB)
{
    const Value_t cost = getDistance(new_city, cityA) +
                         getDistance(new_city, cityB) -
                         getDistance(cityA, cityB);
    const Value_t cost_rounded{ utils::getRound3(cost) };
    return utils::isEqual(cost_rounded, VALUE_ZERO) ? VALUE_ZERO : cost_rounded;
}

bool isInside(const City& city, const City& cityA, const City& cityB,
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

bool hasIntersection(const City& cityA1, const City& cityA2,
                            const City& cityB1, const City& cityB2)
{
    // if (sameCoords(cityA1, cityB1) or sameCoords(cityA1, cityB2) or
    //     sameCoords(cityA2, cityB1) or sameCoords(cityA2, cityB2)) {
    //     return false;
    // }
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

void fitPointsInWindow(MinMaxCoords& minmax_coords_copy, const sf::Vector2u& window_size,
                       Value_t margin)
{
    const auto [min_x_copy, max_x_copy, min_y_copy, max_y_copy] = minmax_coords_copy.value();

    Value_t width{ max_x_copy - min_x_copy };
    Value_t height{ max_y_copy - min_y_copy };
    // Calculate scale factors for x and y to fit the plot within the window, considering margins
    Value_t scaleX = (window_size.x - 2 * margin) / width;
    Value_t scaleY = (window_size.y - 2 * margin) / height;

    // Apply scale and translation to the points
    min_x_copy = (min_x_copy - min_x_copy) * scaleX + margin;
    min_y_copy = (min_y_copy - min_y_copy) * scaleY + margin;
    max_x_copy = (max_x_copy - min_x_copy) * scaleX + margin;
    max_y_copy = (max_y_copy - min_y_copy) * scaleY + margin;
}

void drawPath(const Indices_t& path, const Cities_t& cities, bool show_coords,
              const std::string& title, float close_time, int highlight, int)
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
    fitPointsInWindow(points, window.getSize(), minmax_coords, 50);
    Cities_t cities_copy = cities;
    fitPointsInWindow(cities_copy, window.getSize(), minmax_coords, 50);
    sf::ConvexShape polygon;
    polygon.setPointCount(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        polygon.setPoint(i, points[i]);
    }

    polygon.setFillColor(sf::Color::Transparent); // Transparent fill
    polygon.setOutlineThickness(2); // Line thickness
    polygon.setOutlineColor(sf::Color::Red); // Line color

    // auto minmax_coords_copy = minmax_coords;
    // const auto [min_x, max_x, min_y, max_y] = minmax_coords_copy.value();
    // min_x = min_y = static_cast<Value_t>(0);
    // max_x = max_y = max_coord;
    // fitPointsInWindow(minmax_coords_copy, window.getSize(), 50);
    MinMaxCoords minmax_coords_copy;
    minmax_coords_copy.update(cities_copy);
    const auto [min_x, max_x, min_y, max_y] = minmax_coords_copy.value();
    auto mid_x{(max_x + min_x)/2}, mid_y{(max_y + min_y)/2};
    // Define a vector of pairs to hold point coordinates
    std::vector<std::pair<sf::Vector2f, sf::Vector2f>> linePoints = {
        {{min_x, mid_y}, {max_x, mid_y}},
        {{mid_x, min_y}, {mid_x, max_y}},
        {{mid_x/2, min_y}, {mid_x/2, max_y}},
        {{mid_x - mid_x/2, min_y}, {mid_x - mid_x/2, max_y}},
        {{mid_x + mid_x/2, min_y}, {mid_x + mid_x/2, max_y}},
        {{min_x, mid_y - mid_y/2}, {max_x, mid_y - mid_y/2}},
        {{min_x, mid_y + mid_y/2}, {max_x, mid_y + mid_y/2}},
    };

    sf::Clock clock;
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
        // Check the clock if 10 seconds have passed
        if (clock.getElapsedTime().asSeconds() >= close_time) {
            std::cout << "10 seconds have elapsed. Closing window."
                      << std::endl;
            window.close();
        }

        window.clear(sf::Color::White); // Clear the window
        window.setTitle(title);

        // Draw each line segment
        int change_color{0};
        for (const auto& line : linePoints) {
            sf::VertexArray lineSegment(sf::Lines, 2);
            lineSegment[0].position = line.first;
            lineSegment[1].position = line.second;
            if (change_color < 2) {
                lineSegment[0].color = sf::Color::Black;
                lineSegment[1].color = sf::Color::Black;
            } else {
                lineSegment[0].color = sf::Color::Magenta;
                lineSegment[1].color = sf::Color::Magenta;
            }

            window.draw(lineSegment);
            ++change_color;
        }

        // for (const auto& point : points) {
        //     sf::CircleShape marker(5); // Small circle with radius 5 pixels
        //     marker.setFillColor(sf::Color::Blue); // Use a contrasting color
        //     marker.setPosition(
        //         point.x,
        //         point.y); // Adjust position to center the marker on the vertex

        //     window.draw(marker);
        // }
        for (const City& city : cities_copy) {
            if (show_coords) {
                sf::Text label;
                label.setFont(font); // Set the font
                label.setString(std::to_string(
                    city.id)); // Set the point's ID as the label text
                label.setCharacterSize(14); // Set the character size
                label.setFillColor(sf::Color::Black); // Set the text color
                if (city.id == highlight) {
                    label.setFillColor(sf::Color::Magenta);
                }
                label.setPosition(city.x + 10,
                                  city.y -
                                      10); // Position the label near the marker
                window.draw(label); // Draw the label
            }

            sf::CircleShape marker(5); // Small circle with radius 5 pixels
            if (city.on_stack) {
                marker.setFillColor(
                    sf::Color::Green); // Use a contrasting color
            } else {
                marker.setFillColor(sf::Color::Blue); // Use a contrasting color
            }
            if (city.id == highlight) {
                marker.setFillColor(sf::Color::Magenta);
            }
            marker.setPosition(
                city.x,
                city.y); // Adjust position to center the marker on the vertex

            window.draw(marker);
        }
        window.draw(polygon); // Draw the polygon

        window.display(); // Display what was drawn
    }
}

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
