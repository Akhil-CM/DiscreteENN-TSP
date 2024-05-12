// -*- C++ -*-

#include "tsp_discrete_enn.hpp"

#include <SFML/Window/WindowStyle.hpp>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "utils.hpp"

#ifdef TSP_DRAW
#include <SFML/Graphics.hpp>
#endif

void parseCities(Cities_t& cities, const std::string& filename)
{
    cities.clear();
    cities.reserve(1000);
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
        line = utils::trimLeft(line);
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
        std::string word{ line.substr(
            0, line.find_first_of(utils::Whitespace_Str)) };
        start_parse = utils::isNumber(word);
        if (start_parse) {
            std::istringstream line_stream{ line };
            line_stream >> line_num >> x >> y;
            // utils::printInfoFmt("The coords #%i are (%f, %f)", "parseCities", line_num, x, y);
            City city(line_count, x, y);
            cities.push_back(city);
            ++line_count;
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

auto makeStackLayers(DiscreteENN_TSP& enn_tsp, const Cities_t& quadrant,
                      const MinMaxCoords& minmax_coords, int depth)
{
    Cities_t& cities{ enn_tsp.cities() };
    CityLayers_t& stack{ enn_tsp.stack() };
    if (quadrant.empty()) {
        return std::make_pair(0, -1);
    }
    const auto [min_x, max_x, min_y, max_y] = minmax_coords.value();
    Value_t mid_x = (min_x + max_x) / 2.0;
    Value_t mid_y = (min_y + max_y) / 2.0;

    constexpr int parts{ 4 };
    Cities_t quadrants[parts];
    for (const City& coord : quadrant) {
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
        stack[depth].push_back(city.id);
        cities[city.id].layer = depth;
        ++added;
    }
    const int depth_next = depth + 1;
    int depth_max{depth};
    {
        auto [tmp_add, tmp_depth] = makeStackLayers(enn_tsp,
            quadrants[0], { min_x, mid_x, min_y, mid_y }, depth_next);
        added += tmp_add;
        depth_max = std::max(tmp_depth, depth_max);
    }
    {
        auto [tmp_add, tmp_depth] = makeStackLayers(enn_tsp,
            quadrants[1], { mid_x, max_x, min_y, mid_y }, depth_next);
        added += tmp_add;
        depth_max = std::max(tmp_depth, depth_max);
    }
    {
        auto [tmp_add, tmp_depth] = makeStackLayers(enn_tsp,
            quadrants[2], { min_x, mid_x, mid_y, max_y }, depth_next);
        added += tmp_add;
        depth_max = std::max(tmp_depth, depth_max);
    }
    {
        auto [tmp_add, tmp_depth] = makeStackLayers(enn_tsp,
            quadrants[3], { mid_x, max_x, mid_y, max_y }, depth_next);
        added += tmp_add;
        depth_max = std::max(tmp_depth, depth_max);
    }
    return std::make_pair(added, depth_max);
}

int createStack(DiscreteENN_TSP& enn_tsp, const MinMaxCoords& minmax_coords)
{
    enn_tsp.stack().clear();
    Cities_t cities_copy = enn_tsp.cities();
    const std::size_t num_cities = cities_copy.size();
    const auto& [cities_added, depth] =
        makeStackLayers(enn_tsp, cities_copy, minmax_coords, 1);
    if (static_cast<std::size_t>(cities_added) != num_cities) {
        const std::string& error_msg{ "number of cities in city_layers (" +
                                      std::to_string(cities_added) +
                                      ") doesn't match cities size (" +
                                      std::to_string(num_cities) +
                                      ").\nExiting." };
        utils::printErr(error_msg, "createStack");
        return -1;
    }
    return depth;
}

void makeGridInfo(DiscreteENN_TSP& enn_tsp, int depth, const MinMaxCoords& minmax_coords)
{
    Grid_t& grid{ enn_tsp.grid() };
    GridEdge_t& grid_edges{ enn_tsp.gridEdges() };
    const auto [min_x, max_x, min_y, max_y] = minmax_coords.value();
    int& grid_size{ enn_tsp.gridSize() };
    grid_size = 0x01 << depth;
    [[maybe_unused]] const int num_cells = grid_size * grid_size;
    grid.clear();
    grid = Grid_t(grid_size, std::vector<Indices_t>(grid_size, Indices_t{}));
    grid_edges.clear();
    grid_edges = GridEdge_t(grid_size, std::vector<Edges_t>(grid_size, Edges_t{}));
    const Value_t step_sizeX{ std::abs(max_x - min_x)/grid_size };
    const Value_t step_sizeY{ std::abs(max_y - min_y)/grid_size };
    enn_tsp.gridStart() = std::make_pair(min_x, min_y);
    enn_tsp.gridSteps() = std::make_pair(step_sizeX, step_sizeY);
    enn_tsp.gridStepMinMax() = std::make_pair(std::min(step_sizeX, step_sizeY), std::max(step_sizeX, step_sizeY));
    Cities_t& cities{ enn_tsp.cities() };
    for (City& city: cities) {
        const int x_cell = std::abs(city.x - min_x)/step_sizeX;
        const int y_cell = std::abs(city.y - min_y)/step_sizeY;
        city.x_cell = x_cell == grid_size ? grid_size - 1 : x_cell;
        city.y_cell = y_cell == grid_size ? grid_size - 1 : y_cell;
        grid[city.x_cell][city.y_cell].push_back(city.id);
        city.on_stack = true;
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

void fitPointsInWindow(std::vector<std::pair<sf::Vector2f, sf::Vector2f>>& lines, const sf::Vector2u& window_size,
                       const MinMaxCoords& minmax_coords, Value_t margin)
{
    if (lines.empty())
        return; // Check if the vector is not empty

    const auto [min_x, max_x, min_y, max_y] = minmax_coords.value();

    Value_t width{ max_x - min_x };
    Value_t height{ max_y - min_y };
    // Calculate scale factors for x and y to fit the plot within the window, considering margins
    Value_t scaleX = (window_size.x - 2 * margin) / width;
    Value_t scaleY = (window_size.y - 2 * margin) / height;

    // Apply scale and translation to the points
    for (auto& line : lines) {
        auto& [point1, point2] = line;
        point1.x = (point1.x - min_x) * scaleX + margin;
        point1.y = (point1.y - min_y) * scaleY + margin;
        point2.x = (point2.x - min_x) * scaleX + margin;
        point2.y = (point2.y - min_y) * scaleY + margin;
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

void drawThickLineSegment(sf::RenderWindow& window, const sf::Vector2f& start, const sf::Vector2f& end, float thickness, const sf::Color& color) {
    sf::Vector2f direction = end - start;
    float length = std::sqrt(direction.x * direction.x + direction.y * direction.y);
    direction /= length;  // Normalize

    // Perpendicular vector for thickness
    sf::Vector2f perpendicular(-direction.y * thickness / 2, direction.x * thickness / 2);

    // Four points of the thick line rectangle
    sf::VertexArray quad(sf::Quads, 4);
    quad[0].position = start + perpendicular;
    quad[1].position = start - perpendicular;
    quad[2].position = end - perpendicular;
    quad[3].position = end + perpendicular;

    // Set colors
    for (int i = 0; i < 4; ++i) {
        quad[i].color = color;
    }

    window.draw(quad);
}

void drawThickDashedLine(sf::RenderWindow& window, const sf::Vector2f& start, const sf::Vector2f& end, float dashLength, float spaceLength, float thickness, const sf::Color& color) {
    sf::Vector2f direction = end - start;
    float lineLength = std::sqrt(direction.x * direction.x + direction.y * direction.y);
    direction /= lineLength;  // Normalize

    sf::Vector2f currentPos = start;
    bool drawDash = true;
    while (lineLength > 0) {
        float segmentLength = drawDash ? dashLength : spaceLength;
        if (segmentLength > lineLength) segmentLength = lineLength;
        sf::Vector2f segmentEnd = currentPos + direction * segmentLength;

        if (drawDash) {
            drawThickLineSegment(window, currentPos, segmentEnd, thickness, color);
        }

        currentPos = segmentEnd;
        lineLength -= segmentLength;
        drawDash = !drawDash;  // Alternate between dash and space
    }
}

void drawPath(const DiscreteENN_TSP& enn_tsp, bool draw_coords,
              const std::string& title, float close_time, int highlight, int pos_start, int pos_end)
{
    const Indices_t& path{ enn_tsp.path() };
    const Cities_t& cities{enn_tsp.cities()};
    // Get the current desktop video mode
    const Value_t margin_size = 50;
    [[maybe_unused]] sf::VideoMode desktopMode =
        sf::VideoMode::getDesktopMode();
    // Create a full-screen window using the current desktop resolution
    MinMaxCoords minmax_coords;
    minmax_coords.update(cities);
    [[maybe_unused]] const auto [min_coord, max_coord] = minmax_coords.minmax();
    const auto [min_x, max_x, min_y, max_y] = minmax_coords.value();
    // min_x = min_y = static_cast<Value_t>(0);
    min_x = min_y = min_coord;
    max_x = max_y = max_coord;
    // sf::RenderWindow window(sf::VideoMode(minmax_coords.max_x,
    //                                       minmax_coords.max_y),
    //                         "TSP Route. Press Any Key to Close");
    // sf::RenderWindow window(
    // sf::VideoMode(minmax_coords.max_x, minmax_coords.max_y),
    // "TSP Route. Press Any Key to Close", sf::Style::Fullscreen);
    // sf::RenderWindow window(sf::VideoMode(800, 600), "TSP Route. Press Any Key to Close", sf::Style::Fullscreen);
    sf::RenderWindow window(desktopMode, "TSP Route. Press Any Key to Close",
                            sf::Style::Fullscreen);
    // sf::RenderWindow window(sf::VideoMode(max_coord + margin_size, max_coord + margin_size),
    //                         "TSP Route. Press Any Key to Close");
    // sf::RenderWindow window(sf::VideoMode(max_coord + margin_size, max_coord + margin_size),
    //                         "TSP Route. Press Any Key to Close", sf::Style::Fullscreen);

    // Main loop
    VectSF_t points;
    points.reserve(path.size());
    path2SFVector(path, cities, points);
    fitPointsInWindow(points, window.getSize(), minmax_coords, margin_size);
    Cities_t cities_copy = cities;
    fitPointsInWindow(cities_copy, window.getSize(), minmax_coords, margin_size);
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
    auto mid_x{(max_x + min_x)/2}, mid_y{(max_y + min_y)/2};
    // Define a vector of pairs to hold point coordinates
    std::vector<std::pair<sf::Vector2f, sf::Vector2f>> line_points = {
        {{min_x, mid_y}, {max_x, mid_y}},
        {{mid_x, min_y}, {mid_x, max_y}},
        {{mid_x/2, min_y}, {mid_x/2, max_y}},
        {{mid_x - mid_x/2, min_y}, {mid_x - mid_x/2, max_y}},
        {{mid_x + mid_x/2, min_y}, {mid_x + mid_x/2, max_y}},
        {{min_x, mid_y - mid_y/2}, {max_x, mid_y - mid_y/2}},
        {{min_x, mid_y + mid_y/2}, {max_x, mid_y + mid_y/2}},
    };
    const int layers = std::pow(2, 6);
    const Value_t step_x = (max_x - min_x)/layers;
    const Value_t step_y = (max_y - min_y)/layers;
    for (int idx{0}; idx != layers; ++idx) {
        line_points.push_back({{min_x + (idx*step_x), min_y}, {min_x + (idx*step_x), max_y}});
        line_points.push_back({{min_x, min_y + (idx*step_y)}, {max_x, min_y + (idx*step_y)}});
    }
    fitPointsInWindow(line_points, window.getSize(), minmax_coords, margin_size);
    std::vector<std::pair<sf::Vector2f, sf::Vector2f>> line_points2;
    const auto [grid_step_min, grid_step_max] = enn_tsp.gridStepMinMax();
    if (pos_start != -1 and pos_end != -1) {
        const City& city_start{ cities[pos_start] };
        const City& city_end{ cities[pos_end] };
        const Value_t p1_x = city_start.x;
        const Value_t p1_y = city_start.y;
        const Value_t p2_x = city_end.x;
        const Value_t p2_y = city_end.y;
        const Value_t width_x = p2_x - p1_x;
        const Value_t width_y = p2_y - p1_y;
        const Value_t length = std::sqrt(width_x*width_x + width_y*width_y) ;
        Value_t step_startX = p1_x;
        Value_t step_startY = p1_y;
        utils::printInfoFmt("{step_startX, step_startY} : {%f, %f}", "drawPath", step_startX, step_startY);
        utils::printInfoFmt("{step_endY, step_endY} : {%f, %f}", "drawPath", p2_x, p2_y);
        const int multiplierX = (p1_x < p2_x) ? 1 : -1;
        const int multiplierY = (p1_y < p2_y) ? 1 : -1;
        while ((multiplierX*step_startX < multiplierX*p2_x or utils::isEqual(step_startX, p2_x)) and (multiplierY*step_startY < multiplierY*p2_y or utils::isEqual(step_startY, p2_y))) {
            const Value_t step_startX_prev{ step_startX };
            const Value_t step_startY_prev{ step_startY };
            step_startX += utils::getRound2((width_x/length)*grid_step_min/2.);
            step_startY += utils::getRound2((width_y/length)*grid_step_min/2.);
            utils::printInfoFmt("{step_startX_prev, step_startY_prev} and {step_startX, step_startY} : {%f, %f} and {%f, %f}", "drawPath", step_startX_prev, step_startY_prev, step_startX, step_startY);
            line_points2.push_back({{step_startX_prev, step_startY_prev}, {step_startX, step_startY}});
        }
        utils::printInfoFmt("Drawing (pos_start, pos_end) : (%i, %i)", "drawPath", pos_start, pos_end);
    }
    fitPointsInWindow(line_points2, window.getSize(), minmax_coords, margin_size);

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
        for (const auto& line : line_points) {
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
            const auto color = city.on_stack ? sf::Color::Red : sf::Color::Blue;
            if (draw_coords or (city.id == highlight)) {
                sf::Text label;
                label.setFont(font); // Set the font
                label.setString(std::to_string(
                    city.id)); // Set the point's ID as the label text
                label.setCharacterSize(24); // Set the character size
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
            marker.setFillColor(color); // Use a contrasting color
            if (city.id == highlight) {
                marker.setFillColor(sf::Color::Magenta);
            }
            marker.setPosition(
                city.x,
                city.y); // Adjust position to center the marker on the vertex

            window.draw(marker);
        }
        window.draw(polygon); // Draw the polygon

        if (line_points2.size() != 0) {
            for (const auto& line : line_points2) {
                sf::VertexArray lineSegment(sf::Lines, 2);
                lineSegment[0].position = line.first;
                lineSegment[1].position = line.second;
                lineSegment[0].color = sf::Color::Black;
                lineSegment[1].color = sf::Color::Black;

                window.draw(lineSegment);
                drawThickDashedLine(window, line.first, line.second, 20, 10, 5, sf::Color::Black);
            }
        }

        window.display(); // Display what was drawn
    }
}

int images_count{1};
void savePath(const DiscreteENN_TSP& enn_tsp, bool draw_coords,
              const std::string& title, int highlight, int pos_start, int pos_end)
{
    const Indices_t& path{ enn_tsp.path() };
    const Cities_t& cities{enn_tsp.cities()};
    // Get the current desktop video mode
    const Value_t margin_size = 50;
    [[maybe_unused]] sf::VideoMode desktopMode =
        sf::VideoMode::getDesktopMode();
    // Create a full-screen window using the current desktop resolution
    MinMaxCoords minmax_coords;
    minmax_coords.update(cities);
    [[maybe_unused]] const auto [min_coord, max_coord] = minmax_coords.minmax();
    const auto [min_x, max_x, min_y, max_y] = minmax_coords.value();
    // min_x = min_y = static_cast<Value_t>(0);
    // min_x = min_y = min_coord;
    // max_x = max_y = max_coord;

    sf::RenderWindow window(sf::VideoMode(max_coord + margin_size, max_coord + margin_size),
                            "TSP Route. Press Any Key to Close", sf::Style::None);

    // Main loop
    VectSF_t points;
    points.reserve(path.size());
    path2SFVector(path, cities, points);
    // fitPointsInWindow(points, window.getSize(), minmax_coords, margin_size);
    Cities_t cities_copy = cities;
    // fitPointsInWindow(cities_copy, window.getSize(), minmax_coords, margin_size);
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
    auto mid_x{(max_x + min_x)/2}, mid_y{(max_y + min_y)/2};
    // Define a vector of pairs to hold point coordinates
    std::vector<std::pair<sf::Vector2f, sf::Vector2f>> line_points = {
        {{min_x, mid_y}, {max_x, mid_y}},
        {{mid_x, min_y}, {mid_x, max_y}},
        {{mid_x/2, min_y}, {mid_x/2, max_y}},
        {{mid_x - mid_x/2, min_y}, {mid_x - mid_x/2, max_y}},
        {{mid_x + mid_x/2, min_y}, {mid_x + mid_x/2, max_y}},
        {{min_x, mid_y - mid_y/2}, {max_x, mid_y - mid_y/2}},
        {{min_x, mid_y + mid_y/2}, {max_x, mid_y + mid_y/2}},
    };
    // fitPointsInWindow(line_points, window.getSize(), minmax_coords, margin_size);
    const int layers = std::pow(2, 6);
    const Value_t step_x = (max_x - min_x)/layers;
    const Value_t step_y = (max_y - min_y)/layers;
    for (int idx{0}; idx != layers; ++idx) {
        line_points.push_back({{min_x + (idx*step_x), min_y}, {min_x + (idx*step_x), max_y}});
        line_points.push_back({{min_x, min_y + (idx*step_y)}, {max_x, min_y + (idx*step_y)}});
    }

    std::vector<std::pair<sf::Vector2f, sf::Vector2f>> line_points2;
    const auto [grid_step_min, grid_step_max] = enn_tsp.gridStepMinMax();
    if (pos_start != -1 and pos_end != -1) {
        const City& city_start{ cities[pos_start] };
        const City& city_end{ cities[pos_end] };
        const Value_t p1_x = city_start.x;
        const Value_t p1_y = city_start.y;
        const Value_t p2_x = city_end.x;
        const Value_t p2_y = city_end.y;
        const Value_t width_x = p2_x - p1_x;
        const Value_t width_y = p2_y - p1_y;
        const Value_t length = std::sqrt(width_x*width_x + width_y*width_y) ;
        Value_t step_startX = p1_x;
        Value_t step_startY = p1_y;
        utils::printInfoFmt("{step_startX, step_startY} : {%f, %f}", "savePath", step_startX, step_startY);
        utils::printInfoFmt("{step_endY, step_endY} : {%f, %f}", "savePath", p2_x, p2_y);
        while ((step_startX < p2_x or utils::isEqual(step_startX, p2_x)) and (step_startY < p2_y or utils::isEqual(step_startY, p2_y))) {
            const Value_t step_startX_prev{ step_startX };
            const Value_t step_startY_prev{ step_startY };
            step_startX += utils::getRound2((width_x/length)*grid_step_min/2.);
            step_startY += utils::getRound2((width_y/length)*grid_step_min/2.);
            utils::printInfoFmt("{step_startX_prev, step_startY_prev} and {step_startX, step_startY} : {%f, %f} and {%f, %f}", "savePath", step_startX_prev, step_startY_prev, step_startX, step_startY);
            line_points2.push_back({{step_startX_prev, step_startY_prev}, {step_startX, step_startY}});
        }
        utils::printInfoFmt("Drawing (pos_start, pos_end) : (%i, %i)", "savePath", pos_start, pos_end);
    }

    sf::Clock clock;
    sf::Font font;
    font.loadFromFile("/usr/share/fonts/TTF/Roboto-Regular.ttf");
    const bool save_image{true};
    const std::string& image_filename{ "./Output/Images/" + title + "/DrawPlane%d.png" };
    // int image_count{1};
    if (save_image) {
        const auto image_dir = utils::getCleanPath(utils::stdfs::path(image_filename).parent_path());
        if (not utils::stdfs::is_directory(image_dir)) {
            utils::stdfs::create_directories(image_dir);
        }
    }
    // Activate the window to allow off-screen rendering
    window.setActive(true);
    window.clear(sf::Color::White); // Clear the window
    window.setTitle(title);

    // Draw each line segment
    int change_color{0};
    for (const auto& line : line_points) {
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
        const auto color = city.on_stack ? sf::Color::Red : sf::Color::Blue;
        if (draw_coords or (city.id == highlight)) {
            sf::Text label;
            label.setFont(font); // Set the font
            label.setString(std::to_string(
                city.id)); // Set the point's ID as the label text
            label.setCharacterSize(42); // Set the character size
            label.setFillColor(sf::Color::Black); // Set the text color
            label.setPosition(city.x + 10,
                              city.y -
                                  10); // Position the label near the marker
            window.draw(label); // Draw the label
        }

        sf::CircleShape marker(5); // Small circle with radius 5 pixels
        marker.setFillColor(color); // Use a contrasting color
        if (city.id == highlight) {
            marker.setFillColor(sf::Color::Magenta);
        }
        marker.setPosition(
            city.x,
            city.y); // Adjust position to center the marker on the vertex
        // utils::printInfoFmt("The coords are (id : %i, %f, %f)", "savePath", city.id, city.x, city.y);

        window.draw(marker);
    }
    window.draw(polygon); // Draw the polygon

    if (line_points2.size() != 0) {
        for (const auto& line : line_points2) {
            sf::VertexArray lineSegment(sf::Lines, 2);
            lineSegment[0].position = line.first;
            lineSegment[1].position = line.second;
            lineSegment[0].color = sf::Color::Black;
            lineSegment[1].color = sf::Color::Black;

            window.draw(lineSegment);
            drawThickDashedLine(window, line.first, line.second, 20, 10, 5, sf::Color::Black);
        }
    }

    sf::Texture texture;
    texture.create(window.getSize().x, window.getSize().y);
    texture.update(window);
    sf::Image screenshot = texture.copyToImage();

    if (!screenshot.saveToFile(utils::stringFmt(image_filename, images_count))) {
        std::cerr << "Failed to save screenshot." << std::endl;
    }
    ++images_count;
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
