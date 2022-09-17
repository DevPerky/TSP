#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <iostream>
#include "tsp.hpp"
#include <random>
#include <thread>
#include <string>

struct Params {
    TSP::Params tspParams;
    int generations;
    size_t numPoints;
    int seed;

    static Params fromArgs(int argc, const char **argv) {
        Params params;
        params.generations = 100;
        params.numPoints = 10;
        params.seed = 0;
        params.tspParams.killFactor = 0.8f;
        params.tspParams.mixFactor = 0.2f;
        params.tspParams.mutationFactor = 0.2f;
        params.tspParams.populationSize = 100;

        for(int i = 0; i < argc - 1; i++) {
            std::string parameter(argv[i]);
            std::string value(argv[i + 1]);
            try {
                if(parameter == "--population") {
                    params.tspParams.populationSize = std::stoi(value);
                }
                else if(parameter == "--mutation") {
                    params.tspParams.mutationFactor = std::stof(value);
                }
                else if(parameter == "--kill") {
                    params.tspParams.killFactor = std::stof(value);
                }
                else if(parameter == "--mix") {
                    params.tspParams.mixFactor = std::stof(value);
                }
                else if(parameter == "--generations") {
                    params.generations = std::stoi(value);
                }
                else if(parameter == "--points") {
                    params.numPoints = std::stoi(value);
                }
                else if(parameter == "--seed") {
                    params.seed = std::stoi(value);
                }
            }
            catch (std::exception e) {
                std::string errorMessage = "Failed parsing parameter " + parameter + "\n";
                throw errorMessage;
            }
        }

        return params;

    }

    const std::string toString() const {
        return "generations:\t" + std::to_string(generations) + "\n" +
            "points:\t\t" + std::to_string(numPoints) + "\n" +
            "seed:\t\t" + std::to_string(seed) + "\n" +
            "population:\t" + std::to_string(tspParams.populationSize) + "\n" +
            "mutation:\t" + std::to_string(tspParams.mutationFactor) + "\n" + 
            "kill:\t\t" + std::to_string(tspParams.killFactor) + "\n" +
            "mix:\t\t" + std::to_string(tspParams.mixFactor) + "\n";
    }
};

void randomizePoints(std::vector<glm::vec2> &points, const int xMax, const int yMax) {
    for(auto &p : points) {
        p.x = rand() % xMax;
        p.y = rand() % yMax;
    }
}

std::string generateWindowTitleString(float timeSeconds, int generations, float bestSolution, float worstSolution) {
    return std::string("TLP")
    + " | "
    + "Time: " + std::to_string(timeSeconds) + " seconds"
    + " | " 
    + "generation: " + std::to_string(generations)
    + " | "
    + "best: " + std::to_string(bestSolution)
    + " | "
    + "worst: " + std::to_string(worstSolution);
}

GLFWwindow *createWindow(int width, int height, const char *title) {
    GLFWwindow *window = glfwCreateWindow(width, height, title, NULL, NULL);
    glfwSetWindowAttrib(window, GLFW_RESIZABLE, GLFW_FALSE);
    return window;
}

void drawTriangles(const std::vector<glm::vec2> &vertices) {
    assert(vertices.size() % 3 == 0);
    
    glBegin(GL_TRIANGLES);
    
    for(const auto &v : vertices) {
        glVertex2f(v.x, v.y);
    }

    glEnd();
}

void drawLine(const glm::vec2 &from, const glm::vec2 &to, float width) {
    auto lineDirection = glm::normalize(to - from);
    auto perpendicular = glm::vec2(-lineDirection.y, lineDirection.x);

    const std::vector<glm::vec2> vertices = {
        from + perpendicular * (width * 0.5f),
        to + perpendicular * (width * 0.5f),
        to - perpendicular * (width * 0.5f),

        from + perpendicular * (width * 0.5f),
        to - perpendicular * (width * 0.5f),
        from - perpendicular * (width * 0.5f)
    };

    drawTriangles(vertices);
}

void drawRectangle(float left, float top, float right, float bottom) {
    const std::vector<glm::vec2> vertices = {
        glm::vec2(right, top),
        glm::vec2(left, top),
        glm::vec2(left, bottom),
        
        glm::vec2(right, top),
        glm::vec2(left, bottom),
        glm::vec2(right, bottom)
    };

    drawTriangles(vertices);
}

void drawSolution(const TSP::Solution &solution, const std::vector<glm::vec2> &points, float lineWidth) {
    assert(points.size() == solution.indices.size());

    for(size_t i = 0; i < solution.indices.size(); i++) {
            size_t indexFrom = i;
            size_t indexTo = (i + 1) % solution.indices.size();
            drawLine(points[solution.indices[i]], points[solution.indices[indexTo]], lineWidth);
    }
}

int main(int argc, const char **argv) {


    Params params;
    try {
        params = Params::fromArgs(argc, argv);
    }
    catch (std::exception e) {
        std::cout << "Usage error: " << e.what() << "\n";
        std::cout << "type anything to exit" << std::endl;
        std::cin.get();
        exit(EXIT_FAILURE);
    }

    std::cout << "Parameters:\n" << params.toString() << std::endl;
    bool displayWorst = false;

    std::default_random_engine randomEngine(params.seed);
    srand(params.seed);
    
    const int windowWidth = 1280;
    const int windowHeight = 720;
    auto projectionMatrix = glm::ortho((float)0, (float)windowWidth, (float)windowHeight, (float)0);

    glfwInit();
    auto window = createWindow(windowWidth, windowHeight, "Traveling salesman problem");
    glfwMakeContextCurrent(window);
    
    glMatrixMode(GL_PROJECTION_MATRIX);
    glLoadMatrixf(glm::value_ptr(projectionMatrix));

    std::vector<glm::vec2> points;
    points.resize(params.numPoints);
    randomizePoints(points, windowWidth, windowHeight);
    
    TSP tsp(params.tspParams, points, randomEngine);
    float timeStart = 0;
    float lastTime = 0;

    int currentGeneration = 0;
    while(!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        glClear(GL_COLOR_BUFFER_BIT);
    
        if(currentGeneration < params.generations) {
            
            if(timeStart == 0) {
                timeStart = glfwGetTime();
            }

            tsp.generateOneGeneration();
            lastTime = glfwGetTime();

            std::string windowTitle = generateWindowTitleString(
                lastTime - timeStart,
                currentGeneration + 1,
                tsp.measureSolutionLength(tsp.getCurrentBestSolution()),
                tsp.measureSolutionLength(tsp.getCurrentWorstSolution())
            );

            currentGeneration++;
            glfwSetWindowTitle(window, windowTitle.c_str());
        }
        else {
            std::this_thread::sleep_for(std::chrono::milliseconds(33));
        }
        
        const auto &solution = tsp.getCurrentBestSolution();
        const auto &worstSolution = tsp.getCurrentWorstSolution();
        
        if(displayWorst) {
            glColor4f(1, 0, 0, 0.1f);
            drawSolution(worstSolution, points, 1 );
        }
        
        glColor3f(0, 1, 0);
        drawSolution(solution, points, 5);
        
        glColor3f(1, 1, 0);
        for(auto &p: points) {
            drawRectangle(p.x - 5, p.y - 5, p.x + 5, p.y + 5);
        }

        glfwSwapBuffers(window);
    }

    std::cout << "Results: " << std::endl;
    std::cout << "Solution length:\t" << tsp.measureSolutionLength(tsp.getCurrentBestSolution()) << std::endl;
    std::cout << "Time taken:\t\t" << (lastTime - timeStart) << " seconds" << std::endl;

    
}
