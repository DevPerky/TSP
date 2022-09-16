#include <vector>
#include <glm/vec2.hpp>
#include <random>

class TSP {
public:
    struct Solution {
        std::vector<size_t> indices;
        Solution mixedWith(const Solution &other, float mixFactor) const;
        void mutate();
    };

    struct Params {
        size_t  populationSize;
        float   mutationFactor;
        float   killFactor;
        float   mixFactor;
    };

private:
    std::default_random_engine &m_randomEngine;
    const std::vector<glm::vec2> &m_points;
    const Params &m_params;
    std::vector<Solution> m_population;

    Solution generateRandomSolution() const;

public:

    TSP(const Params &params, const std::vector<glm::vec2> &points, std::default_random_engine &randomEngine);
    
    inline const Solution getCurrentBestSolution() const {
        return m_population.front();
    }

    inline const Solution getCurrentWorstSolution() const {
        return m_population.back();
    }

    void generateOneGeneration();

    const float measureSolutionLength(const Solution &solution) const;

};