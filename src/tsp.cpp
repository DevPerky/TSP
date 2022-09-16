#include "tsp.hpp"
#include <stack>
#include <algorithm>
#include <random>
#include <iostream>
#include <functional>
#include <glm/geometric.hpp>

TSP::Solution TSP::Solution::mixedWith(const Solution &other, float mixFactor) const{
    assert(mixFactor < 1 && mixFactor > 0);
    assert(other.indices.size() == indices.size());

    TSP::Solution result;
    result.indices.resize(indices.size());

    size_t sliceSize = (size_t)std::ceil(((float)indices.size() * mixFactor));
    size_t sliceStart = rand() % (indices.size() - sliceSize);

    std::vector<size_t> slice, unusedIndices;
    slice.resize(sliceSize);

    for(size_t i = 0; i < sliceSize; i++) {
        size_t sliceIndex = sliceStart + i;

        slice[i] = other.indices[sliceIndex];
        result.indices[sliceIndex] = slice[i];
    }

    for(const auto &index : indices) {
        if(std::find(slice.begin(), slice.end(), index) == slice.end()) {
            unusedIndices.push_back(index);
        }
    }

    for(size_t i = 0; i < unusedIndices.size(); i++) {
        size_t resultIndex = (sliceStart + sliceSize + i) % indices.size();
        result.indices[resultIndex] = unusedIndices[i];
    }

    return result;
}



void TSP::Solution::mutate() {
    int swapIndex1 = rand() % indices.size();
    int swapIndex2 = (swapIndex1 + (rand() % (indices.size() - 1)))  % indices.size();

    assert(swapIndex1 != swapIndex2);

    size_t valueAtIndex1 = indices[swapIndex1];
    indices[swapIndex1] = indices[swapIndex2];
    indices[swapIndex2] = valueAtIndex1;
}


TSP::TSP(const Params &params, const std::vector<glm::vec2> &points, std::default_random_engine &randomEngine)
: m_params(params), m_points(points), m_randomEngine(randomEngine) {
    
    // Generate a random initial population
    m_population.resize(m_params.populationSize);
    for(int i = 0; i < m_population.size(); i++) {
        m_population[i] = generateRandomSolution();
    }

}

TSP::Solution TSP::generateRandomSolution() const {

    Solution solution;
    solution.indices.resize(m_points.size());

    for(size_t i = 0; i < m_points.size(); i++) {
        solution.indices[i] = i;
    }

    std::shuffle(solution.indices.begin(), solution.indices.end(), m_randomEngine);

    return solution;
}

void TSP::generateOneGeneration() {

    std::sort(m_population.begin(), m_population.end(), [this](Solution &a, Solution &b) -> bool {
        return measureSolutionLength(a) < measureSolutionLength(b);
    });
        
    const Solution &bestSolution = m_population.front();
    const Solution &worstSolution = m_population.back();
    
    int numSuccessfulSolutions = m_population.size() * (1.0 - m_params.killFactor);
    
    for(size_t j = m_params.killFactor * m_population.size(); j < m_population.size(); j++) {
        
        size_t indexOfParentOne = (rand() % numSuccessfulSolutions); 
        size_t indexOfParentTwo = (rand() % numSuccessfulSolutions); 

        const Solution &parentOne = m_population[indexOfParentOne];
        const Solution &parentTwo = m_population[indexOfParentTwo];

        m_population[j] = parentOne.mixedWith(parentTwo, m_params.mixFactor);
    }

    for(int i = 1; i < numSuccessfulSolutions; i++) {
        float roll = (float)rand() / (float)RAND_MAX;
        bool shouldMutate = roll <= m_params.mutationFactor;

        if(!shouldMutate) {
            continue;
        }
        
        m_population[i].mutate();
    }

}

const float TSP::measureSolutionLength(const Solution &solution) const {
    assert(solution.indices.size() == m_points.size());

    float total = 0;
    
    for(size_t i = 0; i < solution.indices.size(); i++) {
        size_t indexFrom = solution.indices[i];
        size_t indexTo = solution.indices[(i + 1) % solution.indices.size()];

        const glm::vec2 &from = m_points[indexFrom];
        const glm::vec2 &to = m_points[indexTo];

        float distance = glm::length(to - from);
        total += distance;
    }

    return total;
}
