#include <iostream>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <set>

struct Point
{
    int x;
    int y;
};

using DistanceMatrix = std::vector<std::vector<int>>;

int euclideanDistance(const Point& p1, const Point& p2)
{
    return std::round(std::sqrt(std::pow(p1.x - p2.x, 2) + std::pow(p1.y - p2.y, 2)));
}

class Instance
{
public:

    void load(std::filesystem::path path)
    {
        std::ifstream file{ path };
        if (!file.is_open())
        {
            std::cerr << "Error opening file: " << path << "\n";
        }

        std::string line;

        // skip to positions
        while (std::getline(file, line))
        {
            if (line.find("NODE_COORD_SECTION") != std::string::npos)
            {
                break;
            }
        }

        // read positions
        int index, x, y;
        while (std::getline(file, line))
        {
            if (line == "EOF")
            {
                break;
            }
            std::istringstream lineStream{ line };

            lineStream >> index >> x >> y;
            points.push_back({ x,y });
        }

        // Calculate distances
        size_t dim = points.size();
        M = DistanceMatrix(dim, std::vector<int>(dim));
        for (int i = 0; i < dim; ++i)
        {
            for (int j = 0; j < dim; ++j)
            {
                M[i][j] = euclideanDistance(points[i], points[j]);
            }
        }
    }

    std::vector<Point> points;
    DistanceMatrix M;
};

class Solution
{
public:
    void dump(std::filesystem::path path)
    {
        std::ofstream file{ path };
        if (!file.is_open())
        {
            std::cerr << "Error opening file: " << path << "\n";
        }

        for (const auto& point : path1)
        {
            file << point.x << ' ' << point.y << '\n';
        }
        file << '\n';
        for (const auto& point : path2)
        {
            file << point.x << ' ' << point.y << '\n';
        }
    }

    std::vector<Point> path1;
    std::vector<Point> path2;
};

class TSPSolver
{
public:
    virtual Solution run(const Instance& instance) = 0;
};

class GreedyNN : public TSPSolver
{
public:
    Solution run(const Instance& instance)
    {
        Solution sol;
        const auto& M = instance.M;
        const auto& points = instance.points;
        unsigned int dim = points.size();

        std::set<int> unVisited;
        for (int i = 0; i < dim; ++i)
        {
            unVisited.insert(i);
        }

        std::vector<int> path1;
        std::vector<int> path2;

        int start1Index = 0;
        int start2Index = 1;
        path1.push_back(start1Index);
        path2.push_back(start2Index);
        unVisited.erase(start1Index);
        unVisited.erase(start2Index);

        while (!unVisited.empty())
        {
            int minDistance = std::numeric_limits<int>::max();
            int bestPoint = -1;
            for (int pointIndex : unVisited)
            {
                int distance = M[path1[path1.size() - 1]][pointIndex];
                if (distance < minDistance)
                {
                    minDistance = distance;
                    bestPoint = pointIndex;
                }
            }
            unVisited.erase(bestPoint);
            path1.push_back(bestPoint);

            minDistance = std::numeric_limits<int>::max();
            bestPoint = -1;
            for (int pointIndex : unVisited)
            {
                int distance = M[path2[path2.size() - 1]][pointIndex];
                if (distance < minDistance)
                {
                    minDistance = distance;
                    bestPoint = pointIndex;
                }
            }
            unVisited.erase(bestPoint);
            path2.push_back(bestPoint);
        }

        for (int pointIndex : path1)
        {
            sol.path1.push_back(points[pointIndex]);
        }
        for (int pointIndex : path2)
        {
            sol.path2.push_back(points[pointIndex]);
        }

        return sol;
    }

private:
    Point getNearest(const Point& point, const DistanceMatrix& M)
    {

    }
};

class GreedyCycle : public TSPSolver
{
public:
    Solution run(const Instance& instance)
    {
        // TODO:
        return Solution{};
    }
};

int main(int argc, char* argv[])
{
    std::filesystem::path workDir = argv[1];

    Instance instance;
    instance.load(workDir / "kroA100.tsp.txt");

    auto greadyNNSolver = std::make_unique<GreedyNN>();
    Solution solution = greadyNNSolver->run(instance);
    solution.dump(workDir / "kroA100-greadyNN-solution.txt");

    auto greadyCycleSolver = std::make_unique<GreedyCycle>();
    solution = greadyCycleSolver->run(instance);
    solution.dump(workDir / "kroA100-greadyCycle-solution.txt");

    return 0;
}