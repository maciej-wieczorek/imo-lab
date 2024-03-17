#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <vector>
#include <math.h>
#include <filesystem>
#include <random>
#include <chrono>
#include <array>

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

int getRandomNumber(int min, int max) {
    static unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    static std::mt19937 generator(seed);
    std::uniform_int_distribution<int> distribution(min, max);
    return distribution(generator);
}

int getFurthestindex(int index, const DistanceMatrix& M)
{
	int maxdistance = 0;
	int maxindex = -1;
	int dim = M[index].size();

	for (int i = 0; i < dim; ++i)
	{
		if (M[index][i] > maxdistance)
		{
			maxdistance = M[index][i];
			maxindex = i;
		}
	}
	return maxindex;
}

int getNearestIndex(int index, const DistanceMatrix& M)
{
    int minDinstance = std::numeric_limits<int>::max();
	int minIndex = -1;
	int dim = M[index].size();

	for (int i = 0; i < dim; ++i)
	{
        if (i != index)
        {
            if (M[index][i] < minDinstance)
            {
                minDinstance = M[index][i];
                minIndex = i;
            }
        }
	}
	return minIndex;
}

int getLenDiff(int placementIndex, int pointIndex, const std::vector<int>& path, const DistanceMatrix& M)
{
    int l = path[placementIndex];
    int r = placementIndex == path.size() - 1 ? 0 : path[placementIndex + 1];
    return M[l][pointIndex] + M[pointIndex][r] - M[l][r];
}

int getSecondMin(const std::vector<int>& vec) {
    if (vec.size() < 2)
        throw std::invalid_argument("Vector size is less than 2");

    int min1 = std::numeric_limits<int>::max();
    int min2 = std::numeric_limits<int>::max();

    for (int num : vec) {
        if (num < min1) {
            min2 = min1;
            min1 = num;
        }
        else if (num < min2 && num != min1) {
            min2 = num;
        }
    }

    if (min2 == std::numeric_limits<int>::max())
    {
        return min1;
    }

    return min2;
}

std::pair<int,int> getRegret(int pointIndex, const std::vector<int>& path, const DistanceMatrix& M)
{

    std::vector<int> lenDiffs;
    int bestPlacement = -1;
    int bestLenDiff = std::numeric_limits<int>::max();
    for (int i = 0; i < path.size(); ++i)
    {
        int lenDiff = getLenDiff(i, pointIndex, path, M);
        if (lenDiff < bestLenDiff)
        {
            bestLenDiff = lenDiff;
            bestPlacement = i;
        }
        lenDiffs.push_back(lenDiff);
    }

    int secondBestLenDiff = getSecondMin(lenDiffs);

    return std::make_pair(bestPlacement, secondBestLenDiff - bestLenDiff);
}

double calculateMean(const std::vector<int>& vec)
{
    if (vec.empty())
    {
        return 0.0;
    }

    double sum = 0.0;
    for (const auto& element : vec) {
        sum += element;
    }

    return sum / vec.size();
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

        name = path.stem().string();

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
    std::string name;
    int startIndex;
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

    int getScore()
    {
        int score = 0;
        for (int i = 0; i < path1.size() - 1; ++i)
        {
            score += euclideanDistance(path1[i], path1[i + 1]);
        }
        for (int i = 0; i < path2.size() - 1; ++i)
        {
            score += euclideanDistance(path1[i], path1[i + 1]);
        }

        return score;
    }

    std::vector<Point> path1;
    std::vector<Point> path2;
};

class TSPSolver
{
public:
    virtual const char* getName() = 0;
    virtual Solution run(const Instance& instance) = 0;
};

class GreedyNN : public TSPSolver
{
public:
    const char* getName() { return "GreedyNN"; }
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

        int start1Index = instance.startIndex;
        int start2Index = getFurthestindex(start1Index, M);

        std::array<std::vector<int>, 2> paths;
        paths[0].push_back(start1Index);
        paths[1].push_back(start2Index);
        unVisited.erase(start1Index);
        unVisited.erase(start2Index);

        while (!unVisited.empty())
        {
            for (auto& path : paths)
            {
				int minDistance = std::numeric_limits<int>::max();
				int bestPlacement = -1;
				int bestPoint = -1;

                for (int i = 0; i < path.size(); ++i)
                {
					for (int pointIndex : unVisited)
					{
						int distance = M[path[i]][pointIndex];
						if (distance < minDistance)
						{
							minDistance = distance;
							bestPoint = pointIndex;
                            bestPlacement = i;
						}
					}
                }

				unVisited.erase(bestPoint);
                path.insert(path.begin() + bestPlacement + 1, bestPoint);
            }
        }

        for (int pointIndex : paths[0])
        {
            sol.path1.push_back(points[pointIndex]);
        }
        for (int pointIndex : paths[1])
        {
            sol.path2.push_back(points[pointIndex]);
        }

        return sol;
    }
};

class GreedyCycle : public TSPSolver
{
public:
    const char* getName() { return "GreedyCycle"; }

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

        int start1Index = instance.startIndex % points.size();
        int start2Index = getFurthestindex(start1Index, M);

        std::array<std::vector<int>, 2> paths;
        paths[0].push_back(start1Index);
        paths[1].push_back(start2Index);
        unVisited.erase(start1Index);
        unVisited.erase(start2Index);

        while (!unVisited.empty())
        {
            for (auto& path : paths)
            {
				int minLenDiff = std::numeric_limits<int>::max();
				int bestPlacement = -1;
				int bestPoint = -1;

                for (int i = 0; i < path.size(); ++i)
                {
					for (int pointIndex : unVisited)
					{
                        int lenDiff = getLenDiff(i, pointIndex, path, M);
						if (lenDiff < minLenDiff)
						{
							minLenDiff = lenDiff;
							bestPoint = pointIndex;
                            bestPlacement = i;
						}
					}
                }

				unVisited.erase(bestPoint);
                path.insert(path.begin() + bestPlacement + 1, bestPoint);
            }
        }

        for (int pointIndex : paths[0])
        {
            sol.path1.push_back(points[pointIndex]);
        }
        for (int pointIndex : paths[1])
        {
            sol.path2.push_back(points[pointIndex]);
        }

        return sol;
    }
};

class GreedyRegret : public TSPSolver
{
public:
    const char* getName() { return "GreedyRegret"; }

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

        int start1Index = instance.startIndex;
        int start2Index = getFurthestindex(start1Index, M);

        int second1Index = getNearestIndex(start1Index, M);
        int second2Index = getNearestIndex(start2Index, M);

        std::array<std::vector<int>, 2> paths;
        paths[0].push_back(start1Index);
        paths[0].push_back(second1Index);
        paths[1].push_back(start2Index);
        paths[1].push_back(second2Index);
        unVisited.erase(start1Index);
        unVisited.erase(start2Index);
        unVisited.erase(second1Index);
        unVisited.erase(second2Index);

        while (!unVisited.empty())
        {
            for (auto& path : paths)
            {
                int maxRegret = -1;
				int bestPlacement = -1;
				int bestPoint = -1;

				for (int pointIndex : unVisited)
				{
                    auto regret = getRegret(pointIndex, path, M);
                    int placement = regret.first;
                    int regretValue = regret.second;
					if (regretValue > maxRegret)
					{
                        maxRegret = regretValue;
                        bestPlacement = placement;
                        bestPoint = pointIndex;
					}
				}

				unVisited.erase(bestPoint);
                path.insert(path.begin() + bestPlacement + 1, bestPoint);
            }
        }

        for (int pointIndex : paths[0])
        {
            sol.path1.push_back(points[pointIndex]);
        }
        for (int pointIndex : paths[1])
        {
            sol.path2.push_back(points[pointIndex]);
        }

        return sol;
    }
};

int main(int argc, char* argv[])
{
    std::filesystem::path workDir = argc > 1 ? argv[1] : "workdir";


    std::vector<Instance> instances;
    std::string instanceNames[] = { "kroA100.tsp.txt", "kroB100.tsp.txt" };
    for (const auto& instanceName : instanceNames)
    {
        instances.emplace_back();
        instances[instances.size() - 1].load(workDir / instanceName);
    }

    std::cout << "solver, instance, score\n";

    TSPSolver* solvers[] = { new GreedyNN, new GreedyCycle, new GreedyRegret };
    for (TSPSolver* solver : solvers)
    {
        for (auto& instance : instances)
        {
            std::vector<int> scores;
			int bestScore = std::numeric_limits<int>::max();
			Solution bestSolution;

            for (int i = 0; i < 100; ++i)
            {
                instance.startIndex = i;
                Solution solution = solver->run(instance);
                int score = solution.getScore();
                scores.push_back(score);
                if (score < bestScore);
                {
                    bestScore = score;
                    bestSolution = solution;
                }
            }

            int worstScore = *std::max_element(scores.begin(), scores.end());
            int avgScore = calculateMean(scores);

			std::string solFileName = instance.name + '-' + solver->getName() + "-solution.txt";

			bestSolution.dump(workDir / solFileName);

            std::cout << solver->getName() << ", " << instance.name << ", " << avgScore << " (" << bestScore << '-' << worstScore << ")\n";
        }
    }

    return 0;
}