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
#include <algorithm>
#include <functional>

struct Point
{
    int x;
    int y;
};

using DistanceMatrix = std::vector<std::vector<int>>;
using Path = std::vector<int>;

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
    int r = placementIndex == path.size() - 1 ? path[0] : path[placementIndex + 1];
    return M[l][pointIndex] + M[pointIndex][r] - M[l][r];
}

int getSecondMin(const std::vector<int>& vec) {
    if (vec.size() < 2)
        return vec[0];

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
    int regretVal = secondBestLenDiff - bestLenDiff;

    return std::make_pair(bestPlacement, regretVal - 0.5 * bestLenDiff);
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
    void dump(std::filesystem::path path, const Instance& instance)
    {
        std::ofstream file{ path };
        if (!file.is_open())
        {
            std::cerr << "Error opening file: " << path << "\n";
        }

        Point point{};
        for (const auto& pointIndex : path1)
        {
            point = instance.points[pointIndex];
            file << point.x << ' ' << point.y << '\n';
        }
        file << '\n';
        for (const auto& pointIndex : path2)
        {
            point = instance.points[pointIndex];
            file << point.x << ' ' << point.y << '\n';
        }
    }

    int getScore(const Instance& instance)
    {
        if (path1.size() == 0)
        {
            return std::numeric_limits<int>::max();
        }

        const auto& M = instance.M;
        int score = 0;
        for (int i = 0; i < path1.size() - 1; ++i)
        {
            score += M[path1[i]][path1[i + 1]];
        }
        for (int i = 0; i < path2.size() - 1; ++i)
        {
            score += M[path2[i]][path2[i + 1]];
        }
		score += M[path1[path1.size()-1]][path1[0]];
		score += M[path2[path2.size()-1]][path2[0]];

        return score;
    }

    std::vector<int> path1;
    std::vector<int> path2;
    int score;
};

class TSPSolver
{
public:
    virtual ~TSPSolver() {}
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
                        if (i==0 || i == path.size()-1)
                        {
                            int distance = M[path[i]][pointIndex];
                            if (distance < minDistance)
                            {
                                minDistance = distance;
                                bestPoint = pointIndex;
                                bestPlacement = i;
                            }
                        }
                        else
                        {
                            int distance = M[path[i]][pointIndex] + M[path[i+1]][pointIndex] - M[path[i]][path[i+1]];
                            if (distance < minDistance)
                            {
                                minDistance = distance;
                                bestPoint = pointIndex;
                                bestPlacement = i;
                            }
                        }

					}
                }

				unVisited.erase(bestPoint);
                path.insert(path.begin() + bestPlacement + 1, bestPoint);
            }
        }

        sol.path1 = paths[0];
        sol.path2 = paths[1];

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

        sol.path1 = paths[0];
        sol.path2 = paths[1];

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

        std::array<std::vector<int>, 2> paths;
        paths[0].push_back(start1Index);
        paths[1].push_back(start2Index);
        unVisited.erase(start1Index);
        unVisited.erase(start2Index);

        while (!unVisited.empty())
        {
            for (auto& path : paths)
            {
                int maxRegret = std::numeric_limits<int>::min();
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

        sol.path1 = paths[0];
        sol.path2 = paths[1];

        return sol;
    }
};

class RandomSolver : public TSPSolver
{
public:
    const char* getName() { return "RandomSolver"; }

    Solution run(const Instance& instance)
    {
        Solution sol;

        auto points = instance.points;

        while (!points.empty())
        {
            int randomPointIndex = getRandomNumber(0, points.size() - 1);
            sol.path1.push_back(randomPointIndex);
            points.erase(points.begin() + randomPointIndex);

            randomPointIndex = getRandomNumber(0, points.size() - 1);
            sol.path2.push_back(randomPointIndex);
            points.erase(points.begin() + randomPointIndex);
        }

        return sol;
    }
};


class LocalSearch
{
public:
    virtual ~LocalSearch() {}
    virtual const char* getName() = 0;
    virtual Solution run(const Instance& instance, const Solution& initialSolution) = 0;
};

class RandomLocalSearch : public LocalSearch
{
public:
    const char* getName() { return "RandomLocalSearch"; }

    Solution run(const Instance& instance, const Solution& initialSolution)
    {
        Solution sol;
        const auto& M = instance.M;
        const auto& points = instance.points;
        unsigned int dim = points.size();

        return sol;
    }
};

struct ScoredMove
{
    int vertex1;
    int vertex2;
    int distanceDelta;
    int pathIndex;
};

ScoredMove interPathVertexSwap(const Instance& instance, const Solution& sol, bool greedy)
{
    const auto& M = instance.M;

    ScoredMove bestInterMove;
    bestInterMove.distanceDelta = 0;
    const int n1 = sol.path1.size();
    const int n2 = sol.path2.size();
    for (int i = 0; i < n1; ++i)
    {
        for (int j = 0; j < n2; ++j)
        {
            int v1 = sol.path1[i], v1Before = i == 0 ? sol.path1[n1 - 1] : sol.path1[i - 1], v1After = sol.path1[(i + 1) % n1];
            int v2 = sol.path2[j], v2Before = j == 0 ? sol.path2[n2 - 1] : sol.path2[j - 1], v2After = sol.path2[(j + 1) % n2];
            int distanceNow = M[v1][v1After] + M[v1][v1Before] + M[v2][v2Before] + M[v2][v2After];
            int distanceAfter = M[v1][v2After] + M[v1][v2Before] + M[v2][v1After] + M[v2][v1Before];

            int distanceDelta = distanceAfter - distanceNow;

            if (distanceDelta < bestInterMove.distanceDelta)
            {
                bestInterMove = ScoredMove{ i, j, distanceDelta, -1 };
                if (greedy)
                {
                    return bestInterMove;
                }
            }
        }
    }

    return bestInterMove;
}

ScoredMove intraPathVertexSwap(const Instance& instance, const Solution& sol, bool greedy)
{
    const auto& M = instance.M;

    ScoredMove bestIntraMove;
    bestIntraMove.distanceDelta = 0;
    int pathIndex = 0;
    for (auto* pathPtr : { &sol.path1, &sol.path2 })
    {
        auto& path = *pathPtr;
        const int n = path.size();
        for (int i = 1; i < n; ++i)
        {
            for (int j = i + 2; j < n; ++j)
            {
                int v1 = path[i], v1Before = i == 0 ? path[n - 1] : path[i - 1], v1After = path[(i + 1) % n];
                int v2 = path[j], v2Before = j == 0 ? path[n - 1] : path[j - 1], v2After = path[(j + 1) % n];
                int distanceNow = M[v1][v1After] + M[v1][v1Before] + M[v2][v2Before] + M[v2][v2After];
                int distanceAfter = M[v1][v2After] + M[v1][v2Before] + M[v2][v1After] + M[v2][v1Before];

                int distanceDelta = distanceAfter - distanceNow;

                if (distanceDelta < bestIntraMove.distanceDelta)
                {
                    bestIntraMove = ScoredMove{ i, j, distanceDelta, pathIndex };

                    if (greedy)
                    {
                        return bestIntraMove;
                    }
                }
            }
        }

        ++pathIndex;
    }

    return bestIntraMove;
}

ScoredMove edgeSwap(const Instance& instance, const Solution& sol, bool greedy)
{
    const auto& M = instance.M;

    ScoredMove bestMove;
    bestMove.distanceDelta = 0;
    int pathIndex = 0;
    for (auto* pathPtr : { &sol.path1, &sol.path2 })
    {
        auto& path = *pathPtr;
        const int n = path.size();
        for (int i = 0; i < n - 2; ++i)
        {
            for (int j = i + 1; j < n - 1; ++j)
            {
                int v1 = path[i];
                int v2 = path[j];

                int distanceDelta = M[v1][(v1 + 1) % n] - M[v2][(v2 + 1) % n] + M[v1][v2] + M[(v1 + 1) % n][(v2 + 1) % n];

                if (distanceDelta < bestMove.distanceDelta)
                {
                    bestMove = ScoredMove{ i, j, distanceDelta, pathIndex };

                    if (greedy)
                    {
                        return bestMove;
                    }
                }
            }
        }

        ++pathIndex;
    }

    return bestMove;
}

class GreedyVertexLocalSearch : public LocalSearch
{
public:
    const char* getName() { return "GreedyVertexLocalSearch"; }

    Solution run(const Instance& instance, const Solution& initialSolution)
    {
        Solution sol = initialSolution;
        sol.score = sol.getScore(instance);

        const auto& M = instance.M;

        while (true)
        {
            ScoredMove bestInterMove = interPathVertexSwap(instance, sol, true);
            ScoredMove bestIntraMove = intraPathVertexSwap(instance, sol, true);

            // Apply best move
            if (bestInterMove.distanceDelta < bestIntraMove.distanceDelta)
            {
                std::swap(sol.path1[bestInterMove.vertex1], sol.path2[bestInterMove.vertex2]);
                sol.score += bestInterMove.distanceDelta;
            }
            else if (bestIntraMove.distanceDelta < bestInterMove.distanceDelta)
            {
                Path* pathForBestIntraMove = &sol.path1;

                if (bestIntraMove.pathIndex == 1)
                {
                    pathForBestIntraMove = &sol.path2;
                }

                std::swap(pathForBestIntraMove->at(bestIntraMove.vertex1), pathForBestIntraMove->at(bestIntraMove.vertex2));
                sol.score += bestInterMove.distanceDelta;
            }
            else
            {
                break;
            }

        }

        return sol;
    }
};

class GreedyEdgeLocalSearch : public LocalSearch
{
public:
    const char* getName() { return "GreedyEdgeLocalSearch"; }

    Solution run(const Instance& instance, const Solution& initialSolution)
    {
        Solution sol = initialSolution;
        sol.score = sol.getScore(instance);

        const auto& M = instance.M;

        while (true)
        {
            ScoredMove bestMove = edgeSwap(instance, sol, true);
            // Apply best move
            if (bestMove.distanceDelta < 0)
            {
                Path* pathForBestMove = &sol.path1;
                if (bestMove.pathIndex == 1)
                {
                    pathForBestMove = &sol.path2;
                }

                std::reverse(pathForBestMove->begin() + bestMove.vertex1, pathForBestMove->begin() + bestMove.vertex2 + 1);
                sol.score += bestMove.distanceDelta;
            }
            else
            {
                break;
            }

        }

        return sol;
    }
};

class SteepVertexLocalSearch : public LocalSearch
{
public:
    const char* getName() { return "SteepVertexLocalSearch"; }

    Solution run(const Instance& instance, const Solution& initialSolution)
    {
        Solution sol = initialSolution;
        sol.score = sol.getScore(instance);

        const auto& M = instance.M;

        while (true)
        {

            ScoredMove bestInterMove = interPathVertexSwap(instance, sol, false);
            ScoredMove bestIntraMove = intraPathVertexSwap(instance, sol, false);

            // Apply best move
            if (bestInterMove.distanceDelta < bestIntraMove.distanceDelta)
            {
                std::swap(sol.path1[bestInterMove.vertex1], sol.path2[bestInterMove.vertex2]);
                sol.score += bestInterMove.distanceDelta;
            }
            else if (bestIntraMove.distanceDelta < bestInterMove.distanceDelta)
            {
                Path* pathForBestIntraMove = &sol.path1;
                if (bestIntraMove.pathIndex == 1)
                {
                    pathForBestIntraMove = &sol.path2;
                }

                std::swap(pathForBestIntraMove->at(bestIntraMove.vertex1), pathForBestIntraMove->at(bestIntraMove.vertex2));
                sol.score += bestInterMove.distanceDelta;
            }
            else
            {
                break;
            }

        }

        return sol;
    }
};

class SteepEdgeLocalSearch : public LocalSearch
{
public:
    const char* getName() { return "SteepEdgeLocalSearch"; }

    Solution run(const Instance& instance, const Solution& initialSolution)
    {
        Solution sol = initialSolution;
        sol.score = sol.getScore(instance);

        const auto& M = instance.M;

        while (true)
        {
            ScoredMove bestMove = edgeSwap(instance, sol, false);
            // Apply best move
            if (bestMove.distanceDelta < 0)
            {
                Path* pathForBestMove = &sol.path1;
                if (bestMove.pathIndex == 1)
                {
                    pathForBestMove = &sol.path2;
                }

                std::reverse(pathForBestMove->begin() + bestMove.vertex1, pathForBestMove->begin() + bestMove.vertex2 + 1);
                sol.score += bestMove.distanceDelta;
            }
            else
            {
                break;
            }

        }

        return sol;
    }
};

void test1(const std::filesystem::path& workDir, std::vector<Instance>& instances)
{
    std::cout << "solver, initializer, instance, score\n";

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
                int score = solution.getScore(instance);
                scores.push_back(score);
                if (score < bestScore)
                {
                    
                    bestScore = score;
                    bestSolution = solution;
                }
                 
            }
            
            int worstScore = *std::max_element(scores.begin(), scores.end());
            int avgScore = calculateMean(scores);

			std::string solFileName = instance.name + '-' + solver->getName() + "-solution.txt";

			bestSolution.dump(workDir / solFileName, instance);

            std::cout << solver->getName() << ", " << instance.name << ", " << avgScore << " (" << bestScore << '-' << worstScore << ")\n";
        }
    }
}

void test2(const std::filesystem::path& workDir, std::vector<Instance>& instances)
{
    std::cout << "algorithm, initializer, instance, score\n";

    LocalSearch* solvers[] = { new GreedyVertexLocalSearch, new GreedyEdgeLocalSearch, new SteepVertexLocalSearch, new SteepEdgeLocalSearch };
    TSPSolver* initializers[] = { new RandomSolver, new GreedyRegret };

    for (auto* solver : solvers)
    {
        for (auto* initializer : initializers)
        {
            for (auto& instance : instances)
            {
                std::vector<int> scores;
                int bestScore = std::numeric_limits<int>::max();
                Solution bestSolution;

                for (int i = 0; i < 100; ++i)
                {
                    instance.startIndex = i;
                    Solution initialSolution = initializer->run(instance);
                    Solution solution = solver->run(instance, initialSolution);
                    int initialScore = initialSolution.getScore(instance);
                    int score = solution.getScore(instance);
                    scores.push_back(score);
                    if (score < bestScore)
                    {
                        bestScore = score;
                        bestSolution = solution;
                    }
                }

                int worstScore = *std::max_element(scores.begin(), scores.end());
                int avgScore = calculateMean(scores);

                std::string solFileName = instance.name + '-' + solver->getName() + '-' + initializer->getName() + "-solution.txt";

                bestSolution.dump(workDir / solFileName, instance);

                std::cout << solver->getName() << ", " << initializer->getName() << ", " << instance.name << ", " << avgScore << " (" << bestScore << '-' << worstScore << ")\n";
            }
        }
    }
}

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

    //test1(workDir, instances);
    test2(workDir, instances);


    return 0;
}