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
#include <queue>
#include <numeric>

class Timer
{
public:
    void start()
    {
        m_startTime = std::chrono::high_resolution_clock::now();
    }

    void stop()
    {
        m_endTime = std::chrono::high_resolution_clock::now();
    }

    double elapsedMilliseconds() const
    {
        std::chrono::duration<double, std::milli> elapsed = m_endTime - m_startTime;
        return elapsed.count();
    }

    double elapsedMillisecondsNow() const
    {
        std::chrono::duration<double, std::milli> elapsed = std::chrono::high_resolution_clock::now() - m_startTime;
        return elapsed.count();
    }

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> m_startTime;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_endTime;
};

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

std::mt19937& getRandomGenerator()
{
    static unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //static std::mt19937 generator(seed);
    static std::mt19937 generator(0);

    return generator;
}

int getRandomNumber(int min, int max)
{
    std::uniform_int_distribution<int> distribution(min, max);
    return distribution(getRandomGenerator());
}

template <typename T>
void removeNRandomElements(std::vector<T>& vec, int n)
{
    for (int i = 0; i < n; ++i)
    {
        vec.erase(vec.begin() + getRandomNumber(0, vec.size()));
    }
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

template <typename T>
double calculateMean(const std::vector<T>& vec)
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
        score = 0;
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

        auto points = std::vector<int>(0, 100);
        for (int i = 0; i < instance.points.size(); ++i)
        {
            points.push_back(i);
        }

        while (!points.empty())
        {
            int randomPointIndex = getRandomNumber(0, points.size() - 1);
            sol.path1.push_back(points[randomPointIndex]);
            points.erase(points.begin() + randomPointIndex);

            randomPointIndex = getRandomNumber(0, points.size() - 1);
            sol.path2.push_back(points[randomPointIndex]);
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

        Timer timer;
        timer.start();

        while (timer.elapsedMillisecondsNow() >= 1000)
        {
            int move = getRandomNumber(0, 2);
            if (move == 0) // inter vertex swap
            {
				int pathIndex = getRandomNumber(0, 1);
                if (pathIndex == 0)
                {
					int v1 = getRandomNumber(0, sol.path1.size());
					int v2 = v1;
					while (v2 == v1)
					{
						v2 = getRandomNumber(0, sol.path1.size());
					}

					std::swap(sol.path1.at(v1), sol.path1.at(v2));
                }
                else
                {
					int v1 = getRandomNumber(0, sol.path2.size());
					int v2 = v1;
					while (v2 == v1)
					{
						v2 = getRandomNumber(0, sol.path2.size());
					}

					std::swap(sol.path1.at(v1), sol.path1.at(v2));
                }
            }
            else if (move == 1) // intra vertex swap
            {
                int v1 = getRandomNumber(0, sol.path1.size());
                int v2 = v1;
				while (v2 == v1)
				{
					v2 = getRandomNumber(0, sol.path1.size());
				}
            }
            else if (move == 0) // edge swap
            {
                int v1 = getRandomNumber(0, sol.path1.size());
                int v2 = getRandomNumber(0, sol.path2.size());
            }
        }


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
        for (int i = 1; i < n-1; ++i)
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
    int flag = 0;
    ScoredMove bestMove;
    bestMove.distanceDelta = 0;
    int pathIndex = 0;
    for (auto* pathPtr : { &sol.path1, &sol.path2 })
    {
        auto& path = *pathPtr;
        const int n = path.size();
        for (int i = 0; i < n - 2; ++i)
        {
            for (int j = i + 2; j < n - 1; ++j)
            {
                int v1 = path[i], v1After = path[(i+1)%n];
                int v2 = path[j], v2After = path[(j+1)%n];

                int distanceDelta = M[v1][v2] + M[v1After][v2After] - M[v1][v1After] - M[v2][v2After];

                if(distanceDelta == -2180)
                {   
                    // std::cout << "edge swap: " << distanceDelta << ' ' << v1 << ' ' << v1After << ' ' << v2 << ' ' << v2After <<' ' << pathIndex << '\n';
                    flag = 1;
                }
            
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

    if (flag)
    {
        // std::cout << "best move: " << bestMove.distanceDelta << ' ' << bestMove.vertex1 << ' ' << bestMove.vertex2 << ' ' << bestMove.pathIndex << '\n';
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

                std::reverse(pathForBestMove->begin() + bestMove.vertex1 + 1, pathForBestMove->begin() + bestMove.vertex2 + 1);
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
            ScoredMove bestInterMove = interPathVertexSwap(instance, sol, false);
            ScoredMove bestIntraMove = edgeSwap(instance, sol, false);
            // Apply best move
            if (bestInterMove.distanceDelta==bestIntraMove.distanceDelta && bestInterMove.distanceDelta == 0)
            {
                break;
            }
            else if (bestInterMove.distanceDelta <= bestIntraMove.distanceDelta)
            {
               // std::cout << "vertex swap: " << bestInterMove.distanceDelta << ' ' << bestInterMove.vertex1 << ' ' << bestInterMove.vertex2 << '\n';

                std::swap(sol.path1[bestInterMove.vertex1], sol.path2[bestInterMove.vertex2]);
                sol.score += bestInterMove.distanceDelta;
            }
            else if (bestIntraMove.distanceDelta < bestInterMove.distanceDelta)
            {
              //  std::cout << "edge swap: " << bestIntraMove.distanceDelta << ' ' << bestIntraMove.vertex1 << ' ' << bestIntraMove.vertex2 << ' ' << bestIntraMove.pathIndex << '\n';

                Path* pathForBestMove = &sol.path1;
                if (bestIntraMove.pathIndex == 1)
                {
                    pathForBestMove = &sol.path2;
                }

                std::reverse(pathForBestMove->begin() + bestIntraMove.vertex1 + 1, pathForBestMove->begin() + bestIntraMove.vertex2 + 1);
                sol.score += bestIntraMove.distanceDelta;
             
            }
            else
            {
                break;
            }

        }

        sol.getScore(instance);
        return sol;
    }
};

class SteepEdgeLocalSearchWithLM : public LocalSearch
{
public:
    const char* getName() { return "SteepEdgeLocalSearchWithLM"; }


    struct ScoredMoveLM
    {   
        struct EdgeData
        {
            int v1; // v1 -> v2; v1After -> v2After
            int v1After;
            int v2;
            int v2After;
            int pathIndex;
        };
        struct VertexData
        {
            int path1Vertices[3]; // before, vertextoswap, after
            int path2Vertices[3];
        };

        bool isedgeswap;
        int distanceDelta;
        union
        {
            EdgeData edgeData;
            VertexData vertexData;
        };

        bool operator<(const ScoredMoveLM& rhs)
        {
            return distanceDelta < rhs.distanceDelta;
        }
   };

   int cfind(const std::vector<int>& path, int value, int diff)
   {
       auto it = std::find(path.begin(), path.end(), value);
       if (it == path.end())
       {
          return -1;
       }

       int index = it - path.begin();

       if (diff == 0)
       {
           return *it;
       }
       else if (diff < 0)
       {
           if (index == 0)
           {
               return path[path.size() - 1];
           }
           
           return path[index + diff];
       }
       else if (diff > 0)
       {
           if (index == path.size() - 1)
           {
               return path[0];
           }

           return path[index + diff];
       }
       
       throw std::exception{};
   }

    using LMQueue = std::vector<ScoredMoveLM>;

    int getDistanceDeltaEdgeSwap(const DistanceMatrix& M, int v1, int v1After, int v2, int v2After)
    {
        if (v1 == v2)
        {
            return 0;
        }
        int distanceDelta = M[v1][v2] + M[v1After][v2After] - M[v1][v1After] - M[v2][v2After];
        return distanceDelta;
    }

    void initializeLM(LMQueue& LM, const Instance& instance, const Solution& sol)
    {
        const auto& M = instance.M;
        int pathIndex = 0;
        for (auto* pathPtr : { &sol.path1, &sol.path2 })
        {
            auto& path = *pathPtr;
            const int n = path.size();
            for (int i = 0; i < n - 2; ++i)
            {
                for (int j = i + 2; j < n - 1; ++j)
                {
                    int v1 = path[i], v1After = path[(i + 1) % n];
                    int v2 = path[j], v2After = path[(j + 1) % n];

                    int distanceDelta = getDistanceDeltaEdgeSwap(M, v1, v1After, v2, v2After);
                    
                    if (distanceDelta < 0)
                    {
                        ScoredMoveLM m;
                        m.isedgeswap = true;
                        m.edgeData = {v1, v1After, v2, v2After, pathIndex};
                        m.distanceDelta = distanceDelta;
                        LM.push_back(m);
                        m.edgeData = {v1After, v1, v2After, v2, pathIndex};
                        LM.push_back(m);
                    }
                    int distanceDelta2 = getDistanceDeltaEdgeSwap(M, v1, v1After, v2After, v2);
                    if (distanceDelta2 < 0)
                    {
                        ScoredMoveLM m;
                        m.isedgeswap = true;
                        m.edgeData = {v1, v1After, v2After, v2, pathIndex};
                        m.distanceDelta = distanceDelta2;
                        LM.push_back(m);
                        m.edgeData = {v1After, v1, v2, v2After, pathIndex};
                        LM.push_back(m);
                    }
                }
            }

            ++pathIndex;
        }

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

                if (distanceDelta < 0)
                {
                    ScoredMoveLM m;
                    m.isedgeswap = false;
                    m.vertexData = {{v1Before, v1, v1After}, {v2Before, v2, v2After}};
                    m.distanceDelta = distanceDelta;
                    LM.push_back(m);
                }

            }
        }
        
    }

    void updateLM(LMQueue& LM, const ScoredMoveLM& appliedMove, const Instance& instance, const Solution& sol)
    {   
        using namespace std;

        const auto& M = instance.M;
        if(appliedMove.isedgeswap)
        {
          
        
            auto& tmp = appliedMove.edgeData; 
            auto& path = tmp.pathIndex == 0 ? sol.path1 : sol.path2;
            int n = path.size();
            for (int i=0; i<2; i++)
            {
                //add new edge swap move for each swapped edge
                for (int j =0; j < n - 1; ++j)
                {
                    int v1, v2, v1After, v2After;
                    if(i)
                    {
                        v1 = tmp.v1After;
                        v1After = tmp.v2After;                      
                    }
                    else
                    {
                        v1 = tmp.v1;
                        v1After = tmp.v2;      
                    }
                    
                    v2 = path[j];
                    v2After = path[(j + 1) % n];

                    int distanceDelta = getDistanceDeltaEdgeSwap(M, v1, v1After, v2, v2After);

                    ScoredMoveLM m;
                    m.isedgeswap = true;
                    m.distanceDelta = distanceDelta;

                    if (distanceDelta < 0)
                    {
                            
                        m.edgeData = {v1, v1After, v2, v2After, tmp.pathIndex};
                        LM.push_back(m);
                        m.edgeData = {v1After, v1, v2After, v2, tmp.pathIndex};
                        LM.push_back(m);
                       
                    }
                    distanceDelta = getDistanceDeltaEdgeSwap(M, v1, v1After, v2After, v2);
                    if (distanceDelta < 0)
                    {
                        m.distanceDelta = distanceDelta;

                      
                        m.edgeData = {v1, v1After, v2After, v2, tmp.pathIndex};
                        LM.push_back(m);
                        m.edgeData = {v1After, v1, v2, v2After, tmp.pathIndex};
                        LM.push_back(m);
                    }
                    v2After = j==0 ? path[n-1] : path[j-1];
                    distanceDelta = getDistanceDeltaEdgeSwap(M, v1, v1After, v2, v2After);
                    if (distanceDelta < 0)
                    {
                        m.distanceDelta = distanceDelta;
                                            
                        
                        m.edgeData = {v1, v1After, v2, v2After, tmp.pathIndex};
                        
                        LM.push_back(m);
                        m.edgeData = {v1After, v1, v2After, v2, tmp.pathIndex};
                        LM.push_back(m);     
                    }
                    distanceDelta = getDistanceDeltaEdgeSwap(M, v1, v1After, v2After, v2);
                    if (distanceDelta < 0)
                    {

                      
                        m.distanceDelta = distanceDelta;
                        m.edgeData = {v1, v1After, v2After, v2, tmp.pathIndex};
                        LM.push_back(m);
                        m.edgeData = {v1After, v1, v2, v2After, tmp.pathIndex};
                        LM.push_back(m);
                    }
                }   


            }
            //add new vertex swap moves for each vertex in the swapped edge
            const auto& path2 = tmp.pathIndex == 0 ? sol.path2 : sol.path1;
            int n2 = path2.size();
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < n2; ++j)
                {
                    int v1, v1After, v1Before, index;
                    switch (i)
                    {
                    case 0:
                        v1 = tmp.v1;
                        
                        v1Before = cfind(path, v1, -1);
                        v1After = cfind(path, v1, 1);
                        break;

                    case 1:
                        v1 = tmp.v1After;
                        v1Before = cfind(path, v1, -1);
                        v1After = cfind(path, v1, 1);
                        break;

                    case 2:
                        v1 = tmp.v2;
                        v1Before = cfind(path, v1, -1);
                        v1After = cfind(path, v1, 1);
                        break;

                    case 3:
                        v1 = tmp.v2After;
                        v1Before = cfind(path, v1, -1);
                        v1After = cfind(path, v1, 1);
                        break;

                    }
                    

                    int v2 = sol.path2[j], v2Before = j == 0 ? sol.path2[n2 - 1] : sol.path2[j - 1], v2After = sol.path2[(j + 1) % n2];
                    int distanceNow = M[v1][v1After] + M[v1][v1Before] + M[v2][v2Before] + M[v2][v2After];
                    int distanceAfter = M[v1][v2After] + M[v1][v2Before] + M[v2][v1After] + M[v2][v1Before];

                    int distanceDelta = distanceAfter - distanceNow;

                    if (distanceDelta < 0)
                    {
                        ScoredMoveLM m;
                        m.isedgeswap = false;
                        if (tmp.pathIndex == 0)
                        {
                            m.vertexData = { {v1Before, v1, v1After}, {v2Before, v2, v2After} };
                        }
                        else
                        {
                            m.vertexData = { {v2Before, v2, v2After}, {v1Before, v1, v1After} };
                        }
                        m.distanceDelta = distanceDelta;
                        LM.push_back(m);
                    }
                }
            }
        }
        else if(!appliedMove.isedgeswap)
        {
            int Path1VertexAfterSwap= appliedMove.vertexData.path2Vertices[1];
            int Path2VertexAfterSwap= appliedMove.vertexData.path1Vertices[1];

            int P1v1 = Path1VertexAfterSwap;
            int P1v1After = cfind(sol.path1, Path1VertexAfterSwap, 1);
            int P1v12 = cfind(sol.path1, Path1VertexAfterSwap, -1);
            int P1v1After2 = P1v1;

            int P2v1 = Path2VertexAfterSwap;
            int P2v1After = cfind(sol.path2, Path2VertexAfterSwap, 1);
            int P2v12 = cfind(sol.path2, Path2VertexAfterSwap, -1);
            int P2v1After2 = P2v1;

            int v1, v1After;
            auto* pPath = &sol.path1;
            for(int i=0; i<4; i++)
            {
                switch(i)
                {
                    case 0:
                        v1 = P1v1;
                        v1After = P1v1After;
                        pPath = &sol.path1;
                        break;
                    case 1:
                        v1 = P1v1After;
                        v1After = P1v1After2;
                        pPath = &sol.path1;
                        break;
                    case 2:
                        v1 = P2v1;
                        v1After = P2v1After;
                        pPath = &sol.path2;
                        break;
                    case 3:
                        v1 = P2v1After;
                        v1After = P2v1After2;
                        pPath = &sol.path2;
                        break;
                }
                auto& path = *pPath;
                int n = path.size();
                for (int j = 0; j < n - 1; ++j)
                {
                    int v2 = path[j], v2After = path[(j + 1) % n];

                    int distanceDelta = getDistanceDeltaEdgeSwap(M, v1, v1After, v2, v2After);
                    if (distanceDelta < 0)
                    {
                        ScoredMoveLM m;
                        m.isedgeswap = true;
                        m.edgeData = {v1, v1After, v2, v2After, i < 2 ? 0 : 1};
                        m.distanceDelta = distanceDelta;
                        LM.push_back(m);
                        m.edgeData = {v1After, v1, v2After, v2, i < 2 ? 0 : 1};
                        LM.push_back(m);
                    }
                    distanceDelta = getDistanceDeltaEdgeSwap(M, v1, v1After, v2After, v2);
                    if (distanceDelta < 0)
                    {
                        ScoredMoveLM m;
                        m.isedgeswap = true;
                        m.edgeData = {v1, v1After, v2After, v2, i < 2 ? 0 : 1};
                        m.distanceDelta = distanceDelta;
                        LM.push_back(m);
                        m.edgeData = {v1After, v1, v2, v2After, i < 2 ? 0 : 1};
                        LM.push_back(m);
                    }
                    v2After = j == 0 ? path[n - 1] : path[j - 1];
                    distanceDelta = getDistanceDeltaEdgeSwap(M, v1, v1After, v2, v2After);
                    if (distanceDelta < 0)
                    {
                        ScoredMoveLM m;
                        m.isedgeswap = true;
                        m.edgeData = {v1, v1After, v2, v2After, i < 2 ? 0 : 1};
                        m.distanceDelta = distanceDelta;
                        LM.push_back(m);
                        m.edgeData = {v1After, v1, v2After, v2, i < 2 ? 0 : 1};
                        LM.push_back(m);
                    }
                    distanceDelta = getDistanceDeltaEdgeSwap(M, v1, v1After, v2After, v2);
                    if (distanceDelta < 0)
                    {
                        ScoredMoveLM m;
                        m.isedgeswap = true;
                        m.edgeData = {v1, v1After, v2After, v2, i < 2 ? 0 : 1};
                        m.distanceDelta = distanceDelta;
                        LM.push_back(m);
                        m.edgeData = {v1After, v1, v2, v2After, i < 2 ? 0 : 1};
                        LM.push_back(m);
                    }

                }
            }
            auto& tmp = appliedMove.vertexData;

            int verticesToCheck[6] ={tmp.path1Vertices[0],tmp.path2Vertices[1], tmp.path1Vertices[2], tmp.path2Vertices[0], tmp.path1Vertices[1], tmp.path2Vertices[2]};
            
            for(int i=0; i<6; i++)
            {   // 4 8 15 3 19 18
                auto& path = i>=3 ? sol.path2 : sol.path1;
                int v1 = verticesToCheck[i];
                int v1Index = std::find(path.begin(), path.end(), v1) - path.begin();

                int v1Before = cfind(path, v1, -1);
                int v1After = cfind(path, v1, 1);
                int n = i>=3 ? sol.path1.size() : sol.path2.size();

                for(int j=0; j<n; j++)
                {
                    int v2 = i>=3 ? sol.path1[j] : sol.path2[j];
                    int v2Before = j==0? (i>=3 ? sol.path1[sol.path1.size()-1] : sol.path2[sol.path2.size()-1]) : (i>=3 ? sol.path1[j-1] : sol.path2[j-1]);
                    int v2After = i>=3 ? sol.path1[(j+1)%n] : sol.path2[(j+1)%n];
                    int distanceNow = M[v1][v1After] + M[v1][v1Before] + M[v2][v2Before] + M[v2][v2After];
                    int distanceAfter = M[v1][v2After] + M[v1][v2Before] + M[v2][v1After] + M[v2][v1Before];

                    int distanceDelta = distanceAfter - distanceNow;

                    if (distanceDelta < 0)
                    {
                        ScoredMoveLM m;
                        m.isedgeswap = false;
                        if(i>=3)
                        {
                            m.vertexData = {{v2Before, v2, v2After}, {v1Before, v1, v1After}};
                        }
                        else
                        {
                            m.vertexData = {{v1Before, v1, v1After}, {v2Before, v2, v2After}};
                        }
                    
                        m.distanceDelta = distanceDelta;

                        LM.push_back(m);
                    }
                }

            }

        }
        else
        {
            std::cerr << "Invalid move type\n";
            throw std::exception{};
        }
    }
    
    ScoredMoveLM findBestApplicableMove(LMQueue& LM, const Instance& instance, const Solution& sol)
    {
        std::sort(LM.begin(), LM.end());

        ScoredMoveLM bestMove;
        bestMove.distanceDelta = 0;

        for (auto move = LM.begin(); move != LM.end(); )
        {
            if (move->isedgeswap)
            {
                if(move->distanceDelta == -2180)
                {
                    // std::cout << "here\n";
                }
               
                auto& tmp = move->edgeData;
                auto& path1 = tmp.pathIndex == 0 ? sol.path1 : sol.path2;
                auto p1V1After = cfind(path1, tmp.v1, 1);
                auto p2V1After = cfind(path1, tmp.v2, 1);
                if(p1V1After == tmp.v2 || p2V1After == tmp.v1)
                {
                    move = LM.erase(move);
                    
                }
                else if (p1V1After == tmp.v1After && p2V1After == tmp.v2After)
                {
                    bestMove = *move;
                    move = LM.erase(move);
                    break;
                }
                else if (cfind(path1, tmp.v1, 1) == tmp.v1After && cfind(path1, tmp.v2After, 1) == tmp.v2 ||
                    cfind(path1, tmp.v1After, 1) == tmp.v1 && cfind(path1, tmp.v2, 1) == tmp.v2After)
                {
                    ++move;
                }
                else if (cfind(path1,tmp.v1After,1) == tmp.v1 && cfind(path1, tmp.v2After, 1) == tmp.v2)
                {
                    move = LM.erase(move);
                }
                else
                {
                    move = LM.erase(move);
                }
            }
            else
            {
                int P1=move->vertexData.path1Vertices[1];
                int P2=move->vertexData.path2Vertices[1];
                int P1Before = cfind(sol.path1, P1, -1);
                int P1After = cfind(sol.path1, P1, 1);
                int P2Before = cfind(sol.path2, P2, -1);
                int P2After = cfind(sol.path2, P2, 1);
                
                if(cfind(sol.path1,P1,0)==-1 || cfind(sol.path2,P2,0)==-1)
                {
                   
                    move = LM.erase(move);
                    
                }
                else if ((P1Before == move->vertexData.path1Vertices[0] && P1After == move->vertexData.path1Vertices[2] ||
                    P1Before == move->vertexData.path1Vertices[2] && P1After == move->vertexData.path1Vertices[0]) &&
                    (P2Before == move->vertexData.path2Vertices[0] && P2After == move->vertexData.path2Vertices[2] ||
                    P2Before == move->vertexData.path2Vertices[2] && P2After == move->vertexData.path2Vertices[0]))
                {
                    bestMove = *move;
                    move = LM.erase(move);
                    break;
                }
                else 
                {
                   
                    move = LM.erase(move);
                }
            }
        }

        return bestMove;
    }

    Solution run(const Instance& instance, const Solution& initialSolution)
    {
        Solution sol = initialSolution;
        sol.score = sol.getScore(instance);
        const auto& M = instance.M;

        LMQueue LM;
        initializeLM(LM, instance, sol);

        while (true)
        {
            ScoredMoveLM bestMove = findBestApplicableMove(LM, instance, sol);
            if(bestMove.isedgeswap)
            {
               // std::cout << "edge swap: " << bestMove.distanceDelta << ' ' << bestMove.edgeData.v1 << ' ' << bestMove.edgeData.v2 << ' ' << bestMove.edgeData.pathIndex << '\n';
            }
            else
            {
               // std::cout << "vertex swap: " << bestMove.distanceDelta << ' ' << bestMove.vertexData.path1Vertices[1] << ' ' << bestMove.vertexData.path2Vertices[1] << '\n';
            }

            // Apply move
            if (bestMove.isedgeswap && bestMove.distanceDelta < 0)
            {
                Path* pathForBestMove = &sol.path1;
                if (bestMove.edgeData.pathIndex == 1)
                {
                    pathForBestMove = &sol.path2;
                }

                int i = std::find(pathForBestMove->begin(), pathForBestMove->end(), bestMove.edgeData.v1) - pathForBestMove->begin();
                int j = std::find(pathForBestMove->begin(), pathForBestMove->end(), bestMove.edgeData.v2) - pathForBestMove->begin();
                if (i > j)
                {
                    std::swap(i, j);
                }
                // std::reverse(std::find(pathForBestMove->begin(), pathForBestMove->end(), bestMove.edgeData.v1) + 1, std::find(pathForBestMove->begin(), pathForBestMove->end(), bestMove.edgeData.v2) + 1);
                std::reverse(pathForBestMove->begin() + i + 1, pathForBestMove->begin() + j + 1);
                int actualScore = sol.getScore(instance);
                int actualDiff = actualScore - sol.score;
                if (actualDiff == bestMove.distanceDelta)
                {
                    sol.score += bestMove.distanceDelta;
                }
                else
                {
                    // std::cout << "edge swap: " << bestMove.distanceDelta << ' ' << bestMove.edgeData.v1 << ' ' << bestMove.edgeData.v1After << ' ' << bestMove.edgeData.v2 << ' ' << bestMove.edgeData.v2After << '\n';
                    for(auto& v : sol.path1)
                    {
                        std::cout << v << ' ';
                    }
                    // std::cout << '\n';
                    for(auto& v : sol.path2)
                    {
                        std::cout << v << ' ';
                    }
                    // std::cout << '\n';
                    return sol;
                    throw std::exception{};
                }

                updateLM(LM, bestMove, instance, sol);
            }
            else if (bestMove.distanceDelta < 0)
            {
                std::swap(*std::find(sol.path1.begin(), sol.path1.end(), bestMove.vertexData.path1Vertices[1]), *std::find(sol.path2.begin(), sol.path2.end(), bestMove.vertexData.path2Vertices[1]));
                // std::swap(sol.path1[cfind(sol.path1,bestMove.vertexData.path1Vertices[1])], sol.path2[bestMove.vertexData.path1Vertices[1]]);

                int actualScore = sol.getScore(instance);
                int actualDiff = actualScore - sol.score;
                if (actualDiff == bestMove.distanceDelta)
                {
                    sol.score += bestMove.distanceDelta;
                }
                else
                {
                    throw std::exception{};
                }

                updateLM(LM, bestMove, instance, sol);
            }
            else
            {
                break;
            }

        }

        return sol;
    }
};

class CandidateLocalSearch : public LocalSearch
{
public:
    const char* getName() { return "CandidateLocalSearch"; }

    static constexpr int NumCandidates = 10;

    using ClosestVertices = std::vector<std::vector<int>>;

    // vertex swap
    ScoredMove candidateInter(const Instance& instance, const Solution& sol, const ClosestVertices& cv)
    {
        const auto& M = instance.M;
        ScoredMove bestInterMove;
        bestInterMove.distanceDelta = 0;
        const int n1 = sol.path1.size();
        const int n2 = sol.path2.size();
        int pathIndex = 0;
        for (auto* pathPtr : { &sol.path1, &sol.path2 })
        {
            for (int i = 0; i < pathPtr->size(); ++i)
            {
                int v1 = pathPtr->at(i), v1Before = i == 0 ? pathPtr->at(n1 - 1) : pathPtr->at(i - 1), v1After = pathPtr->at((i + 1) % n1);

                for (int v2 : cv[v1])
                {
                    auto* path2Ptr = pathPtr == &sol.path1 ? &sol.path2 : &sol.path1;
                    auto it = std::find(path2Ptr->begin(), path2Ptr->end(), v2);
                    if (it == path2Ptr->end())
                    {
                        continue;
                    }

                    int j = it - path2Ptr->begin();
                    int v2Before = j == 0 ? path2Ptr->at(n2 - 1) : path2Ptr->at(j - 1), v2After = path2Ptr->at((j + 1) % n2);
                    int distanceNow = M[v1][v1After] + M[v1][v1Before] + M[v2][v2Before] + M[v2][v2After];
                    int distanceAfter = M[v1][v2After] + M[v1][v2Before] + M[v2][v1After] + M[v2][v1Before];

                    int distanceDelta = distanceAfter - distanceNow;

                    if (distanceDelta < bestInterMove.distanceDelta)
                    {
                        bestInterMove = ScoredMove{ i, j, distanceDelta, pathIndex };
                    }
                }
            }

            ++pathIndex;
        }

        return bestInterMove;
    }

    // edge swap
    ScoredMove candidateIntra(const Instance& instance, const Solution& sol, const ClosestVertices& cv)
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
                int v1 = path[i], v1After = path[(i+1)%n];

                for (int v2 : cv[v1])
                {
                    auto it = std::find(pathPtr->begin(), pathPtr->end(), v2);
                    if (it == pathPtr->end())
                    {
                        continue;
                    }

                    int j = it - pathPtr->begin();
                    int v2After = path[(j+1)%n];

                    int distanceDelta = M[v1][v2] + M[v1After][v2After] - M[v1][v1After] - M[v2][v2After];

                    if (distanceDelta < bestMove.distanceDelta)
                    {
                        bestMove = ScoredMove{ i, j, distanceDelta, pathIndex };
                    }
                }
            }

            ++pathIndex;
        }

        return bestMove;
    }

    std::vector<size_t> getSmallestIndices(const std::vector<int>& vec, size_t n)
    {
        std::vector<int> sortedVec = vec;
        std::sort(sortedVec.begin(), sortedVec.end(), std::less<int>());
        sortedVec.resize(std::min(n, sortedVec.size()));
        std::vector<size_t> indices;
        for (auto& val : sortedVec)
        {
            auto it = std::find(vec.begin(), vec.end(), val);
            if (it != vec.end())
            {
                indices.push_back(std::distance(vec.begin(), it));
            }
        }

        return indices;
    }

    Solution run(const Instance& instance, const Solution& initialSolution)
    {
        Solution sol = initialSolution;
        sol.score = sol.getScore(instance);

        const auto& M = instance.M;
        ClosestVertices CV = ClosestVertices(M.size(), std::vector<int>(NumCandidates, 0));

        for (int i = 0; i < M.size(); ++i)
        {
            const auto smallestIndices = getSmallestIndices(M[i], NumCandidates+1);
            for (int j = 0; j < NumCandidates; ++j)
            {
                CV[i][j] = smallestIndices[j+1]; // +1 to skip distance to self
            }
        }

        while (true)
        {
            ScoredMove bestInterMove = candidateInter(instance, sol, CV);
            ScoredMove bestIntraMove = candidateIntra(instance, sol, CV);

            // Apply best move
            if (bestInterMove.distanceDelta==bestIntraMove.distanceDelta && bestInterMove.distanceDelta == 0)
            {
                break;
            }
            else if (bestInterMove.distanceDelta <= bestIntraMove.distanceDelta)
            {
                if (bestInterMove.pathIndex == 0)
                {
                    std::swap(sol.path1[bestInterMove.vertex1], sol.path2[bestInterMove.vertex2]);
                }
                else
                {
                    std::swap(sol.path2[bestInterMove.vertex1], sol.path1[bestInterMove.vertex2]);
                }

                int actualScore = sol.getScore(instance);
                int actualDiff = actualScore - sol.score;
                if (actualDiff == bestInterMove.distanceDelta)
                {
                    sol.score += bestInterMove.distanceDelta;
                }
                else
                {
                    throw std::exception{};
                }
            }
            else if (bestIntraMove.distanceDelta < bestInterMove.distanceDelta)
            {
                Path* pathForBestMove = &sol.path1;
                if (bestIntraMove.pathIndex == 1)
                {
                    pathForBestMove = &sol.path2;
                }

                if (bestIntraMove.vertex1 < bestIntraMove.vertex2)
                {
                    std::reverse(pathForBestMove->begin() + bestIntraMove.vertex1 + 1, pathForBestMove->begin() + bestIntraMove.vertex2 + 1);
                }
                else
                {
                    std::reverse(pathForBestMove->begin() + bestIntraMove.vertex2 + 1, pathForBestMove->begin() + bestIntraMove.vertex1 + 1);
                }

                int actualScore = sol.getScore(instance);
                int actualDiff = actualScore - sol.score;
                if (actualDiff == bestIntraMove.distanceDelta)
                {
                    sol.score += bestIntraMove.distanceDelta;
                }
                else
                {
                    throw std::exception{};
                }
            }
            else
            {
                break;
            }

        }

        return sol;
    }
};

class MultipleStartLocalSearch : public TSPSolver
{
public:
    const char* getName() { return "MultipleStartLocalSearch"; }

    using LocalSearchAlgorithm = SteepEdgeLocalSearch;
    using LocalSearchInitializer = RandomSolver;
    static constexpr int iterCount = 100;

    Solution run(const Instance& instance)
    {
        LocalSearchInitializer lsInitializer;
        LocalSearchAlgorithm lsAlgo;

        Solution bestSol;

        for (int i = 0; i < iterCount; ++i)
        {
            Solution startSol = lsInitializer.run(instance);
            Solution sol = lsAlgo.run(instance, startSol);

            if (bestSol.score == 0 || sol.score < bestSol.score)
            {
                bestSol = sol;
            }

            std::cout << "MSLS iter: " << i+1 << '/' << iterCount << "   \r";
        }

        return bestSol;
    }
};

class IteratedLocalSearch1 : public TSPSolver
{
public:
    const char* getName() { return "IteratedLocalSearch1"; }

    using LocalSearchAlgorithm = SteepEdgeLocalSearch;
    using LocalSearchInitializer = RandomSolver;
    static constexpr int iterCount = 1800;

    void edgeSwap(Path* path, int v1, int v2)
    {
        if (v1 > v2)
        {
            std::swap(v1, v2);
        }
        std::reverse(path->begin() + v1 + 1, path->begin() + v2 + 1);
    }

    void vertexSwap(Solution& sol, int v1, int v2)
    {
        std::swap(sol.path1[v1], sol.path2[v2]);
    }

    void perturb(const Instance& instance, Solution& solution)
    {
        int numOfEdgeSwaps = getRandomNumber(2, 4);
        int numOfVertexSwaps = getRandomNumber(2, 4);
        Path* paths[] = { &solution.path1, &solution.path2 };

        // random order of vertex or edge swaps
        while (numOfEdgeSwaps > 0 || numOfVertexSwaps > 0)
        {
            if ((getRandomNumber(0, 1) && numOfEdgeSwaps > 0) || numOfVertexSwaps == 0)
            {
                // edge swap
                int pathIndex = getRandomNumber(0, 1);
                Path* path = paths[getRandomNumber(0, 1)];
                int v1 = getRandomNumber(0, path->size() - 1);
                int v2 = v1;
                while (v2 == v1) // avoid chosing v1 == v2
                {
                    v2 = getRandomNumber(0, path->size() - 1);
                }
                edgeSwap(path, v1, v2);
                --numOfEdgeSwaps;
            }
            else
            {
                // vertex swap
                int v1 = getRandomNumber(0, paths[0]->size() - 1);
                int v2 = getRandomNumber(0, paths[1]->size() - 1);
                vertexSwap(solution, v1, v2);
                --numOfVertexSwaps;
            }
        }
    }

    Solution run(const Instance& instance)
    {
        LocalSearchInitializer lsInitializer;
        LocalSearchAlgorithm lsAlgo;

        Solution startSol = lsInitializer.run(instance);

        Solution bestSol = lsAlgo.run(instance, startSol);
        bestSol.getScore(instance);

        for (int i = 0; i < iterCount; ++i)
        {
            Solution sol = bestSol;
            perturb(instance, sol);

            sol = lsAlgo.run(instance, sol);
            sol.getScore(instance);

            if (sol.score < bestSol.score)
            {
                bestSol = sol;
            }

            std::cout << "ILS1 iter: " << i+1 << '/' << iterCount << "   \r";
        }

        return bestSol;
    }
};

class IteratedLocalSearch2 : public TSPSolver
{
public:
    const char* getName() { return "IteratedLocalSearch2"; }

    using LocalSearchAlgorithm = SteepEdgeLocalSearch;
    using LocalSearchInitializer = RandomSolver;
    static constexpr int iterCount = 1400;
    std::set<int> unVisited;
    std::array<std::vector<int>, 2> paths;
    

    void destroy(Solution& solution, const Instance& instance)
    {   
        const auto& M = instance.M;
        // TODO: also remove edges, select better places for removal
        int n = getRandomNumber(20, 30);
    
        for (int i = 0; i < n; ++i)
        {
            int v = getRandomNumber(0, solution.path1.size() - 1);
            unVisited.insert(solution.path1[v]);
            solution.path1.erase(solution.path1.begin() + v);
            unVisited.insert(solution.path2[v]);
            solution.path2.erase(solution.path2.begin() + v);
        }
        
    }

    void repair(const Instance& instance, Solution& solution)
    {
        paths[0] = solution.path1;
        paths[1] = solution.path2;
        const auto& M = instance.M;
        const auto& points = instance.points;
        
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
        solution.path1 = paths[0];
        solution.path2 = paths[1];
    }

    Solution run(const Instance& instance)
    {
        LocalSearchInitializer lsInitializer;
        LocalSearchAlgorithm lsAlgo;

        Solution startSol = lsInitializer.run(instance);

        Solution bestSol = lsAlgo.run(instance, startSol);

        for (int i = 0; i < iterCount; ++i)
        {
            Solution sol = bestSol;
            destroy(sol,instance);
            repair(instance, sol);

            sol = lsAlgo.run(instance, sol);

            if (sol.score < bestSol.score)
            {
                bestSol = sol;
            }

            std::cout << "ILS2 iter: " << i+1 << '/' << iterCount << "   \r";
        }

        return bestSol;
    }
};

class HAE : public TSPSolver
{
public:
    const char* getName() { return "HAE"; }

    using LocalSearchAlgorithm = SteepEdgeLocalSearch;
    using LocalSearchInitializer = RandomSolver;
    static constexpr int populationSize = 20;
    static constexpr int iterCount = 1400;

    static bool solutionExists(const std::vector<Solution>& solutions, const Solution& solution)
    {
        for (const auto& sol : solutions)
        {
            if (sol.score == solution.score)
            {
                return true;
            }
        }
        
        return false;
    }

    static Solution recombine(const Solution& sol1, const Solution& sol2, const Instance& instance)
    {
        // randomly choose parent solution
        Solution s1;
        const Solution* s2;
        if (getRandomNumber(0, 1))
        {
            s1 = sol1;
            s2 = &sol2;
        }
        else
        {
            s1 = sol2;
            s2 = &sol1;
        }

        std::vector<Path*> s1Paths = { &s1.path1, &s1.path2 };


        // remove not shared parts

        // first find common parts
        std::vector<std::vector<bool>> commonParts(s1Paths.size(), std::vector<bool>(s1.path1.size()));

        for (int path1Idx = 0; path1Idx < s1Paths.size(); ++path1Idx)
        {
            auto* path1 = s1Paths[path1Idx];
            std::vector<bool>& commonPart = commonParts[path1Idx];

            for (int i = 0; i < path1->size() - 1; ++i)
            {
                for (const auto* path2 : std::vector<const Path*>{ &s2->path1, &s2->path2 })
                {
                    if (i >= path1->size() - 1)
                    {
                        break;
                    }

                    for (int j = 0; j < path2->size() - 1; ++j)
                    {
                        if (path1->at(i) == path2->at(j) && path1->at(i + 1) == path2->at(j + 1))
                        {
                            while (i < path1->size() && j < path2->size() && path1->at(i) == path2->at(j))
                            {
                                commonPart[i] = true;
                                ++i;
                                ++j;
                            }
                            
                            break;
                        }
                    }
                }
            }
        }

        std::set<int> unVisited;

        // then set -1 in places that are not common
        for (int pathIdx = 0; pathIdx < s1Paths.size(); ++pathIdx)
        {
            auto* path = s1Paths[pathIdx];
            std::vector<bool>& commonPart = commonParts[pathIdx];

            for (int i = 0; i < path->size() - 1; ++i)
            {
                if (!commonPart[i])
                {
                    unVisited.insert(path->at(i));
                    path->at(i) = -1;
                }
            }
        }

        // remove isolated points
        for (auto* path : std::vector<Path*>{ &s1.path1, &s1.path2 })
        {
            for (int i = 0; i < path->size(); ++i)
            {
                if (path->at(i) == -1)
                {
                    continue;
                }

                int n1, n2;
                if (i == 0)
                {
                    n1 = path->size() - 1;
                    n2 = i + 1;
                }
                else if (i == path->size() - 1)
                {
                    n1 = i - 1;
                    n2 = 0;
                }
                else
                {
                    n1 = i - 1;
                    n2 = i + 1;
                }

                if (path->at(n1) == -1 && path->at(n2) == -1)
                {
                    unVisited.insert(path->at(i));
                    path->at(i) = -1;
                }
            }
        }

        // remove the empty spaces
        s1.path1.erase(std::remove(s1.path1.begin(), s1.path1.end(), -1), s1.path1.end());
        s1.path2.erase(std::remove(s1.path2.begin(), s1.path2.end(), -1), s1.path2.end());

        // repair
        const int desiredPathLen = sol2.path1.size();
        while (!unVisited.empty())
        {
            for (auto* path : std::vector<Path*>{ &s1.path1, &s1.path2 })
            {
                if (path->size() == 0)
                {
                    path->push_back(*unVisited.begin());
                    unVisited.erase(unVisited.begin());
                }

                if (path->size() == desiredPathLen)
                {
                    continue;
                }

				int minLenDiff = std::numeric_limits<int>::max();
				int bestPlacement = -1;
				int bestPoint = -1;

                for (int i = 0; i < path->size(); ++i)
                {
					for (int pointIndex : unVisited)
					{
                        int lenDiff = getLenDiff(i, pointIndex, *path, instance.M);
						if (lenDiff < minLenDiff)
						{
							minLenDiff = lenDiff;
							bestPoint = pointIndex;
                            bestPlacement = i;
						}
					}
                }

				unVisited.erase(bestPoint);
                path->insert(path->begin() + bestPlacement + 1, bestPoint);
            }
        }

        s1.getScore(instance);
        return s1;
    }

    static int getWorstSolIdx(const std::vector<Solution>& solutions)
    {
        int worstScoreIdx = 0;
        int worstScore = solutions[worstScoreIdx].score;

        for (int i = 0; i < solutions.size(); ++i)
        {
            if (solutions[i].score > worstScore)
            {
                worstScore = solutions[i].score;
                worstScoreIdx = i;
            }
        }

        return worstScoreIdx;
    }

    static int getBestSolIdx(const std::vector<Solution>& solutions)
    {
        int bestScoreIdx = 0;
        int bestScore = solutions[bestScoreIdx].score;

        for (int i = 0; i < solutions.size(); ++i)
        {
            if (solutions[i].score < bestScore)
            {
                bestScore = solutions[i].score;
                bestScoreIdx = i;
            }
        }

        return bestScoreIdx;
    }

    static bool isDifferentEnough(const std::vector<Solution>& solutions, const Solution& solution)
    {
        for (const auto& sol : solutions)
        {
            if (sol.score == solution.score)
            {
                return false;
            }
        }

        /*
        std::vector<const Path*> s1Paths = { &solution.path1, &solution.path2 };
        for (const auto* path1 : s1Paths)
        {
            for (const auto& sol : solutions)
            {
                std::vector<const Path*> s2Paths = { &sol.path1, &sol.path2 };
                for (const auto* path2 : s2Paths)
                {
                    int countDiff = 0;
                    for (const auto& point : *path1)
                    {
                        if (std::find(path2->begin(), path2->end(), point) == path2->end())
                        {
                            ++countDiff;
                            if (countDiff > 0.05 * path1->size())
                            {
                                break;
                            }
                        }

                    }

                    if (countDiff <= 0.05 * path1->size())
                    {
                        return false;
                    }
                }
            }
        }
        */

        return true;
    }
   
    Solution run(const Instance& instance)
    {
        static std::vector<int> iters;
        auto start = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::milliseconds(1800);
        if (instance.name == "kroB200.tsp")
        {
            duration = std::chrono::milliseconds(1700);
        }

        std::vector<Solution> population;
        // initial population
        while (population.size() < populationSize)
        {
            Solution newSolution = LocalSearchAlgorithm().run(instance, LocalSearchInitializer().run(instance));
            if (!solutionExists(population, newSolution))
            {
                population.push_back(newSolution);
            }
        }

        int diffCount = 0;
        int sameCount = 0;

        int i = 0;

        while (std::chrono::high_resolution_clock::now() - start < duration)
        {
            int sol1Idx = getRandomNumber(0, population.size() - 1);
            int sol2Idx = sol1Idx;
            while (sol2Idx == sol1Idx)
            {
                sol2Idx = getRandomNumber(0, population.size() - 1);
            }

            const Solution& sol1 = population[sol1Idx];
            const Solution& sol2 = population[sol2Idx];

            Solution newSol = recombine(sol1, sol2, instance);
            newSol = LocalSearchAlgorithm().run(instance, newSol);
    
            int worstSolIdx = getWorstSolIdx(population);

            if (newSol.getScore(instance) < population[worstSolIdx].score && isDifferentEnough(population, newSol))
            {
                population.erase(population.begin() + worstSolIdx);
                population.push_back(newSol);
                ++diffCount;
            }
            else
            {
                ++sameCount;
            }

            ++i;

            std::cout << "Accepted/Rejected: " << diffCount << "/" << sameCount << "     \r";
        }

        iters.push_back(i);
        int sum = std::accumulate(iters.begin(), iters.end(), 0);
        double average = static_cast<double>(sum) / iters.size();
        std::cout << "\niters average: " << average << "\n";

        return population[getBestSolIdx(population)];
    }
};

void test1(const std::filesystem::path& workDir, std::vector<Instance>& instances)
{
    std::cout << "solver, initializer, instance, score\n";

    TSPSolver* solvers[] = { new GreedyRegret };
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
    std::cout << "algorithm, initializer, instance, score, avg time\n";

    Timer timer;
    LocalSearch* solvers[] = { new GreedyVertexLocalSearch, new GreedyEdgeLocalSearch, new SteepVertexLocalSearch, new SteepEdgeLocalSearch };
    TSPSolver* initializers[] = { new RandomSolver, new GreedyRegret };

    for (auto* solver : solvers)
    {
        for (auto* initializer : initializers)
        {
            for (auto& instance : instances)
            {
                std::vector<int> scores;
                std::vector<double> times;
                int bestScore = std::numeric_limits<int>::max();
                Solution bestSolution;

                for (int i = 0; i < 100; ++i)
                {
                    instance.startIndex = i;
                    Solution initialSolution = initializer->run(instance);
                    timer.start();
                    Solution solution = solver->run(instance, initialSolution);
                    timer.stop();
                    times.push_back(timer.elapsedMilliseconds());
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
                double bestTime = *std::min_element(times.begin(), times.end());
                double worstTime = *std::max_element(times.begin(), times.end());
                double avgTime = calculateMean(times);

                std::string solFileName = instance.name + '-' + solver->getName() + '-' + initializer->getName() + "-solution.txt";

                bestSolution.dump(workDir / solFileName, instance);

                std::cout << solver->getName() << ", " << initializer->getName() << ", " << instance.name << ", " << avgScore << " (" << bestScore << '-' << worstScore << "), " << avgTime << " (" << bestTime << '-' << worstTime << ") ms\n";
            }
        }
    }
}

void test3(const std::filesystem::path& workDir, std::vector<Instance>& instances)
{
    std::cout << "algorithm, initializer, instance, score, avg time\n";

    Timer timer;
    LocalSearch* solvers[] = { new SteepEdgeLocalSearchWithLM };
    TSPSolver* initializers[] = { new RandomSolver, new GreedyRegret };

    Solution debugInitialSolution = initializers[0]->run(instances[0]);
    for (auto* solver : solvers)
    {
        for (auto* initializer : initializers)
        {
            for (auto& instance : instances)
            {
                std::vector<int> scores;
                std::vector<double> times;
                int bestScore = std::numeric_limits<int>::max();
                Solution bestSolution;

                for (int i = 0; i < 100; ++i)
                {   
                    instance.startIndex = i;
                    Solution initialSolution = initializer->run(instance);
                    //Solution initialSolution = debugInitialSolution;
                    timer.start();
                    Solution solution = solver->run(instance, initialSolution);
                    timer.stop();
                    times.push_back(timer.elapsedMilliseconds());
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
                double bestTime = *std::min_element(times.begin(), times.end());
                double worstTime = *std::max_element(times.begin(), times.end());
                double avgTime = calculateMean(times);

                std::string solFileName = instance.name + '-' + solver->getName() + '-' + initializer->getName() + "-solution.txt";

                bestSolution.dump(workDir / solFileName, instance);

                std::cout << solver->getName() << ", " << initializer->getName() << ", " << instance.name << ", " << avgScore << " (" << bestScore << '-' << worstScore << "), " << avgTime << " (" << bestTime << '-' << worstTime << ") ms\n";
            }
        }
    }
}

void test4(const std::filesystem::path& workDir, std::vector<Instance>& instances)
{
    std::cout << "algorithm, instance, score, avg time\n";

    Timer timer;
    TSPSolver* solvers[] = { new MultipleStartLocalSearch, new IteratedLocalSearch1, new IteratedLocalSearch2 };

    for (auto* solver : solvers)
    {
        for (auto& instance : instances)
        {
            std::vector<int> scores;
            std::vector<double> times;
            int bestScore = std::numeric_limits<int>::max();
            Solution bestSolution;

            for (int i = 0; i < 10; ++i)
            {   
                instance.startIndex = i;
                timer.start();
                Solution solution = solver->run(instance);
                timer.stop();
                times.push_back(timer.elapsedMilliseconds());
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
            double bestTime = *std::min_element(times.begin(), times.end());
            double worstTime = *std::max_element(times.begin(), times.end());
            double avgTime = calculateMean(times);

            std::string solFileName = instance.name + '-' + solver->getName() + "-solution.txt";

            bestSolution.dump(workDir / solFileName, instance);

            std::cout << solver->getName() << ", " << instance.name << ", " << avgScore << " (" << bestScore << '-' << worstScore << "), " << avgTime << " (" << bestTime << '-' << worstTime << ") ms\n";
        }
    }
}

void test5(const std::filesystem::path& workDir, std::vector<Instance>& instances)
{
    std::cout << "algorithm, instance, score, avg time\n";

    Timer timer;
    TSPSolver* solvers[] = { new GreedyCycle, new HAE };

    for (auto* solver : solvers)
    {
        for (auto& instance : instances)
        {
            std::vector<int> scores;
            std::vector<double> times;
            int bestScore = std::numeric_limits<int>::max();
            Solution bestSolution;

            for (int i = 0; i < 10; ++i)
            {   
                instance.startIndex = i;
                timer.start();
                Solution solution = solver->run(instance);
                timer.stop();
                times.push_back(timer.elapsedMilliseconds());
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
            double bestTime = *std::min_element(times.begin(), times.end());
            double worstTime = *std::max_element(times.begin(), times.end());
            double avgTime = calculateMean(times);

            std::string solFileName = instance.name + '-' + solver->getName() + "-solution.txt";

            bestSolution.dump(workDir / solFileName, instance);

            std::cout << solver->getName() << ", " << instance.name << ", " << avgScore << " (" << bestScore << '-' << worstScore << "), " << avgTime << " (" << bestTime << '-' << worstTime << ") ms\n";
        }
    }
}

int main(int argc, char* argv[])
{
    std::filesystem::path workDir = argc > 1 ? argv[1] : "workdir";


    std::vector<Instance> instances;
    std::string instanceNames[] = { "kroA200.tsp.txt", "kroB200.tsp.txt" };
    for (const auto& instanceName : instanceNames)
    {
        instances.emplace_back();
        instances[instances.size() - 1].load(workDir / instanceName);
    }

   // test3(workDir, instances);
     test5(workDir, instances);


    return 0;
}