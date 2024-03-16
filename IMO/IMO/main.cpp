#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <vector>
#include <math.h>
#include <filesystem>

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
        int start2Index = getFurthestindex(start1Index, M, dim);
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

    int getFurthestindex(int index, const DistanceMatrix& M, int dim)
    {
        int maxdistance = 0;
        int maxindex = -1;

        for (int i =0; i < dim; i++)
        {
            if (M[index][i] > maxdistance)
            {
                maxdistance = M[index][i];
                maxindex = i;
            }
        }
        return maxindex;
    }
};


class  Edge
{
    int from;
    int to;
    int distance;
    public:
    Edge(int from, int to, int distance) : from(from), to(to), distance(distance) {}
    getFrom() {return from;}
    getTo() {return to;}
    getDistance() {return distance;}
};

class GreedyCycle : public TSPSolver
{
public:
    Solution run(const Instance& instance)
    {
        Solution sol;
        const auto& M = instance.M;
        const auto& points = instance.points;
        unsigned int dim = points.size();
        int sumdistance1=0;
        int sumdistance2=0;
        std::vector<Edge> path1;
        std::vector<Edge> path2;

        std::set<int> unVisited;

         for (int i = 0; i < dim; ++i)
        {
            unVisited.insert(i);
        }

        int start1Index = 0;
        int start2Index = getFurthestindex(start1Index, M, dim);
        unVisited.erase(start1Index);
        unVisited.erase(start2Index);

        int bestPoint=getNearestindex(start1Index, M, unVisited);
        int minDistance = M[start1Index][bestPoint];

        Edge e1(start1Index, bestPoint, minDistance);
        path1.push_back(e1);
        sumdistance1+=minDistance;
        unVisited.erase(bestPoint);

        bestPoint=getNearestindex(start2Index, M, unVisited);
        minDistance = M[start2Index][bestPoint];

        Edge e2(start2Index, bestPoint, minDistance);
        path2.push_back(e2);
        sumdistance2+=minDistance;
        unVisited.erase(bestPoint);

        bestPoint=getNearestindex2(path1[0].getFrom(),path1[0].getTo(), M, unVisited);
        e1=Edge(path1[0].getTo(), bestPoint, M[path1[0].getTo()][bestPoint]);
        path1.push_back(e1);
        sumdistance1+=e1.getDistance();
        unVisited.erase(bestPoint);
        e1=Edge(bestPoint, path1[0].getFrom(), M[bestPoint][path1[0].getFrom()]);
        path1.push_back(e1);
        sumdistance1+=e1.getDistance();

        bestPoint=getNearestindex2(path2[0].getFrom(),path2[0].getTo(), M, unVisited);
        e2=Edge(path2[0].getTo(), bestPoint, M[path2[0].getTo()][bestPoint]);
        path2.push_back(e2);
        sumdistance2+=e2.getDistance();
        unVisited.erase(bestPoint);
        e2=Edge(bestPoint, path2[0].getFrom(), M[bestPoint][path2[0].getFrom()]);
        path2.push_back(e2);
        sumdistance2+=e2.getDistance();


        while(!unVisited.empty())
        {
            int bestgain1=std::numeric_limits<int>::max();
            int bestgain2=std::numeric_limits<int>::max();
            int pointtoadd1=-1;
            int edgetoremoveindex1=-1;
            int pointtoadd2=-1;
            int edgetoremoveindex2=-1;


            for(int j=0; j<path1.size(); j++)
            {
                for (int i: unVisited)
                {
                    int gain = M[path1[j].getFrom()][i]+M[path1[j].getTo()][i]-path1[j].getDistance();
                    if (gain < bestgain1)
                    {
                        bestgain1=gain;
                        pointtoadd1=i;
                        edgetoremoveindex1=j;
                    }
                }
            }
            
            path1.insert(path1.begin()+edgetoremoveindex1+1, Edge(path1[edgetoremoveindex1].getFrom(), pointtoadd1, M[path1[edgetoremoveindex1].getFrom()][pointtoadd1]));
            path1.insert(path1.begin()+edgetoremoveindex1+2, Edge(pointtoadd1, path1[edgetoremoveindex1].getTo(), M[pointtoadd1][path1[edgetoremoveindex1].getTo()]));
            path1.erase(path1.begin()+edgetoremoveindex1);
            sumdistance1+=M[path1[edgetoremoveindex1].getFrom()][pointtoadd1]+M[pointtoadd1][path1[edgetoremoveindex1].getTo()];
            unVisited.erase(pointtoadd1);

            
            for(int j=0; j<path2.size(); j++)
            {
                for (int i: unVisited)
                {
                    int gain = M[path2[j].getFrom()][i]+M[path2[j].getTo()][i]-path2[j].getDistance();
                    if (gain < bestgain2)
                    {
                        bestgain2=gain;
                        pointtoadd2=i;
                        edgetoremoveindex2=j;
                    }
                }
            }
            
            path2.insert(path2.begin()+edgetoremoveindex2+1, Edge(path2[edgetoremoveindex2].getFrom(), pointtoadd2, M[path2[edgetoremoveindex2].getFrom()][pointtoadd2]));
            path2.insert(path2.begin()+edgetoremoveindex2+2, Edge(pointtoadd2, path2[edgetoremoveindex2].getTo(), M[pointtoadd2][path2[edgetoremoveindex2].getTo()]));
            path2.erase(path2.begin()+edgetoremoveindex2);
            sumdistance2+=M[path2[edgetoremoveindex2].getFrom()][pointtoadd2]+M[pointtoadd2][path2[edgetoremoveindex2].getTo()];
            unVisited.erase(pointtoadd2);
            
           

        }

        for (Edge e: path1)
        {
            sol.path1.push_back(points[e.getFrom()]);
        }
        for (Edge e: path2)
        {
            sol.path2.push_back(points[e.getFrom()]);
        }

        return sol;

      
    }

private:
     int getFurthestindex(int index, const DistanceMatrix& M, int dim)
    {
        int maxdistance = 0;
        int maxindex = -1;

        for (int i =0; i < dim; i++)
        {
            if (M[index][i] > maxdistance)
            {
                maxdistance = M[index][i];
                maxindex = i;
            }
        }
        return maxindex;
    }
    int getNearestindex(int index, const DistanceMatrix& M, std::set<int> unVisited)
    {
        int mindistance = std::numeric_limits<int>::max();
        int minindex = -1;

        for (int i: unVisited)
        {
            if (M[index][i] < mindistance)
            {
                mindistance = M[index][i];
                minindex = i;
            }
        }
        return minindex;
    }

    int getNearestindex2(int index1, int index2, const DistanceMatrix& M, std::set<int> unVisited)
    {
        int minsum = std::numeric_limits<int>::max();
        int minindex = -1;

        for (int i : unVisited)
        {
            if (M[index1][i]+ M[index2][i] < minsum)
            {
                minsum = M[index1][i]+ M[index2][i];
                minindex = i;
            }

        }
        return minindex;
    }
};

int main(int argc, char* argv[])
{
    std::filesystem::path workDir = "workdir";

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