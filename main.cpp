#include <iostream>
#include <iostream>
#include <limits.h>
#include <string.h>
#include <queue>
#include <vector>
#include <fstream>
#include <windows.h>

using namespace std;

void output(LONGLONG time, string algo_name, int max_flow)
{
    ofstream out("../output.txt", std::ios::app);
    if(out.is_open())
    {
        out << "Algorithm's name: " << algo_name << ". Time consumed: " << time << " ns." << " Max flow is " << max_flow << '\n';
    }
    out.close();
}

void output(string& filename)
{
    ofstream out("../output.txt", std::ios::app);
    if(out.is_open())
    {
        out << "File name: " << filename << '\n';
    }
    out.close();
}

void output(string text, bool append)
{
    ofstream out;
    if(append)
        out.open("../output.txt", std::ios::app);
    else
        out.open("../output.txt");
    if(out.is_open())
        out << text << '\n';
    out.close();
}

int readInput(string& input, vector<vector<int>>& graph)
{
    int V = 0;
    string line;

    ifstream in(input);

    if(in.is_open())
    {
        while(getline(in, line))
        {
            graph.emplace_back(vector<int>());
            int i = 0;
            string tmp;
            while(i < line.size())
            {
                if(line[i] != '\t')
                    tmp += line[i];
                else
                {
                    graph[V].push_back(stoi(tmp));
                    i++;
                    tmp = "";
                    continue;
                }
                i++;
            }
            graph[V].push_back(stoi(tmp));
            V++;
        }
        in.close();
    }
    return V;
}

bool dinic_bfs(vector<vector<int>>& graph, vector<vector<int>>& f, vector<int>& shortest_path, size_t V, int s, int t)
{
    for (int i = 0; i < shortest_path.size(); ++i)
    {
        shortest_path[i] = INT_MAX;
    }
    shortest_path[s] = 0;
    queue<int> q;
    q.push(s);

    while(!q.empty())
    {
        int u = q.front();
        q.pop();
        for (int i = 0; i < V; ++i)
        {
            if(graph[u][i] > 0 && (f[u][i] < graph[u][i]) && shortest_path[i] == INT_MAX)
            {
                shortest_path[i] = shortest_path[u] + 1;
                q.push(i);
            }
        }
    }
    return shortest_path[t] != INT_MAX;
}

int dinic_dfs(int u, int minC, vector<vector<int>>& graph, vector<vector<int>>& f, vector<int>& shortest_path, vector<int>& p, int t, size_t V)
{
    if(u == t || minC == 0)
        return minC;
    for (int i = p[u]; i < V; ++i)
    {
        if(shortest_path[i] == shortest_path[u] + 1)
        {
            int flow = dinic_dfs(i, min(minC, graph[u][i] - f[u][i]), graph, f, shortest_path, p, t, V);
            if(flow != 0)
            {
                f[u][i] += flow;
                f[i][u] -= flow;
                return flow;
            }
            p[u]++;
        }
    }
    return 0;
}

int dinic(size_t V, vector<vector<int>>& graph, int s, int t)
{
    int maxFlow = 0;
    vector<int> shortest_path(V);
    vector<vector<int>> f(V);
    for (int i = 0; i < f.size(); ++i)
    {
        f[i] = vector<int>(V);
    }
    vector<int> p(V);

    while(dinic_bfs(graph, f, shortest_path, V, s, t))
    {
        for (int i = 0; i < p.size(); ++i)
        {
            p[i] = 0;
        }
        int flow = dinic_dfs(s, INT_MAX, graph, f, shortest_path, p, t, V);
        while(flow != 0)
        {
            maxFlow += flow;
            flow = dinic_dfs(s, INT_MAX, graph, f, shortest_path, p, t, V);
        }
    }
    return maxFlow;
}

bool find_path(size_t V, vector<vector<int>>& rGraph, int cur, int t, vector<int>& parent, vector<bool>& visited)
{
    visited[cur] = true;

    for (int i = 0; i < V; ++i)
    {
        if(!visited[i] && rGraph[cur][i] > 0)
        {
            parent[i] = cur;
            bool found = find_path(V, rGraph, i, t, parent, visited);
            if(found)
                return true;
        }
    }

    return visited[t];
}

int fordFulkerson(size_t V, vector<vector<int>>& graph, int s, int t)
{
    int u, v;

    vector<vector<int>> rGraph = graph;

    vector<int> parent(V);
    parent[s] = -1;
    int max_flow = 0;

    vector<bool> visited;
    visited.assign(V, 0);

    while (find_path(V, rGraph, s, t, parent, visited))
    {
        int path_flow = INT_MAX;
        for (v=t; v!=s; v=parent[v])
        {
            u = parent[v];
            path_flow = min(path_flow, rGraph[u][v]);
        }

        for (v=t; v != s; v=parent[v])
        {
            u = parent[v];
            rGraph[u][v] -= path_flow;
            rGraph[v][u] += path_flow;
        }
        max_flow += path_flow;
        for (int i = 0; i < V; ++i)
        {
            visited[i] = false;
        }
    }

    return max_flow;
}


// for Edmonds-Karp
bool bfs(size_t V, vector<vector<int>>& rGraph, int s, int t, int parent[])
{
    vector<bool> visited;
    visited.assign(V, false);

    queue <int> q;
    q.push(s);
    visited[s] = true;
    parent[s] = -1;

    while (!q.empty())
    {
        int u = q.front();
        q.pop();

        for (int v=0; v<V; v++)
        {
            if (!visited[v] && rGraph[u][v] > 0)
            {
                q.push(v);
                parent[v] = u;
                visited[v] = true;
            }
        }
    }

    return visited[t];
}


int edmondsKarp(size_t V, vector<vector<int>>& graph, int s, int t)
{
    int u, v;

    vector<vector<int>> rGraph = graph;

    int parent[V];

    int max_flow = 0;

    while (bfs(V, rGraph, s, t, parent))
    {
        int path_flow = INT_MAX;
        for (v=t; v!=s; v=parent[v])
        {
            u = parent[v];
            path_flow = min(path_flow, rGraph[u][v]);
        }

        for (v=t; v != s; v=parent[v])
        {
            u = parent[v];
            rGraph[u][v] -= path_flow;
            rGraph[v][u] += path_flow;
        }
        max_flow += path_flow;
    }

    return max_flow;
}

int main() {
    string ending[] = {"_0.0.txt", "_0.5.txt", "_1.0.txt", "_disco.txt"};
    output("", false);
    for(int amount = 10; amount <= 310; amount+=300)
    {
        for(int k = 0; k < 4; k++)
        {
            output("\n", true);
            vector<vector<int>> graph;
            string input = "../input/input_";
            input += to_string(amount) + ending[k];
            int V = readInput(input, graph);
            // find s and t
            vector<int> s, t;
            for (int i = 0; i < graph.size(); ++i)
            {
                int sum_s = 0, sum_t = 0;
                for (int j = 0; j < graph[i].size(); ++j)
                {
                    sum_s += graph[j][i];
                    sum_t += graph[i][j];
                }
                if (sum_t == 0)
                    t.push_back(i);
                if (sum_s == 0)
                    s.push_back(i);
            }
            int max_flow = 0;
            LARGE_INTEGER beg, end, freq, delta;
            QueryPerformanceFrequency(&freq);
            for (int j = 0; j < 10; j++)
            {
                max_flow = 0;
                QueryPerformanceCounter(&beg);
                for (int i = 0; i < s.size(); ++i)
                {
                    max_flow += fordFulkerson(V, graph, s[i], t[i]);
                }
                QueryPerformanceCounter(&end);
                delta.QuadPart += ((end.QuadPart - beg.QuadPart) * 1000000000) / freq.QuadPart;
            }
            delta.QuadPart /= 10;
            output(input);
            output(delta.QuadPart, "Ford Fulkerson", max_flow);
            delta.QuadPart = 0;

            for (int j = 0; j < 10; j++)
            {
                max_flow = 0;
                QueryPerformanceCounter(&beg);
                for (int i = 0; i < s.size(); ++i)
                {
                    max_flow += edmondsKarp(V, graph, s[i], t[i]);
                }
                QueryPerformanceCounter(&end);
                delta.QuadPart += ((end.QuadPart - beg.QuadPart) * 1000000000) / freq.QuadPart;
            }
            delta.QuadPart /= 10;
            output(delta.QuadPart, "Edmons Karp", max_flow);
            delta.QuadPart = 0;

            for (int j = 0; j < 10; j++)
            {
                max_flow = 0;
                QueryPerformanceCounter(&beg);
                for (int i = 0; i < s.size(); ++i)
                {
                    max_flow += dinic(V, graph, s[i], t[i]);
                }
                QueryPerformanceCounter(&end);
                delta.QuadPart += ((end.QuadPart - beg.QuadPart) * 1000000000) / freq.QuadPart;
            }
            delta.QuadPart /= 10;
            output(delta.QuadPart, "Dinic", max_flow);
        }
    }
    return 0;
}