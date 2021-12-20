#pragma once

#include <vector>
#include <string>

using namespace std;

class SimTable
{
private:
    vector<vector<float>> table;
    string file_path;
    vector<float> index;
    vector<float> columns;
    size_t numRows;
    size_t numCols;

public:
    SimTable(string file_path);
    ~SimTable();
    int initTable();

    // DistTable interpolateOnRows();
    // DistTable interpolateonColumns();
    void displayTable();
    void displayHeader();
    void displayIndex();


};
