#pragma once

#include <vector>
#include <string>

using namespace std;

class DistTable
{
private:
    vector<vector<float>> table;
    string file_path;
    vector<float> index;
    vector<float> columns;
    size_t numRows;
    size_t numCols;

public:
    DistTable(string file_path);
    ~DistTable();
    int initTable();
    float interpolate(float x,float y);
    // DistTable interpolateOnRows();
    // DistTable interpolateonColumns();
    void displayTable();
    void displayHeader();
    void displayIndex();


};
