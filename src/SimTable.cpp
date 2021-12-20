#include "DistTable.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include "math.h"

SimTable::SimTable(string file_path): file_path(file_path) {}

int DistTable::initTable(){
    std::ifstream fin (this->file_path);
    
    if (!fin.is_open()) {
        std::cerr << "error: file open failed '" << this->file_path <<"'.\n";
        return 1;
    }

    int n,m;
    fin >> n >> m;

    index.resize(n);
    columns.resize(m);

    std::vector<std::vector<float>>   data(n,std::vector<float>(m));

    for (size_t i = 0; i < m; i++) {
        fin >> columns[i];
    }
    

    for (size_t i = 0; i < n; i++) {
        fin >> index[i];
        for (size_t j = 0; j < m; j++){
            fin >> data[i][j];
        }
    }

    table = data;

    return 0;
}


float DistTable::interpolate(float x,float y) {
    //Thanks to Noa, no thanks:
    float rowNum = lower_bound(index.begin(), index.end(), x) - index.begin() - 0.5;//abs(x-index[0])/abs(index[1] - index[0]);
    
    if(rowNum > index.size()-1){
        rowNum = rowNum-1;
    }

    vector<float> lowerRow = table[floor(rowNum)];
    vector<float> upperRow = table[ceil(rowNum)];

    vector<float> newRow(lowerRow.size());

    float lowerFactorRow = abs(x-index[floor(rowNum)])/abs(index[ceil(rowNum)] - index[floor(rowNum)]);
    float upperFactorRow = abs(index[ceil(rowNum)]-x)/abs(index[ceil(rowNum)] - index[floor(rowNum)]);

    for (size_t i = 0; i < lowerRow.size(); i++) {
        newRow[i] = lowerFactorRow*lowerRow[i] + upperFactorRow*upperRow[i];
    }
    
    float colNum = lower_bound(columns.begin(), columns.end(), y) - columns.begin() - 0.5;//abs(y-columns[0])/abs(columns[1] - columns[0]);

    float lowerFactorCol = abs(y-columns[floor(colNum)])/abs(columns[ceil(colNum)] - columns[floor(colNum)]);
    float upperFactorCol = abs(columns[ceil(colNum)]-y)/abs(columns[ceil(colNum)] - columns[floor(colNum)]);

    float result = lowerFactorCol*newRow[floor(colNum)] + upperFactorCol*newRow[ceil(colNum)];

    return result;
}

// DistTable DistTable::interpolateOnRows() {


// }

// DistTable DistTable::interpolateonColumns() {

// }


void DistTable::displayTable(){
    for (size_t i = 0; i < index.size(); i++) {
        for (size_t j = 0; j < columns.size(); j++){
            cout << fixed << setprecision(3) << table[i][j] << " ";
        }
        cout << endl;
    }
}


void DistTable::displayHeader() {
    for (size_t i = 0; i < columns.size(); i++){
            cout << fixed << setprecision(3) << columns[i] << " ";
        }
    cout << endl;
}
void DistTable::displayIndex() {
    for (size_t i = 0; i < index.size(); i++){
            cout << fixed << setprecision(3) << index[i] << " ";
        }
    cout << endl;
}


DistTable::~DistTable() {
}

// Thanks to my lovely wife which i adore so much