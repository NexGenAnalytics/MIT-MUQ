#include <random>
#include "readCsv.cpp"

void readCsv(std::string filename, std::vector<double>* a);
void readCsv(std::string filename, std::vector<std::vector<double>>* a);
void readCsv(std::string filename, std::vector<std::vector<double>>* a, int row);
void writeCsv(std::string filename, std::vector<double> a);
void writeCsv(std::string filename, std::vector<std::vector<double>> a);
