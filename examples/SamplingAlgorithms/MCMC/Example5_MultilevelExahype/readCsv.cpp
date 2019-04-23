#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>

void readCsv(std::string filename, std::vector<double>* a){

   std::ifstream infile(filename);
   if(! infile)
       throw std::runtime_error("Could not open file " + filename);
   std::string line;
   int row = 0;
   while (std::getline(infile, line))
   {
     if (line.back() == ';') {
       a->push_back(std::stod(line));
       //std::cout << "number: " << (*a)[row] << std::endl;
       row++;
     }
   }
}

void readCsv(std::string filename, std::vector<std::vector<double>>* a){
    std::cout << "Reading in measurements" << std::endl;
   std::ifstream infile(filename);
   if(! infile)
       throw std::runtime_error("Could not open file " + filename);
   std::string line;
   int row = 0;
   while (std::getline(infile, line))
   {
     if (line.back() == ';') {
         std::stringstream ss(line);
       double d1; ss >> d1;   // read 1
       char cc; ss >> cc;   // read 1
       double d2; ss >> d2;   // read 2
       a->push_back({d1,d2});
       row++;
     }
   }

}


void readCsv(std::string filename, std::vector<std::vector<double>>* a, int row){
    std::ifstream infile(filename);
    if(! infile)
        throw std::runtime_error("Could not open file " + filename);
    std::string line;
    std::getline(infile, line);
    while (std::getline(infile, line))
    {
        std::stringstream ss(line);
        double d1; ss >> d1;   // read 1
        char cc; ss >> cc;   // read 1
        double d2; ss >> d2;   // read 2
        char cc2; ss >> cc2;   // read 1
        double d3; ss >> d3;   // read 2
        char cc3; ss >> cc3;   // read 1
        double d4; ss >> d4;   // read 2
        //std::cout << d1 << " " << d2 << " " << d3 << " " << d4 << std::endl;
        a->push_back({d1,d2,d3,d4});
    }

}

void writeCsv(std::string filename, std::vector<double> a){
    std::ofstream myfile;
    myfile.open (filename);
    for(int i = 0; i < a.size(); i++){
        myfile << a[i] << ";" << std::endl;
    }
    myfile << std::endl;
    myfile.close();
}


void writeCsv(std::string filename, std::vector<std::vector<double>> a){
    std::ofstream myfile;
    myfile.open (filename);
    for(int i = 0; i < a.size(); i++){
        for(int j = 0; j < a[i].size(); j++){
            myfile << a[i][j] << ";";
        }
        myfile <<  std::endl;
    }
    myfile << std::endl;
    myfile.close();
}
