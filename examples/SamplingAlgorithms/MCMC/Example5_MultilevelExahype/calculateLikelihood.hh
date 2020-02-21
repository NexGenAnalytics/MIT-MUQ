#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>

//Read in Output values
//compare to "real" values
//return number


double calculateLikelihood(std::vector<double> solution){
    int numProbes = 16;
    double likelihood =  0.0;  //0.25*(ll[0] + ll[1] + ll[2] + ll[3]);

    //calculate differences at the 4 different probe points
    for(int i_probe=1; i_probe <= numProbes; i_probe++){

        double diff;
        double dh;

        {
        std::string filestr ="Output/Reference/waveheight"+std::to_string(i_probe)+"-rank-0.probe";
        const char* filename = filestr.c_str();
        std::ifstream infile(filename);
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
                diff = d4; 
        }
        }
        /*{
        std::string filestr ="Output/waveheight"+std::to_string(i_probe)+"-rank-0.probe";
        const char* filename = filestr.c_str();
        std::ifstream infile(filename);
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
                dh = d4;
        }
        std::cout << "Probe_ref " << i_probe << " " << diff << std::endl;
        diff -= dh;
        std::cout << "Probe current " << i_probe << " " << dh << std::endl;
        std::cout << "Probe diff " << i_probe << " " << diff << std::endl;
        }*/
        std::cout << "Solution " << i_probe - 1 << " "  << solution[i_probe-1] << std::endl;
        diff-=solution[i_probe-1];

        std::cout << "Diff     " << i_probe - 1 << " " << diff << std::endl;
        
        /*std::string arg = "python Output/TFMisfit/TFMisfit.py -r Output/Reference/waveheight"+std::to_string(i_probe+1)+".probe -s Output/waveheight"+std::to_string(i_probe+1)+"-rank-0.probe > Output/output.csv 2> Output/err";
        const char *command = arg.c_str(); 
        system(command);

        std::string filename = "Output/output.csv";
        std::ifstream infile(filename);
        if(! infile)
            throw std::runtime_error("Could not open file " + filename);
        std::string line;
        std::getline(infile, line);
        std::stringstream ss(line); //read bathymetry misfit
        double d1; ss >> d1;   // read 1
        char cc; ss >> cc;   // read 1
        double d2; ss >> d2;   // read 2
        std::cout << d1 << " " << d2 << std::endl;
        likelihood = std::max(likelihood, d1*d1 + d2*d2);*/  //2 norm
        likelihood -= 1.0/numProbes * .5 * diff*diff * 500;
    }
    std::ofstream ost;
    ost.open("likelihood.log", std::ios::app);
    ost << std::exp(likelihood) << std::endl;
    return likelihood;
}
