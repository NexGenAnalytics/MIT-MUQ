# First, let's update our package manager and install some basic dependencies
apt update
apt install -y git cmake g++

# Optionally, install some advanced dependencies. muq will automatically download these if not found on your system.
apt install -y libhdf5-dev libboost-all-dev libeigen3-dev libsundials-dev libnlopt-dev


# Get muq from git repository
git clone https://bitbucket.org/mituq/muq2.git


# Let's compile! This is just the usual cmake build procedure.
cd muq2/build

# Here you can choose your install directory by setting CMAKE_INSTALL_PREFIX
cmake -DCMAKE_INSTALL_PREFIX=$PWD/install ..

muq_install_dir=$PWD/install # Remember MUQ install directory for later

make -j4

make install


# Now let's try out our shiny new MUQ install with an example!

cd ../examples/SamplingAlgorithms/MCMC/Example1_Gaussian/cpp

# Create a build directory
mkdir -p build

cd build

# Configure the example with cmake, pointing to our MUQ install
cmake -DMUQ_DIR=${muq_install_dir}/lib/cmake/MUQ ..

make

./GaussianSampling
