# target link lib to: OpenMP, Sundials, nlopt, hdf5, boost::system, boost::graph
message(STATUS "SetupmuqUtilities")

target_link_libraries(
    muqUtilities PUBLIC
    OpenMP::OpenMP_CXX
    SUNDIALS::cvodes
    SUNDIALS::idas
    SUNDIALS::kinsol
    SUNDIALS::nvecserial
    NLopt::nlopt
    hdf5::hdf5
    hdf5::hdf5_cpp
    hdf5::hdf5_hl
    Boost::system
    Boost::filesystem
    Boost::graph
)
