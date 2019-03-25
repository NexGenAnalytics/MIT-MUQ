#ifndef PARALLELFLAGS_H_
#define PARALLELFLAGS_H_

#if MUQ_HAS_MPI

namespace muq {
  namespace SamplingAlgorithms {

    const int ControlTag = 1;
    const int WorkgroupTag = 2;

    enum ControlFlag : const int {

      // SamplingProblem
      FINALIZE,
      INIT_PROBLEM,
      LOGDENSITY,
      TEST,
      QOI,
      //WRITE_VTK_FROM_QOI,

      // Assignment
      ASSIGN,
      ASSIGN_COLLECTOR,
      UNASSIGN,
      SAMPLE,
      SAMPLE_BOX,
      QUIT,

      SAMPLE_BOX_DONE,

      MEANS,
      MEANS_DONE,

      // Workgroup phonebook
      GET_WORKGROUP,
      SET_WORKGROUP,
      GET_WORKGROUPS,
      UNSET_WORKGROUP,
      GET_LARGEST_INDEX,
      SET_WORKER_READY

    };

  }
}

#endif

#endif
