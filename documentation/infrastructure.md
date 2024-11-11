\page infrastructure Developer Infrastructure Notes

### Git
- Hosted at https://github.com/NexGenAnalytics/MIT-MUQ

### Testing
- GTest for running c++ unit tests for each compile groups, sequential and MPI parallel
- RunAllNotebooks.sh for automatically running all python notebook examples

### CI
- See https://github.com/NexGenAnalytics/MIT-MUQ/tree/master/.github/workflows in MUQ2 repo

### Website
- Hosted on Github pages.
- Accessible at https://nexgenanalytics.github.io/MIT-MUQ/
- When something is committed to this repository, Github pipelines runs Jekyll and generate the html.
- The site is static, but the slack invitation button requires extra processing that cannot occur on a static site.  To overcome this,  a "cloud function" has been set up on the Google cloud platform to handle POST requests when the invitation modal and interact with the slack API.

### Doxygen
- Built via CI when release tag is defined on master repo, pushed as commit to website repo

### Docker (in progress)
- Dockerfiles are located at https://github.com/NexGenAnalytics/MIT-MUQ-containers
- The docker images are built on dockerhub automatically via pipelines.
- Images are available at: https://github.com/orgs/NexGenAnalytics/packages?repo_name=MIT-MUQ-containers

### Conda
- The MUQ2 conda image lives on conda-forge.
- The recipe for the MUQ2 conda image lives in the [muq-feedstock repository](https://github.com/conda-forge/muq-feedstock) on github.   
- Directions for updating the recipe manually can be found [here](https://conda-forge.org/docs/maintainer/updating_pkgs.html).
- Linus and Matt are currently recipe maintainers.