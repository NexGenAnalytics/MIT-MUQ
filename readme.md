# Overview

Welcome to MUQ (pronounced “muck”), a modular software framework for defining and solving forward and inverse uncertainty quantification problems.

## Purpose

Uncertainty quantification (UQ) is important in many different applications.
MUQ aims to make advanced probabilistic UQ tools easy to use in either c++ or python,
and enable cutting-edge method development through its modular structure.

MUQ has a variety of capabilities, including:

*  Various Markov chain Monte Carlo methods
*  Graphical modeling with a mix of statistical and physical components.
*  Gaussian processes
*  Karhunen Loève expansions.
*  Transport maps
*  Nonlinear Optimization
*  Generalized Polynomial Chaos Expansions

## Installation:

MUQ is available on Linux and OSX as a conda package, docker image, or from source. For many users, getting started can be as easy as running

```
conda install -c conda-forge muq
```

For more installation options, check out the [installation guide](https://nexgenanalytics.github.io/MIT-MUQ/latest/muqinstall.html).

## Getting Started

MUQ is composed of several different modules, which work together to define and solve UQ problems. 
Documentation for each of these modules is included with our doxygen-generated [API documentation](https://nexgenanalytics.github.io/MIT-MUQ/latest/index.html). 
Most applications will require using the [modeling module](https://nexgenanalytics.github.io/MIT-MUQ/latest/group__modeling.html) to define statistical models or interact with user-defined models. Learning the basics of this module is therefore a good place to start.

#### Interested in forward UQ?

- First, get acquainted with the [modeling module](https://nexgenanalytics.github.io/MIT-MUQ/latest/group__modeling.html).  You'll need to use one or more instances of the [ModPiece class](https://nexgenanalytics.github.io/MIT-MUQ/latest/classmuq_1_1Modeling_1_1ModPiece.html) to define the model that will be evaluated by the UQ algorithm.
- Once you have a model, check out the [polynomial chaos module](https://nexgenanalytics.github.io/MIT-MUQ/latest/group__polychaos.html).
- Other examples can be found by selecting the "PCE" examples on the MUQ [webpage](https://nexgenanalytics.github.io/MIT-MUQ/examples.html).

#### Want to tackle Bayesian inverse problems?

- Just like for forward UQ, you'll want to get familiar with the [modeling module](https://nexgenanalytics.github.io/MIT-MUQ/latest/group__modeling.html) module to define a forward model.  The [WorkGraph class](https://nexgenanalytics.github.io/MIT-MUQ/latest/classmuq_1_1Modeling_1_1WorkGraph.html) within the modeling module is also used to combine multiple components (e.g., the prior, forward model, and likelihood function) comprising the Bayesian posterior distribution.
- Look at methods in the [sampling algorithms](https://nexgenanalytics.github.io/MIT-MUQ/latest/group__sampling.html) module to generate samples of your Bayesian posterior.
- Other examples can be found by filtering the "MCMC" examples on the MUQ [webpage](https://nexgenanalytics.github.io/MIT-MUQ/examples.html).

You can also find many [examples](https://nexgenanalytics.github.io/MIT-MUQ/examples.html) using both the c++ and Python interfaces to MUQ.  These examples can provide useful starting places for using MUQ on your own problems.

#### Getting Connected

Join the MUQ Slack channel via our [website](https://nexgenanalytics.github.io/MIT-MUQ) to get in touch with MUQ developers and other users. We are always happy to help!

## Citing

When publishing work based on MUQ, please cite our publication in the Journal of Open Source Software.

<div><pre><code class="language-plaintext">@article{Parno2021,
  doi = {10.21105/joss.03076},
  url = {https://doi.org/10.21105/joss.03076},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {68},
  pages = {3076},
  author = {Matthew Parno and Andrew Davis and Linus Seelinger},
  title = {MUQ: The MIT Uncertainty Quantification Library},
  journal = {Journal of Open Source Software}
}</code></pre></div>

## Contributing

#### Want to help develop MUQ?

Yes, please! We frequently discuss future developments on Slack ([join via our website](https://nexgenanalytics.github.io/MIT-MUQ)), so feel free to drop by!
Then fork the [muq2 repository](https://github.com/NexGenAnalytics/MIT-MUQ) and submit a pull request when ready.
Also check out our [style guide](https://nexgenanalytics.github.io/MIT-MUQ/latest/muqstyle.html).

#### Find a bug?

[Submit the issue](https://github.com/NexGenAnalytics/MIT-MUQ/issues).  Make sure to label the issue as a bug.

#### Want a new feature?

[Submit a request](https://github.com/NexGenAnalytics/MIT-MUQ/pulls).  Label the issue as an enhancement or proposal.


#### Developer Information
- \subpage infrastructure
- \subpage muqstyle
