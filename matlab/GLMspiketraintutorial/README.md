# GLMspiketraintutorial
Simple tutorial on Gaussian and Poisson generalized linear models (GLMs) for spike train data.  
**Author**: [Jonathan Pillow](http://pillowlab.princeton.edu), Nov 2016.

**NEW (Feb 2022)**:  There is now a [python version](https://github.com/pillowlab/GLMspiketraintutorial_python) of this tutorial!

**Slides**: This tutorial was prepared for use in a
"Short Course" on [Data Science and Data Skills for Neuroscientists](https://neuronline.sfn.org/scientific-research/data-science-and-data-skills-for-neuroscientists#:~:text=Data%20science%20is%20fast%2Dgrowing,be%20used%20in%20different%20circumstances) organized at the SFN 2016 meeting. The slides used during the 1-hour
presentation are available [here](https://github.com/pillowlab/GLMspiketraintutorial/blob/master/slides/slides_SFNshortcourse_Nov2016.pdf). 

**Dataset**:  A small dataset required for the tutorial is available
[here](https://pillowlab.princeton.edu/data/data_RGCs.zip). (Simply download this dataset, unzip it, and place it in the directory where you are running the tutorial).

**Description**: This tutorial contains four scripts that introduce various methods for fitting and simulating from Gaussian and Poisson regression models for spike train data. Each script is self-contained, with 'blocks' of code that can be executed interactively, which demonstrate each step in the fitting / analysis / model comparison pipeline:

* **tutorial1_PoissonGLM.m** - illustrates the fitting of a
linear-Gaussian GLM (also known as the 'linear least-squares
regression model') and a Poisson GLM (aka 'linear-nonlinear-Poisson'
model) to single retinal ganglion cell responses to a temporal white
noise stimulus.

* **tutorial2_spikehistcoupledGLM.m** - fitting of an autoregressive
Poisson GLM (i.e., a GLM with spike-history) and a multivariate
autoregressive Poisson GLM (a GLM with spike-history AND coupling
between neurons).

* **tutorial3_regularization_linGauss.m** - regularizing
  linear-Gaussian model  parameters using maximum a posteriori (MAP)
  estimation under two kinds of priors:
  - (1) ridge regression (aka  "L2 penalty"); 
  - (2) L2 smoothing prior (aka "graph Laplacian").  


* **tutorial4_regularization_PoissonGLM.m** - MAP estimation of
  Poisson-GLM parameters using same two priors considered in
  tutorial3.

<br>
 
------- 

**Citation**:  If you wish to cite this tutorial, feel free to acknowledge the paper from which it developed:  
[Pillow et al, *Nature* 2008](http://pillowlab.princeton.edu/pubs/abs_Pillow08_nature.html).
 



 <br>

-------

**Relevance / comparison to other GLM packages**:

This tutorial is designed primarily for pedagogical purposes. The
tutorial scripts are (almost entirely) self-contained, making it easy
to understand the basic steps involved in simulating and fitting. It
is easy to alter these scripts (e.g., to incorporate different kinds
of regressors, or different kinds of priors for
regularization). However, this implementation is not memory-efficient
and does not support some of the advanced features available in other
GLM packages (e.g., smooth basis functions for spike-history filters,
memory-efficient temporal convolutions, different timescales for
stimulus and spike-history components, low-rank parametrization of
spatio-temporal filters, flexible handling of trial-based data).  For
more advanced features and applications, see the following two
repositories:

- [neuroGLM](http://pillowlab.princeton.edu/code_neuroGLM.html) -
  designed for single-neuron, trial-structured data. Supports flexible design matrices with multiple types of
  regressors. Relevant pub: [Park et al, *Nat Neurosci* 2014](http://pillowlab.princeton.edu/pubs/abs_ParkI_NN14.html).

- [GLMspiketools](http://pillowlab.princeton.edu/code_GLM.html) -
  designed for single- and multi-neuron spike trains with flexible
  nonlinearities, multiple timescales, and low-rank parametrization of
  filters.  Relevant pub: [Pillow et al, *Nature* 2008](http://pillowlab.princeton.edu/pubs/abs_Pillow08_nature.html).

