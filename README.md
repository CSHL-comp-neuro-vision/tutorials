# Hello! Welcome to the tutorials for the CSHL Computational Neuroscience: Vision course. 

## Origins
 
These tutorials have been generously donated by various course participants and organizers over the years. For this reason, you’ll find that one tutorial may not follow exactly the same style/organization as another. If you find one style particularly helpful -- or another unhelpful – let us know! We’re always looking to improve.

## Programming language

Currently, the vast majority of these tutorials are in Matlab, but there are some Python tutorials in the `Python` folder. See the bottom of this page for a summary. 

## Where do I start?

Get your dependencies in order: Some tutorials have shared dependencies, so it is easiest to just go ahead and add the “dependencies_shared” folder to your path, with all its subfolders. To do this, you can navigate to the top level tutorials folder and enter this into your command prompt:

`addpath(genpath('dependencies_shared'))`

Other tutorial-specific dependencies are included in the folder with each tutorial and are loaded by each tutorial script.


## Covering some fundamentals

Some of our tutorials cover quite general concepts that should be broadly useful. These may be a good place to start, even just as a refresh.

### Linear Algebra

Dive into these tutorials to learn why Gilbert Strang has called linear algebra a “wonderful adventure.” Or in a practical sense, maybe you want to gain a better understanding of linear equations, vectors, and matrices. A beginner’s tutorial is offered in “leastsqTutorial.”  If you play with the concepts in this tutorial, you may go from knowing next to nothing about linear algebra to an intuition that will help with multiple linear regression. “linearAlgebraTutorial” offers another take on a general introduction to this material. You could then apply some of these concepts in the “PCATutorial” and “pcasvdTutorial” where you’ll be introduced to principal components analysis and singular value decomposition (which are fundamental concepts for dimensionality reduction, among other applications).

### SignalProcessing

The general principles of how to process and manipulate signals provide important foundations for understanding vision. This collection of signal processing tutorials is a great place to start. Check out the “FourierTutorial” for an introduction to the Fourier Transform (that is, learn how to take a function of time or space (or whatever) and express it as a function of frequency. Dig more into the time domain with “linSysTutorial” and then “Linear_1D_Filter_Tutorial,” which will take you through some applications of filtering signals in the time domain. The “samplingTutorial” will take you through signal sampling in multiple domains. You can move into higher dimensional signals by checking out “imageTutorial,” “imageFormationTutorial,” and “pyramidTutorial.” For a little something different, “ICA_Tutorial” will show you how to use independent components analysis to separate about multiple sources that are mixed together in measured signals.

### Bayesian Estimation

Bayesian statistics provide a principled way to incorporate prior knowledge when you’re using noisy measurements to estimate something (anything?). “bayesTutorialOne” and “bayesTutorialTwo” offer two examples of applying Bayesian estimation techniques. The Bayesian framework also provides an excellent model of perceptual inference!

### Diffusion Process

Diffusion processes provide a flexible way to model dynamic systems that are subjected to random fluctuations over time. This tutorial specifically uses a diffusion process to model perceptual decisions that are made by accumulating (noisy) information about a stimulus.

### Modeling Individual Neurons

Here we have an assortment of tutorials that consider the properties and processes that characterize individual neurons, or “nerve cells” if you want to be a bit old fashioned. You might want to start with the “Hodgkin_HuxleyTutorial” and “intAndFireTutorial,” which implement models of how electrical properties of neurons generate action potentials. But real action potentials have important stochastic properties, too. Check out “poissonTutorial” for an introduction to the Poisson model of stochastic neuronal firing. “stochasticProcessesTutorial” introduces a few facts about renewal processes, beginning with the poisson point process and extending intuitions to gamma processes and beyond.

### Or you want to model a lot of neurons…Functional MRI (fMRI)

Functional MRI (fMRI) measures changes in blood flow and oxygenation associated with the underlying neuronal response. Here we have two tutorials on fMRI. Curious about event-related fMRI study designs and deconvolution-based analyses? “ER_fMRI_Tutorial” is a great place to start. You can build on what you learned in that tutorial by trying out “fMRI_Classification_Tutorial,” which demonstrates how to run a classification analysies on multivariate fMRI data. Bonus: the techniques described in this tutorial can be applied nearly identically to populations of neural responses.

## What about putting it all together?

### Modeling Light Detection

Light detection is an awfully important part of vision. The “PhotonDetectionTutorial” explores the process of light detection, and along the way reinforces concepts from Poisson distributions, signal detection theory, Bayesian estimation and more. Meanwhile, “nonlinearPoolingDemo” follows a similar framework, focusing on the role of non-linearities in the retina. A third tutorial engages you in a similar thought experiment (“sdtTutorial”) while also exploring a discrimination task. And for something a little different (that is, a deep-dive into phototransduction) check out “DiffEqTutorial.”

### Choice Probability

This tutorial deals with the relationship between single neuron activity and the behavioral choices that we (or monkeys) make in psychophysics. It incorporates topics from signal detection theory and motion processing. If you’ve covered these topics, this tutorial should provide you with some good pay off because it shows how signal detection theory gives us something to measure in neurons that will ultimately constrain theories of neural coding.

### Generalized Linear Model (GLM) (coming soon)

### Psychophysics

Is it the best type of physics? Decide for yourself with this suite of 3 tutorials on psychophysics that take you through applications of signal detection theory (“Psycho_Tutorial_I_SDT”), estimation of psychophysical thresholds (“Psycho_Tutorial_II_Thresholds”), and bootstrapping techniques (“Psycho_Tutorial_III_Bootstrapping”).

### All About Color

Color – it’s just in your head…or is it? This is a suite of 3 tutorials on color. “colorRenderingTutorial” will introduce you to modelling surfaces and illuminants and how the photoreceptors respond to the resulting light that enters the eye. “colorSpaceTutorial” builds on this to illustrate aspects of the color spaces that we use to represent our perception of color. Finally, the “colorConstancy” tutorial delves in to a simple model of color constancy – that is, how we perceive objects to have consistent colors even under changing illumination. These tutorials are best consumed in the order laid out here (rendering -> space -> constancy)

### Modeling Attention

This tutorial works through the 'normalization model of attention' based on the 2009 paper published in Neuron by John Reynolds and David Heeger. The model is worth going through because it hits on a range of topics associated with this course - linear filtering, normalization, population responses, and – of course – attention.  It's a powerful model because it can predict a range of physiological, behavioral and fMRI results from attention studies.

### Modeling Motion Processing

Much has been learned about vision by studying the visual signals and visual processing associated with motion. Consider starting with “motionTutorial”: yhis tutorial presents some concepts for representing and analyzing visual motion, including motion energy and computer vision. Dig in more deeply to a model of how motion is processed in MT (Rust at al., 2006) in “RustMTModel.” For another take, “STCovTutorial” demonstrates how spike-triggered covariance can be used to recover functional models in white-noise experiments.

## Python tutorials include introductions to:

* Machine learning (SVM)
* FFT
* Signal filtering (low/high pass filtering)
* Mutual Information

All python tutorials are standalone jupyter notebooks (.ipynb files, which can also be uploaded and run using google colab if you don't have a standalone install on your machine). 

## TODO: 

Organization: Not convinced this mixture of shared/tutorial-specific dependencies makes sense, but that's what I started with. If we stick with this organization, need to review the comments at the top of each script and make sure it matches the current actual dependencies

Tutorials currently giving errors in Matlab 2020:
* pcasvdTutorial
* FourierTutorial
* pyramidTutorial
* stochasticProcessesTutorial
* fMRI_Classification_Tutorial
* choiceProbabiltyTutorial
* Psycho_Tutorial_III_Bootstrapping

GLM Tutorial: Pillow needs to look at this and update the guide to describe what it is
