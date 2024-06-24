### Tutorial for basic continuous time recurrent neural network using pytorch. Train models to perform tasks that similar to those used in many neuroscience experiments

* Go-NoGo task
* Delayed Match To Sample (DMTS)
* The context dependent task from [Mante et al., 2013](https://www.nature.com/articles/nature12742)

#### If you already have a full python install from Anaconda, then you can probably stay in your (base) env and may just have to install `torch` (see [here](https://pytorch.org/get-started/locally/) for more details and for information about how to install on your system). Make sure you check that you're installing CUDA specific dependencies if you have a CUDA capable GPU. On a Mac running a M1/M2 chip you can use MPS/Metal (see [here](https://developer.apple.com/metal/pytorch/))

#### If you want to create a custom environment that has minimal other stuff in it, you can do something like this (only tested on a Mac running Sonoma). 

Open a terminal (command window). If you are in the `(base)` env from conda, you can deactivate e.g.:

`conda deactivate` 

Install virtualenv

`pip install virtualenv`

Upate pip for good luck

`python -m pip install --upgrade pip`

Create a virtualenv (calling mine pt for pytorch, but you can pick whatever name is useful for you). 

`python -m venv ~/pt`

Activate the venv

`source ~/pt/bin/activate`

Install pytorch and a few other packages in the venv. We're primarily using torch, but you'll want jupyter notebooks and matplotlib as well

`pip install torch`

`pip install jupyter`

`pip install matplotlib`
