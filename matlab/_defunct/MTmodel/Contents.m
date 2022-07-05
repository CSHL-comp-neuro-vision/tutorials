% ---------------------------------------------------------------------
% MTmodel - implementation of a 2-stage model of neuronal responses
% in primate visual areas V1 and MT.
%
% Version 1.0,  September 2005.
% Authors: Timothy Saint (saint@cns.nyu.edu) and Eero Simoncelli (eero.simoncelli@nyu.edu)
%
% See README file for brief description.
% See ChangeLog file for latest modifications. 
% See HELP subdirectory for demonstrations.
% Type "help <command-name>" for documentation on individual commands.
% ----------------------------------------------------------------------
%
% ----------------- HELP functions ---------------------------------------
%
% shTutorial1.m         Overview of the software, generating tuning curves, etc.
%
% ----------------- MODEL functions ---------------------------------------
%
% mt2sin                    Get the paramters for gratings preferred by given MT neurons.
% shGetDims                 Find the size your stimulus must be to successfully run the MT model.
% shGetNeuron               Extract the response of neuron(s) at one spatial position from shModel outputs.
% shGetScale                Extract all the neuronal responses at a particular scale from shMatrix.
% shGetSubPop               Extract the responses of all the neurons with similar tuning from shModel output
% shMkV1Filter              Make the linear filter that is the front end of a given model V1 neuron.
% shModel                   Run the Simoncelli & Heeger model
% shMtPopulationResponse    Compute the response of a large population of MT neurons to a stimulus.
% shV1PopulationResponse    Compute the response of a large population of V1 neurons to a stimulus.
% v12sin                    Get the paramters of the drifting grating preferred by given V1 neurons.
%
% ------------------ STIM functions -------------------------------------
%
% mkBar             Make a drifting bar stimulus
% mkDots            make a drifting dot stimulus
% mkFract           make a drifting fractal noise stimulus 
% mkGlass           Make a translational dynamic glass pattern movie
% mkPlaid           make a plaid stimulus
% mkSin             make a drifting grating
% mkWedge           make a wedge in the Fourier domain
% mkWin             make a circular window with raised cosine edges
%
% ------------------ TUNE functions -------------------------------------
%
% shTuneBarSpeed                Response vs. speed of a drifting bar 
% shTuneDotCoherence            Response vs. coherence of random dot motion
% shTuneDotDensity              Response vs. dot density
% shTuneDotDirection            Response vs. direction of dot motion
% shTuneDotMaskDirection        Response vs. direction of a dot mask
% shTuneDotSpeed                Response vs. speed of dot motion
% shTuneGratingArea             Response vs. size of drifting grating's window
% shTuneGratingContrast         Response vs. contrast of drifting grating
% shTuneGratingDirection        Response vs. motion of grating drift
% shTuneGratingMaskDir          Response vs. direction of mask grating
% shTuneGratingSf               Response vs. spatial frequency of grating
% shTuneGratingTf               Response vs. temporal frequency of grating
% shTunePlaidDirection          Response vs. direction of plaid motion
% shTunePlaidOriOri             Response vs. orientation of both plaid components
%
% ----------------- PARS functions ---------------------------------------
%
% mkGaussianFilter          Make a 1D gaussian filter with a given SIGMA.
% shPars                    Get a default PARS structure
% shParsScaleFactors        Set scale factors for the pars structure and pick pars.mtalpha
% shParsV1PopulationDirections  Get evenly spread V1 neurons for a population
%
% ------------------ SHOW functions --------------------------------------
% 
% shShowMtPopulationResponse            Show the response of a population of MT neurons          
% shShowV1NeuronSpectrum                Show the fourier spectrum of V1 neuron(s)
% shShowV1PopulationDirectionsDots      Show the neurons in V1 population as dots on a sphere
% shShowV1PopulationDirectionsDotsMovie Show the neurons in V1 population being chosen
% shShowV1PopulationResponse            Show the response of a population of V1 neurons
% 
% ----------------- Bonus functions ---------------------------------------
%
% atan3                         Four quadrant inverse tangent with 0 <= theta <= 2.*pi
% crossAxes                     Put good old fashioned axes through the center of a plot.
% cyl2rec                       Transform [az, h, r] coordinates to [y, x,
% draw3dLevelSurface            Draw a 3D level surface of a 3D function.
% drawCylinder                  Draw a cylinder on the current axes.
% drawPlane                     Draw a plane on the current axes.
% drawSphere                    Draw a sphere on the current axes.
% flipBook                      Show the values of a 3D matrix as a movie
% gray2rgbsc                    Convert a 2D matrix into RGB format for display with IMAGE
% labelCrossAxes                Label crossAxes
% max2                          Overall maximum of a multidimensional matrix
% min2                          Overall minimum of a multidimensional matrix
% rec2sphere                    Transform [y, x, t] coordinates to [az, el, radius]
% sphere2rec                    Transform [az, el, radius] coordinates to [y, x, t]
%
% ----------------- Pyrtools functions ------------------------------------
%
% showIm                Display a 2D matrix as a grayscale image
% pixelAxes             Set the axes of the current plot to avoid aliasing
% range2                Overall min and max of a multidimensional matrix
%
% ----------------- MEX functions ---------------------------------------
%
% blurDn3               Blur and downsample a 3D matrix
% upblur3               Blur and upsample a 3D matrix (NOT YET WORKING)
% validCorrDn3          Do correlation with a 3D matrix
%
