Data set from Pillow et al (2008)
http://www.ncbi.nlm.nih.gov/pubmed/18650810


LICENSE: You may use these data for scientific exploration and teaching. You may not use these data in a publication or for profit.


1) SpTimesRGC.mat - 1 x 27 cell array: each cell is a vector of spike times from a single neuron in response to a 20 minute stimulus.  Importantly, the spike times are in units of STIMULUS FRAMES, where the stimulus frame rate was 119.9820Hz (meaning that one frame is 8.33458 ms, so a spike time of "10" corresponds to a spike time at 0.0833 seconds).  

Cells #1-16 are OFF cells, and #17-27 are ON cells.

2)  Stim_reduced.mat - This 144000 x 100 matrix has the stimulus movie, where each row corresponds to a 10 x 10 spatial stimulus.  (This is the stimulus that elicited the spikes in SpTimesRGC.mat).

Then following two files have spike times and a stimulus for a second experiment consisting of 600 repeats of a 10-second stimulus (useful for computing a PSTH, testing predictions, etc):

3)  MtspRGCrpt.mat - 1 x 27 cell array: each contains a matrix spike times emitted in response to 600 repeats of a 10-second stimulus. (Each column in a matrix is the spike times from a single trial, padded with zeros at the end).

4) Stim_reducedRpt.mat - the repeat stimulus

There's a matlab script in the directory ("loadDataAndComputeSTA.m") that will load the data and compute the spike-triggered average of the first neuron, just in case that's helpful.
