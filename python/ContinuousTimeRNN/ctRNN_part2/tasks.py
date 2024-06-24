#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 8 09:20:16 2024

@author: john serences (jserences@ucsd.edu)

Class to define tasks for training/evaluting continuous time RNN
go/no-go, delayed match to sample, and mante tasks adapted from Kim,Li,Sejnowski 2019 PNAS 
should be easy to make your own following the same general format
"""

#---------------------------------------------------------------
# imports
#---------------------------------------------------------------
import numpy as np
import torch

#---------------------------------------------------------------
# class continuous time tasks (ctTASKS)
#---------------------------------------------------------------
class ctTASKS:
    """
    Firing-rate RNN model for excitatory and inhibitory neurons
    Initialization of the firing-rate model with recurrent connections
    """
    def __init__(self, settings):

        # make a list of expected key-value pairs (params)...
        # raise error if missing... probably a nicer way of doing this, but...
        # and this won't catch bad values...just missing values
        self.task = settings['task']

        if self.task == 'go-nogo':
            # expected param list for this task
            exp_p = ['T', 'stim_on', 'stim_dur', 'n_trials','acc_amp_thresh']
            
            # check for and assign task-relevant params
            self.param_check(settings,exp_p)
            
            # if passes check, assign
            self.T = settings.get('T')
            self.stim_on = settings.get('stim_on')
            self.stim_dur = settings.get('stim_dur')
            self.n_trials = settings.get('n_trials')  
            self.acc_amp_thresh = settings.get('acc_amp_thresh')
            
        elif self.task == 'dmts':
            # expected param list for this task
            exp_p = ['T', 'stim_on', 'stim_dur', 'delay', 'n_trials','acc_amp_thresh']
            
            # check for and assign task-relevant params
            self.param_check(settings,exp_p)
            
            # if passes check, assign
            self.T = settings.get('T')
            self.stim_on = settings.get('stim_on')
            self.stim_dur = settings.get('stim_dur')
            self.delay = settings.get('delay')  
            self.n_trials = settings.get('n_trials')  
            self.acc_amp_thresh = settings.get('acc_amp_thresh')
        
        elif self.task == 'mante':
            # expected param list for this task
            exp_p = ['T', 'stim_on', 'stim_dur','n_trials','acc_amp_thresh']
            
            # check for and assign task-relevant params
            self.param_check(settings,exp_p)
            
            # if passes check, assign
            self.T = settings.get('T')
            self.stim_on = settings.get('stim_on')
            self.stim_dur = settings.get('stim_dur')
            self.n_trials = settings.get('n_trials')  
            self.acc_amp_thresh = settings.get('acc_amp_thresh')

        else:
            raise ValueError(f'{self.task} is not a supported task')            
        
        # success!
        print(f'{self.task} task has been initialized')
        
    #---------------------------------------------------------------
    # check that all neccessary task params are in place
    #---------------------------------------------------------------        
    def param_check(self, settings, exp_p):
        # do the parm check for this task
        for p in exp_p:
            if p not in settings:
                raise ValueError(f'Missing stim settings param {p}')

    #---------------------------------------------------------------
    # Go/No-Go task (i.e. 'respond' if stim present, otherwise no response)
    #---------------------------------------------------------------
    def stim_go_nogo(self):
        '''
        
        Generate the [time x trial] input stimulus matrix 
        for the go-nogo task
    
        INPUT
            settings: dictionary with task specs
                T: duration of a single trial (in steps)
                stim_on: stimulus starting time (in steps)
                stim_dur: stimulus duration (in steps)
                n_trials: number of trials to generate for each batch
                
        OUTPUT
            u: [T x n_trials] stimulus matrix, type torch.Tensor (float). 
            tri_type: vector indicating trial type (1 for a "go" trial or 0 for "nogo" trial)
                tri_type is input to targ_go_nogo method that generates model output target 
                for each trial
            
        '''
        
        # alloc storage for stim time series 'u' and trial type
        # will use the trial labels to generate 'target' outputs
        # for model training in fucntion below...
        u = torch.zeros((self.T,self.n_trials,1)) 
        tri_type = np.full(self.n_trials, torch.nan)
    
        # generate all trials (the full batch of inputs)
        for nt in range(self.n_trials): 
            # determine if go/no-go trial
            # go == 1, nogo has all zeros
            if torch.rand(1) <= 0.5:
                u[self.stim_on:self.stim_on+self.stim_dur,nt,0] = 1
                tri_type[nt] = 1
                
            else:
                tri_type[nt] = 0 
            
        # convert stim to torch tensor type
        # and return stims + tri_type
        #return torch.from_numpy(u).type(torch.float), tri_type
        return u, tri_type
    
    #---------------------------------------------------------------
    # Delayed match to sample task: stim1, delay, stim2: same or different?
    #---------------------------------------------------------------
    def stim_dmts(self):
        """
        Generate the [time x trial] input stimulus matrix 
        for the dmts task
    
        INPUT
            settings: dict containing the following keys
                T: duration of a single trial (in steps)
                stim_on: stimulus starting time (in steps)
                stim_dur: stimulus duration (in steps)
                delay: delay b/w two stimuli (in steps)
                n_trials: number of trials to generate for each batch

        OUTPUT
            u: T x n_trials x 2 stimulus matrix (stim 1 - memory sample, stim 2 - memory probe)
            tri_type: 1 for 'same' or 0 for 'diff'
        """
        
        # alloc storage for stim time series 'u' and trial type
        # will use the trial labels to generate 'target' outputs
        # for model training in fucntion below...
        u = torch.zeros((self.T,self.n_trials,2)) 
        tri_type = np.full(self.n_trials, torch.nan)
        
        # DMTS task
        # generate all trials (the full batch of inputs)
        for nt in range(self.n_trials): 
            
            # first stimulus
            labs = []
            if torch.rand(1) < 0.50:
                u[self.stim_on:self.stim_on+self.stim_dur,nt,0] = 1
                labs.append(1)
            else:
                u[self.stim_on:self.stim_on+self.stim_dur,nt,0] = -1
                labs.append(-1)
        
            # second stimulus after delay period
            if torch.rand(1) < 0.50:
                u[self.stim_on+self.stim_dur+self.delay:self.stim_on+2*self.stim_dur+self.delay,nt,1] = 1
                labs.append(1)
            else:
                u[self.stim_on+self.stim_dur+self.delay:self.stim_on+2*self.stim_dur+self.delay,nt,1] = -1
                labs.append(-1)
        
            if np.prod(labs) == 1:
                tri_type[nt] = 1
            else:
                tri_type[nt] = 0
    
        #return torch.from_numpy(u).type(torch.float), tri_type
        return u, tri_type

    def stim_mante(self):
        """
        Generate the input stimulus matrix for the
        mante task (i.e. response depends on context)
    
        INPUT
            settings: dict containing the following keys
                T: duration of a single trial (in steps)
                stim_on: stimulus starting time (in steps)
                stim_dur: stimulus duration (in steps)

        OUTPUT
            u: T x n_trials x 4 stimulus matrix (first 2 inputs for motion/color and the second
            2 rows for context), type torch.Tensor (float). 
            tri_type: either +1 or -1
        """
        # alloc storage for stim time series 'u' and trial type
        # will use the trial labels to generate 'target' outputs
        # for model training in fucntion below...
        u = torch.zeros((self.T,self.n_trials,4)) 
        tri_type = np.full(self.n_trials, torch.nan)
        
        # Mante task
        # generate all trials (the full batch of inputs)
        for nt in range(self.n_trials): 
            # Color/motion sensory inputs

            u_lab = torch.zeros((2, 1))
            if torch.rand(1) <= 0.50:
                u[self.stim_on:self.stim_on+self.stim_dur,nt,0] = torch.randn(self.stim_dur) + 0.5
                u_lab[0, 0] = 1
            else:
                u[self.stim_on:self.stim_on+self.stim_dur,nt,0] = torch.randn(self.stim_dur) - 0.5
                u_lab[0, 0] = -1
        
            if torch.rand(1) <= 0.50:
                u[self.stim_on:self.stim_on+self.stim_dur,nt,1] = torch.randn(self.stim_dur) + 0.5
                u_lab[1, 0] = 1
            else:
                u[self.stim_on:self.stim_on+self.stim_dur,nt,1] = torch.randn(self.stim_dur) - 0.5
                u_lab[1, 0] = -1
        
            # Context input
            if torch.rand(1) <= 0.50:
                u[:,nt,2] = 1
        
                if u_lab[0, 0] == 1:
                    tri_type[nt] = 1
                elif u_lab[0, 0] == -1:
                    tri_type[nt] = -1
            else:
                u[:,nt,3] = 1
        
                if u_lab[1, 0] == 1:
                    tri_type[nt] = 1
                elif u_lab[1, 0] == -1:
                    tri_type[nt] = -1
    
        #return torch.from_numpy(u).type(torch.float), tri_type
        return u, tri_type

    
    #---------------------------------------------------------------
    # target signal for go-nogo
    #---------------------------------------------------------------
    def target_go_nogo(self, tri_type):
        '''
        Generate the [time x trial] input stimulus matrix 
        for the go-nogo task
    
        INPUT
            settings: dictionary with task specs
                T: duration of a single trial (in steps)
                stim_on: stimulus starting time (in steps)
                stim_dur: stimulus duration (in steps)
                n_trials: number of trials to generate for each batch
            tri_type: marks each trial generated by stim_go_nogo as "go" or "nogo"
                
        OUTPUT
            targs: [T x n_trials] stimulus matrix, type torch.Tensor (float)
            
        '''
        
        # define target outputs for model training
        targs = torch.zeros((self.T, self.n_trials))
        for nt in range(self.n_trials):
            if tri_type[nt] == 1:
                targs[self.stim_on+self.stim_dur:,nt] = 1
    
        # convert to torch tensor type
        # and return targs
        return targs
        #return torch.from_numpy(targs).type(torch.float)
    
    #---------------------------------------------------------------
    # target signal for dmts
    #---------------------------------------------------------------    
    def target_dmts(self, tri_type):
        """
        Method to generate a continuous target signal (z) 
        for the dmts task
    
        INPUT
            settings: dict containing the following keys
                T: duration of a single trial (in steps)
                stim_on: stimulus starting time (in steps)
                stim_dur: stimulus duration (in steps)
                delay: delay b/w two stimuli (in steps)
                n_trial: number of trials in batch
                tri_type: n_trial array of 1 ('same') or 0 ('diff')
        
        OUTPUT
            targs: [T x n_trial] target signal
        """

        # task end time step
        task_end_T = self.stim_on+2*self.stim_dur + self.delay

        # define target outputs for model training
        targs = torch.zeros((self.T, self.n_trials))
        for nt in range(self.n_trials):
            if tri_type[nt] == 1:
                targs[10+task_end_T:10+task_end_T+100,nt] = 1
                
            elif tri_type[nt] == 0:
                targs[10+task_end_T:10+task_end_T+100,nt] = -1
    
        #return torch.from_numpy(targs).type(torch.float)
        return targs

    #---------------------------------------------------------------
    # target signal for mate
    #---------------------------------------------------------------   
    def target_mante(self, tri_type):
        """
        Method to generate a continuous target signal (z) 
        for the MANTE task
    
        INPUT
            settings: dict containing the following keys
                T: duration of a single trial (in steps)
                stim_on: stimulus starting time (in steps)
                stim_dur: stimulus duration (in steps)
                taus: time-constants (in steps)
                DeltaT: sampling rate
            label: either +1 or -1
        OUTPUT
            z: 1xT target signal
        """
    
        targs = torch.zeros((self.T, self.n_trials))
        for nt in range(self.n_trials):
            if tri_type[nt] == 1:
                targs[self.stim_on+self.stim_dur:,nt] = 1
            else:
                targs[self.stim_on+self.stim_dur:,nt] = -1
    
        #return torch.from_numpy(targs).type(torch.float)
        return targs

    #---------------------------------------------------------------
    # Compute task accuracy using defined criteria
    #---------------------------------------------------------------
    def compute_acc(self, settings,outputs,targets): 

        # if task == go-nogo
        if  self.task == 'go-nogo':
        
            # compute mean of model output on each trial in last batch over 
            # last timepoints in each trial
            mean_out = np.mean(np.squeeze(outputs.detach().numpy())[self.stim_on+self.stim_dur:,:],axis=0)
            
            # if output was above threshold, then model guessed 'go' trial
            mean_out[mean_out>self.acc_amp_thresh[0]] = 1
            
            # if output was below threshold, then model guessed 'nogo' trial
            mean_out[mean_out<self.acc_amp_thresh[1]] = 0
            
            # compute max of target to see if each trial was actually go (1) or nogo (0) 
            mean_targ = np.max(targets.detach().numpy(),axis=0)    
        
            # compute mean acc over all trials in the last batch
            acc = np.mean(mean_out == mean_targ)
            
        # if task == dmts
        elif  self.task == 'dmts':
        
            task_end_T = self.stim_on+2*self.stim_dur + self.delay
            
            # compute mean of model output on each trial in last batch over 
            # last timepoints in each trial
            # abs because the target value is either 1 or -1 and abs allows
            # us to just compare to one threshold for correct/incorrect...might yeild 
            # an occasional false positive early in training, but not likely as training 
            # advances.
            mean_out = np.abs( np.mean(np.squeeze(outputs.detach().numpy())[10+task_end_T:10+task_end_T+100,:],axis=0) )
            
            # if output was above threshold, then model correct
            mean_out[mean_out>self.acc_amp_thresh] = 1
            
            # if output was below threshold, then model guessed 'nogo' trial
            mean_out[mean_out<self.acc_amp_thresh] = 0
        
            # compute mean acc over all trials in the last batch
            acc = np.mean(mean_out)      
            
        else:
            task = settings['task']
            raise ValueError(f'{task} task not supported by compute_acc')
        
        return acc
