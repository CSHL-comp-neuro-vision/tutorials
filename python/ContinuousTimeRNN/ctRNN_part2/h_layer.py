#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 12:22:59 2024

@author: johnserences
"""

import torch
import numpy as np

class HLayer(torch.nn.Module):
    '''
    create a custom hidden layer that supports sparse connections and 
    that can enforce strict assignment of exc and inh units (Dale's principle)
    
    '''
    
    def __init__(self, rnn_settings):
        """
        In the constructor we instantiate four parameters and assign them as
        member parameters.
        """
        super().__init__()
       
        #----------------------------
        # basic params defining the layer
        #----------------------------
        self.h_size = rnn_settings.get('h_size',32)                # size of hidden layer 
        
        # set probability of connections in hidden layer (p_rec) and prob that
        # a connection is inhibitory (p_inh)...apply Dale's principle? 
        self.p_rec = rnn_settings.get('p_rec', 0.2)
        self.p_inh = rnn_settings.get('p_inh', 0.2)
        self.apply_dale = rnn_settings.get('apply_dale', True)
        self.w_dist = rnn_settings.get('w_dist', 'gaus')
        self.w_gain = rnn_settings.get('w_gain', 1.5)
        # W_hid (hidden) weights and bias
        self.W_hid_trainable = rnn_settings.get('W_hid_trainable',True)
        self.bias_hid_trainable = rnn_settings.get('bias_hid_trainable',False)
        self.bias_hid = rnn_settings.get('bias_hid', torch.zeros(self.h_size))
        
        #----------------------------
        # create the hidden layer of size [h_size x h_size]
        #----------------------------
        # self.h_layer = torch.nn.Linear(self.h_size, self.h_size)
        
        # define cell types
        self.exc, self.inh, self.exc_size, self.inh_size = self.define_cell_type()
                
        # Then assign weights/requires_grad flag (trainable)
        # to hidden weights/bias using p_rec,p_inh,Dale's principle, etc
        self.weight, self.mask = self.init_W_hid() 
        
        # bias of layer
        self.bias = torch.nn.Parameter(self.bias_hid,requires_grad = rnn_settings.get('bias_hid_trainable',False))
        
    #--------------------------------
    # Define cell types (exc/inh) and set up to either
    # follow Dale's principle or not...
    #--------------------------------
    def define_cell_type(self):

        """
        Randomly assign units as exc or inh based on desired
            proportion of inhibitory cells
            Do so in accordance with Dale's principle (or not)

        Returns
            exc: bool marking excitatory units
            inh: bool marking inhibitory units
            exc_size: number of excitatory units
            inh_size: number of inhibitory units

        """
        
        # If applying Dale's principle
        if self.apply_dale == True:
            
            # index of inh units based on desired proportion (p_inh)
            inh = torch.rand(self.h_size) < self.p_inh
            
            # if not inhibitory, then excitatory
            exc = ~inh
            
            # number of inh units (inh_size) and 
            # number of exc units (exc_size)
            inh_size = len(torch.where(inh == True)[0])
            exc_size = self.h_size - inh_size

        # If not applying Dale's principle
        elif self.apply_dale == False:
            
            # no separate inhibitory units defined
            inh = torch.full((self.h_size,),False) 
            exc = torch.full((self.h_size,),True)
            inh_size = 0
            exc_size = self.h_size

        return exc, inh, exc_size, inh_size


    #--------------------------------
    # Initialize custom weight matrix for hidden layer... 
    # apply Dale's principle if desired
    #--------------------------------
    def init_W_hid(self):

        '''
        Generate a connectivity weight matrix for the hidden layer W_hid
        using either a gaussian or gamma distribution.
        
        INPUTS:
            P_rec: probability of recurrent connections between units in the 
                hidden layer
            inh: [h_size x h_size] matrix indicating which connections should be 
                inhibitory
            w_dist: distribution to determine weights (gaussian or gamma)
            gain: scale factor for the weights
            apply_dale: apply Dale's principle? 
            note: can add more control over the gaussian/gamma distributions
                but here using values from Kim et al PNAS 2019
        
        OUTPUTS:
            w: [h_size x h_size] matrix of weights 
            m: mask of size [h_size x h_size] of 1's (excitatory units)
                  and -1's (for inhibitory units)
            bias: hidden layer bias
        Final weight matrix is w*mask as implemented in the recurrence method
        '''
        
        # Weight matrix [h_size x h_size] matrix
        w_hid = torch.zeros((self.h_size, self.h_size), dtype = torch.float32)
        ind = torch.where(torch.rand(self.h_size, self.h_size) < self.p_rec)
        
        if self.w_dist == 'gamma':
            w_hid[ind[0], ind[1]] = np.random.gamma(2, 0.003, len(ind[0]))
            
        elif self.w_dist == 'gaus':
            w_hid[ind[0], ind[1]] = torch.normal(torch.zeros(len(ind[0])), torch.ones(len(ind[0])))
            w_hid = w_hid/torch.sqrt(torch.tensor(self.h_size) * torch.tensor(self.p_rec)) * self.w_gain 
            
        # if using Dale's law, make all weights initially positive
        # so that we can ensure that the mask of 1's and -1's has a consistent effect
        # (otherwise it might turn some positive values negative and vice versa depending
        # on the initial sign of the weights)
        if self.apply_dale == True:
            
            # abs weights
            w_hid = torch.abs(w_hid)
        
            # mask matrix - set desired proportion of units to be inhibitory
            mask = torch.eye(self.h_size, dtype=torch.float32)
            mask[torch.where(self.inh==True)[0], torch.where(self.inh==True)[0]] = -1
            
            # convert to torch param, assume not trainable (requires_grad = False)
            mask = torch.nn.Parameter(mask, requires_grad = False)
        
        # if not enforcing Dale's then just make the mask all ones (no effect)
        # not going to use it anyway...
        else:
            mask = torch.ones((self.h_size,self.h_size), dtype=torch.float32)

        return torch.nn.Parameter(w_hid, requires_grad = self.W_hid_trainable), mask       
    
    
    #--------------------------------
    # Define operations on forward sweep through
    # hidden layer
    #--------------------------------
    def forward(self, r_t):
        """
        take a tensor of input data and return
        a tensor of output data 
        
        INPUT: r_t current 'firing rates' of hidden layer units
        
        OUTPUT: r_t * W_hid + bias, either under Dale's principle or 
                    not...
                    
        Note: Never re-assign/modify the weights, just use them to compute output.
                Directly modifying the weights messes up auto-grad
        """        
        
        # make sure that exc/inh currents are maintained
        # if enforcing Dale's principle
        if self.apply_dale == True:
            # rectify weights: max(0,weight) or abs (see Song et al., 2016 PLoS Comp Bio)             
            # then multiply by the mask  
            w = torch.matmul(torch.relu(self.weight), self.mask)
        
        else:
            # leave weight unchanged
            w = self.weight

        # compute output for each trial in the current batch
        out = torch.matmul(r_t,w.T) + self.bias
        
        # return...
        return out
    