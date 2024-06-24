#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 14:41:07 2024

@author: john serences (jserences@ucsd.edu)
"""

import torch
#import numpy as np
from h_layer import *

#--------------------------------
# class to define the recurrent layer object
# this specifies interactions between input
# layer and the hiden layer, as well as interactions
# between units in the hidden layer
#--------------------------------
class RLayer(torch.nn.Module):

    
    #----------------------------
    # set up the instance...
    #----------------------------
    def __init__(self, rnn_settings):
        '''
        Initialize the params for the network
            INPUT
                rnn_settings: dict with network parameters
        '''
        #----------------------------
        # init
        #----------------------------
        super().__init__()
        
        #----------------------------
        # get input params for network
        #----------------------------
        self.h_size = rnn_settings.get('h_size',32)                # size of hidden layer
        self.inp_size = rnn_settings.get('inp_size',1)             # number of inputs, last dim of a [time x batch size x num_inputs] matrix (where batch size == n_trials) 
        self.out_size = rnn_settings.get('out_size',self.inp_size) # output size - defaults to inp_size
        self.dt = rnn_settings.get('dt',1)                         # timestep
        self.tau = rnn_settings.get('tau',20)                      # decay
        self.alpha = self.dt/self.tau                              # determines the influence of prior activation on current activation
        
        act_func = rnn_settings.get('act_func','relu')             # activation function relating x_t to r
        if act_func == 'relu':
            self.act_func = torch.relu
        elif act_func == 'sigmoid':
            self.act_func = torch.sigmoid
        else:
            raise ValueError(f'{act_func} is not a currently supported activation function')

        # pre and post activation function noise...
        self.preact_n = rnn_settings.get('preact_n', 0)            # noise before passing x_t to activation function
        self.postact_n = rnn_settings.get('postact_n', 0)          # noise after passing x_t to the activation function

        #----------------------------
        # determine if weight/bias params are 
        # trainable/not-trainable (requires_grad = True/False)
        #----------------------------

        # W_in weights and bias
        self.W_in_trainable = rnn_settings.get('W_in_trainable',False)
        self.bias_in_trainable = rnn_settings.get('bias_in_trainable',False)

        # W_out weights and bais
        self.W_out_trainable = rnn_settings.get('W_out_trainable',True)
        self.bias_in_trainable = rnn_settings.get('bias_in_trainable',True)

        #----------------------------
        # determine the weights/biases for W_in
        # and W_out
        #----------------------------
        self.W_in = rnn_settings.get('W_in')
        self.bias_in = rnn_settings.get('bias_in')
        self.W_out = rnn_settings.get('W_out')
        self.bias_out = rnn_settings.get('bias_out')
        
        #----------------------------
        # create layers
        #----------------------------
        
        # create the input layer, input size x hidden layer size
        # torch.nn.Linear means linear transform of input stimulus (u * W_in + bias_in), where u == stimulus. 
        self.inp_layer = torch.nn.Linear(self.inp_size, self.h_size)
        
        # create the hidden layer of size [h_size x h_size]
        # also linear (inp_layer * W_hid + bias_hid)
        self.h_layer = HLayer(rnn_settings)

        #----------------------------
        # Then assign weights/requires_grad flag (trainable)
        # to input weights/bias
        #----------------------------
        self.inp_layer.weight = torch.nn.Parameter(self.W_in,requires_grad=self.W_in_trainable)
        self.inp_layer.bias = torch.nn.Parameter(self.bias_in,requires_grad=self.bias_in_trainable)
        
    #--------------------------------
    # define the operations for each recurrent 
    # step...stack these to form the model over time (in forward method)
    #--------------------------------
    def recurrence(self, r_t, x_t, u_t):
        ''' 
        Method to control recurrence 
            INPUT
                r_t: state of hidden (recurrent) units at time t
                x_t: input current at time t
                u_t: stimulus at time t

            OUPUT
                r_t1: state of hidden units at next time point (t+1)
                x_t1: voltage at next time point (t+1)
                
            Note: v_t1 and r_t1 are linked via the activation function (i.e. r_t1 = act_func(x_t1).
            
        '''
        # pass stimulus (u) at time t to input layer, performs linear operation
        # u_t * W_in + b
        w_in_u = self.inp_layer(u_t)  
        
        # then eval hidden layer - influence of 
        # each unit on other units
        # via linear operation r_t * W_hid + b
        w_hid_r = self.h_layer(r_t) 
        
        # update state of hidden layer units for next time step
        # (alpha (dt/tau) determines the weight given to current inputs vs prior state
        x_t1 = (1-self.alpha) * x_t + self.alpha * (w_hid_r + w_in_u)
    
        # additive noise applied before passing through activation function
        if self.preact_n > 0:
            preact_n = torch.randn((u_t.size(0), self.h_size)) * self.preact_n
            x_t1 = x_t1 + self.alpha * preact_n
    
        # apply activation function to get firing rate at next time step
        r_t1 = self.act_func(x_t1)
    
        # additive noise applied after passing through activation function
        if self.postact_n > 0:
            postact_n = torch.randn((u_t.size(0), self.h_size)) * self.postact_n
            r_t1 = r_t1 + postact_n

        # return new states, which will form the basis
        # for the next update at the next time step, and so on...
        return r_t1, x_t1

    #--------------------------------
    # define how stimulus inputs proogate through
    # the network...builds a stack of the network states
    #--------------------------------    
    
    def forward(self, inp):
        """
        Define how inputs propogate through the network by making a stack 
        of states...
            INPUT
                inp: [timepoints x batch_size(num trials) x inp_size]
                
            OUTPUT
                stacked_states: [seq_len x batch size x hidden_size], stack of hidden layer status
        """

        # initial states
        x_t = torch.zeros((inp.size(1), self.h_size), device=inp.device)
        r_t = self.act_func(x_t)
        
        # list of updated states after passing input at each timepoint
        # through recurrent func 
        states = []
        for i in range(inp.size(0)):
            
            r_t, x_t = self.recurrence(r_t, x_t, inp[i])
            # append to the list of states
            states.append(r_t)
  
        return torch.stack(states, dim=0)


#--------------------------------
# make the model object by calling the 
# recurrent object and adding an output
# layer. 
#--------------------------------
class ctRNN(torch.nn.Module):
    
    #--------------------------------
    # init an instance with recurrent layer (input + hidden)
    # and output layer
    #--------------------------------
    def __init__(self, rnn_settings):
        super().__init__()

        # define recrrent layer (recurrent processing of input and hidden layer)
        self.recurrent_layer = RLayer(rnn_settings)
        
        #define the output, or readout, layer
        self.output_layer = torch.nn.Linear(rnn_settings.get('h_size'), rnn_settings.get('out_size',None))

    #--------------------------------
    # Instantiate the forward model    
    #--------------------------------
    def forward(self, inputs):
        
        # get the state of the hidden layer units
        hidden_states = self.recurrent_layer(inputs)
        
        # run it through the output layer
        output = self.output_layer(hidden_states.float())
        
        return output, hidden_states
    
    #--------------------------------
    # Define loss function to use during model training 
    # This will determine the distance between model output and 
    # the desired (target) output during supervised learning
    # This is just the mean squared error (MSE), which you could
    # also implement with nn.MSELoss() but this example can be 
    # modified if you want a custom loss func     
    #--------------------------------
    def mse_loss(self, outputs, targets):
        '''
        INPUT
            output: [time x trial_in_batch x output_size] model output of type torch.Tensor (float)
            target: [time x trial_in_batch] target of type torch.Tensor (float)
    
        OUTPUT
            loss, mean squared error (MSE) between outputs and targets. type torch.Tensor
        
        writing it out step by step for clarity, but more compact and harder to read is:
            torch.divide(torch.sum(torch.square(torch.subtract(torch.squeeze(output),target).flatten())),output.shape[0] * output.shape[1])
        '''
    
        # compute difference between output and target (squeeze in case output_size == 1)
        # flatten to vectorize the time x batch size (trials) matrices before passing to 
        # subsequent ops
        o_t_diff = torch.subtract(torch.squeeze(outputs),targets).flatten()
    
        # square the diff cause don't care about sign (and squared penalizes far errors)
        o_t_sq_diff = torch.square(o_t_diff)
    
        # sum of squares
        o_t_ss_diff = torch.sum(o_t_sq_diff)
        
        # divide by number of data points to get mean squared error
        o_t_mse = torch.divide(o_t_ss_diff, o_t_diff.shape[0])
        
        return o_t_mse
        