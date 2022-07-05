
# include the import before the functions so that when we call this module later we load numpy once. 
import numpy as np

def define_filt(fx, filt_freq, type_filt):

    """Define a simple lp or a hp brick wall filter - don't use this in practice, just for tutorial

    Args:
        fx (array): list of frequencies, taken from np.fft.rfftfreq
        filt_freq (float): cutoff frequency
        type_filt (str): low pass or high pass filter, options 'lp', 'hp'

    Returns:
        brick: brick wall filter in frequency domain

    """
    # find the freq in our FFT range that is closest to our desired cutoff point
    cutoff_pnt = np.argmin(np.abs(fx-filt_freq)) # more general way to do it...

    # check filter type...if lp then do nothing, else if hp invert, else return msg

    # init filter with all 0's
    brick = np.zeros(len(fx))

    if type_filt == 'lp':
        # make the filter
        brick[0:cutoff_pnt] = 1
    
    elif type_filt == 'hp':
        # make the filter
        brick[cutoff_pnt:] = 1
    
    else:    
        print('error - specify lp or hp filter')
        return 0
    
    return brick

# Another function to apply the filter. this second function will be in the same module      
def apply_filt(input_sig, input_filter):
    """Apply a filter to an input timeseries (using freq domain multiplication)

    Args:
        input_sig (float): timeseries to be filtered (in time domain)
        input_filter (float): filter to apply to input_sig (in frequency domain)

    Returns:
        filt_sig (float array): filtered signal in time domain

    """
    
    # fft our time domain signal
    fft_sig = np.fft.rfft(input_sig)

    # multiply in freq domain, then ifft to go back into the time domain
    return np.fft.irfft(fft_sig*input_filter)


# Another function to apply the filter. this second function will be in the same module      
def apply_filt_2D(input_sig, input_filter):
    """Apply a filter to an input timeseries (using freq domain multiplication)

    Args:
        input_sig (float): timeseries to be filtered (in time domain). To keep it simple for our sample
        data set, time will run over rows in the data so we will fft with axis=1
        
        input_filter (float): filter to apply to input_sig (in frequency domain)

    Returns:
        filt_sig (float array): filtered signal in time domain

    """
    r,c = input_sig.shape
    
    # fft our time domain signal
    fft_sig = np.fft.rfft(input_sig, axis = 1)
    
    # make our filter the same size as the data. 
    filt_2d = np.tile(input_filter,(r,1))
    
    # multiply in freq domain, then ifft to go back into the time domain
    return np.fft.irfft(fft_sig*filt_2d)
