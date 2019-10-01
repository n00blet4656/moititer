"""Fits/plots data to a sigmoid curve, and determines ratio of virus volume(in uL) to cell count to maintain target MOI.

    This module employs four functions to take fluorescence data generated from a resazurin read, to titer the virus via fitting to a sigmoid curve, according to an MOI described by the user.

Standard workflow of this module is as follows:
input_data -> data_proc -> curve_opt_fit

This module operates on the presumption that negative and positive controls are being measured, as this gives a basis to scale/normalize the dataset.

To use the function input_data, your sheet must consist of two columns: one for the measurement and one for the amount of virus in the well.
"""




def input_data(file, sheet, n, name):
    """Processes data from Excel File
    
    Parameters
    ----------
    file: filepath
    
    sheet: name of the sheet
    
    n: number of data points
    
    name: two-object list that describes your x and y values(let first column be your virus amount, second be the response recorded

    
    Returns
    -------
    out : pandas dataframe
    """
    import pandas as pd
    return pd.read_excel(file, sheet, nrows = n, header = None, names=name)

def sigfunc(x, c):
    """Function for sigmoid curve fitting
    
    Parameters
    ----------
    
    x: int or float
    
    c: EC50 value"""
    import numpy as np
    return (1)/(1+ 10**(np.log10(c)-x))

def data_proc(data, label1, label2, pos, neg):
    """sorts the dataframe from controls and data, and normalizes data based on controls
    
    Parameters
    ----------
    data: resulting dataframe from input data function
    
    label1: name of the first column described in input_data
    
    label2: name of the second column described in input_data
    
    pos: integer/string that identifies the positive control in the column col
    
    neg: integer/string that identifies the negative control in the column col
    
    Returns
    -------
    
    out : pandas dataframe
    """  
    import pandas as pd
    data_rfu = data.loc[~((data[label1] == pos) | (data[label1] == neg))]
    data_rfu = data_rfu.groupby(label1, as_index = False).mean()
    
    data_ctrl = data.loc[((data[label1] == pos) | (data[label1] == neg))]
    data_ctrl = data_ctrl.groupby(label1, as_index = False).mean()
    
    positive = data_ctrl.loc[(data_ctrl[label1] == pos)]
    positive = positive.mean()
    
    negative = data_ctrl.loc[data_ctrl[label1] == neg]
    negative = negative.mean()
    
    data_rfu[label2] = data_rfu[label2].apply(lambda x: (x - negative)/(positive - negative))
    return data_rfu
                                            
def curve_opt_fit(data,label1, label2, moi, cell_ct):
    """ fits data from data_proc to the sigmoid function sigfunc

    Parameters
    ----------
    data: dataframe produced by data_proc function
    
    label1: name of the first column described in input_data
    
    label2: name of the second column described in input_data
    
    moi = target moi
    
    cell_ct = cell count per well
    
    Returns
    -------
    
    out : 
        
        graph of optimized curve fit with data points
        
        return string describing amount of virus to transduce per cell amount to maintain desired MOI"""
    import numpy as np
    from scipy import optimize
    import matplotlib.pyplot as plt
    x = data[label1]
    y = data[label2]
    popt, pcov = optimize.curve_fit(sigfunc, xdata = x, ydata = y)
    x1 = np.linspace(-20, 20, num = 40000)
    plt.plot(x1, sigfunc(x1, *popt))
    plt.plot(x, y, 'ro')
    plt.xlabel('virus')
    plt.ylabel('response')                                     
    virus_transd = (np.log10(popt) - np.log10(1/moi - 1))                                       

    return ("You need to transduce " + str(virus_transd) +" uL of your virus to every " + str(cell_ct) + " cells to maintain an MOI of " + str(moi))                                       
                                              