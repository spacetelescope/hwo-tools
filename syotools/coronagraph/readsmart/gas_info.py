def gas_info():
    """
    Loads and returns HITRAN gas information (description)

    Parameters
    ----------
    None

    Returns
    -------
    load_gases.gas_info(): tuple array (mixed type)

    Notes
    -------
    This program loads the gas codes, formula code, descriptor string
    (formula name), mean molecular masses, and wavenumber coverage  of HITRAN gases 
    1-47 and additional noble and volcanic gases from gases.csv in /fixed_input/.'

    There are three gas code designations:

    'index' - the value that corresponds to the position of the gas info in the tuple array.
        This also corresponds to the HITRAN gas index up until gas number '47'.
    'h_index' - HITRAN gas index. Set to '0' if not a HITRAN gas.
    'old_index' - index value of gas for old smart_interface bash code.
        Matches 'index' and 'h_index' through gas index=39. 
    
    The formula code is realized in three forms:

    'Formula' - molecular formula with proper capitalization 
        (string, e.g. 'CH3Cl')
    'FORMULA' - molecular formula all caps (string, e.g. 'CH3CL')
    'formula' - molecular formula all lower case (string, e.g. 'ch3cl')

    These additional attributes are contained in the tuple:
    'name' - common name for the molecule or atomic gas species
    'mass' - mean molecular mass in g/mol
    'minwn' - minimum wavenumber of transitions cataloged in HITRAN 2012 database
    'maxwn' - maximum wavenumber of transitions cataloged in HITRAN 2012 database
    'n'    - index of refraction of the gas at 589 nm at STP (n=-1 if unknown)

    References
    -------
    HITRAN gas codes from: http://www.cfa.harvard.edu/hitran/molecules.html
    Mean molecular masses from Wikipedia: http://www.wikipedia.org
    Refractive indices from: http://www.kayelaby.npl.co.uk/general_physics/2_5/2_5_7.html

    Examples
    -------
    >>> load_gases.gas_info()[:][1]
    (1, 1, 1, 'H2O', 'H2O', 'h2o', 'water', 18.02, 0.0, 25711.0, 1.000256)
    >>> load_gases.gas_info()['FORMULA'][48]
    'HE'
    >>> load_gases.gas_info()['Formula'][48]
    'He'
    >>> load_gases.gas_info()['Formula'][37]
    'HOBr'

    """

    import numpy as np
    import os
    module_path = os.path.dirname(__file__)
    #read gases csv file
    file_name = '/'.join([module_path,'gases.csv'])
    #file_name = 'fixed_input/gases.csv'
    gas_info = np.genfromtxt(file_name,dtype=None,delimiter=',',comments='#',\
           names=['index', 'h_index', 'old_index', 'Formula', 'FORMULA', 'formula', 'name', \
                  'mass', 'minwn', 'maxwn', 'n'])
       
    return gas_info
