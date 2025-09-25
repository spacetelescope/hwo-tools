import sys
import pickle

import pytest
import numpy as np
from bokeh.models import ColumnDataSource

from main import update_data, initialize_setup

from syotools.spectra.spec_defaults import pysyn_spectra_library


hwo = None
pivotwave = None
template_to_start_with = "Flat (AB)"

# to do: Allow nested values in the arguments (vary over values that are only valid/available
# given some other value)
def VaryMany(*args, **keys):
    inlist = []
    maxlen = 0
    # what is the longest list of parameters to vary over?
    for arg in args:
        #print(arg)
        if len(arg) > maxlen:
            maxlen = len(arg)
    # shuffle (and expand) all parameters to make random combinations
    expandargs = []
    for idx,arg in enumerate(args):
        #print(len(arg))
        # first, we can make repeats of the array as needed by resizing it to
        #  something larger
        arg = np.resize(arg,maxlen)
        #print(len(arg))
        # now we shuffle. We don't want to use randint() because every option must be
        #  picked at least once. This is an in-place operation.
        np.random.seed(idx)
        np.random.shuffle(arg)
        expandargs.append(arg)
    # now we have full lists of arguments and items, resized to the same length. Pair them up.
    for x in range(maxlen):
        calc = []
        for arg in expandargs:
            calc.append(arg[x])
        inlist.append(calc)

    return inlist

def check_relative_diff(actual, expected, rel_tol=0.1):
    """
    Simple function to check if two lists are approximately equal within a relative tolerance.
    Reports percentage differences for values that exceed the tolerance.

    Args:
        actual: List of actual values
        expected: List of expected values
        rel_tol: Relative tolerance (default: 0.1 or 10%)

    Returns:
        True if all values are within tolerance, False otherwise
    """
    if len(actual) != len(expected):
        print(f"Lists have different lengths: actual={len(actual)}, expected={len(expected)}")
        return False

    all_within_tolerance = True
    differences = []

    for i, (a, e) in enumerate(zip(actual, expected)):
        for j, (act, exp) in enumerate(zip(actual[i], expected[i])):
            if exp == 0:
                # Can't calculate relative difference if expected is zero
                if act != 0:
                    all_within_tolerance = False
                    differences.append((i, a, e, "inf"))
            else:
                rel_diff = abs(act - exp) / abs(exp)
                pct_diff = 100 * rel_diff

                if rel_diff > rel_tol:
                    all_within_tolerance = False
                    differences.append((i, j, act, exp, pct_diff))

    if not all_within_tolerance:
        print("\nValues exceeding relative tolerance:")
        for i, j, act, exp, pct in differences:
            if pct == "inf":
                print(f"  Index {i},{j}: actual={act}, expected={exp}, difference=infinite (division by zero)")
            else:
                print(f"  Index {i},{j}: actual={act:.6g}, expected={exp:.6g}, difference={pct:.2f}%")

    return all_within_tolerance


initialize_setup()

class item():
    def __init__(self, value=30):
        self.value = value

exptime = item(30.0)
template = item("Flat (AB)")
magnitude = item(20.0)
aperture = item(6.0)
redshift = item(0.0)
grating = item("G150M (R =30,000)")

exptimes = np.linspace(1, 20, 20)
templates = list(pysyn_spectra_library)
magnitudes = np.linspace(20, 35, 1)
apertures = np.linspace(4,10, 1)
redshifts = np.linspace(0, 3, 20)
gratings = ["G120M (R = 30,000)", "G150M (R = 30,000)", "G180M (R = 30,000)", "G155L (R = 5,000)", "G145LL (R = 500)"]

test_setups = VaryMany(exptimes, templates, apertures, magnitudes, redshifts, gratings)

def do_comparisons(inputs):
    global exptime
    global template
    global magnitude
    global aperture
    global redshift
    global gratings
    global pivotwave

    exptime.value = inputs[0]
    template.value = inputs[1]
    magnitude.value = inputs[2]
    aperture.value = inputs[3]
    redshift.value = inputs[4]
    grating.value = inputs[5]
    snr_results, spectrum_template = update_data(None, None, None)

    print(snr_results.data)

    return snr_results.data["sn"], spectrum_template.data["f"]

def create_comparisons(reset):
    test_results = []
    for setup in test_setups:
        test_results.append((setup, do_comparisons(test_setups)))
    if reset:
        with open("test_uvspec.pickle", "wb") as picklefile:
            pickle.dump(test_results, picklefile)

'''
LOAD IT
'''
try:    
    with open("test_uvspec.pickle", "rb") as picklefile:
        test_results = pickle.load(picklefile)
except FileNotFoundError:
    create_comparisons(True)
    with open("test_uvspec.pickle", "rb") as picklefile:
        test_results = pickle.load(picklefile)

@pytest.mark.parametrize("testset", test_results)
def test_camera(testset):
    
    inputs = testset[0]
    expected = testset[1]
    print(inputs)
    snr_results, spectrum_template = do_comparisons(inputs)

    assert check_relative_diff(expected, (snr_results, spectrum_template))


if __name__ == "__main__":
    reset=False
    if (len(sys.argv) > 1) and (sys.argv[1] == "reset"):
        reset = True
    create_comparisons(reset)
