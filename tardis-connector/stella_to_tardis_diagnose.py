import logger
import numpy as np
import astropy.units as u

from tardis.io.model import read_stella_model
from tardis.util.base import is_valid_nuclide_or_elem


logger = logging.get_logger(__name__)
logging.basicConfig(level=logging.INFO)



def get_stella_model(mdl, x_col, y_col):

    """
    Parameters:
    x_col: first col
    y_col: second_col
    -----------
    mdl: a stella model


    Returns:
    -------------------

    The required data.
    """
    # Ideally we can read any set of columns
    stella_data = read_stella_model(mdl).data

    return stella_data[x_col], stella_data[y_col]



def check_homology(tau_limit):

    """
    tau_limit: limit reading the model upto a given tau

    Return
    --------

    check if the models have inner velocities increasing monotonically.

    """
    
    

     




