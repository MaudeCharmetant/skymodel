# -*- coding: utf-8 -*-
# 
#  This file is part of the ccatp_sky_model.
# 
#  sky_model is free software; you can redistribute it and/or modify
#  it under the terms of the MIT License.
# 
#  sky_model is distributed in the hope that it will be useful,but 
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See 
#  the provided copy of the MIT License for more details.

"""ccatp_sky_model provides a model of the microwave sky that will be used for 
   science forecasts for the CCAT-prime survey. The sky model combines maps from
   numerical simulations of extragalactic sources with the most recent Galactic
   foregrounds maps based on Planck data.
"""

__version__ = "1.0"

__bibtex__ = """

"""

from .ccatp_sky_model import (convert_units, px_size, sample_sphere_uniform, project_maps, simulate_gal_foregrounds, simulate_cib, simulate_radio_ps, simulate_cmb, simulate_tSZ, simulate_kSZ, simulate_white_noise, simulate_atmosphere, simulate_iras_ps, simulate_nvss_ps, ccatp_sky_model)
