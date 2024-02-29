# ArrayDesignOpt
Python files developed to simulate 2D random array patterns, and transform regular uniform arrays to pseudo random through side lobe minimization.<br />
Similar to the LWA (long wavelength array) stations, the stations for the NNNNNN Peru Project were designed using the Leonid Kogan's algorithm. <br />
<br />
<br />
This library is still in development and testing :) ...



## Instalation
Requierement libraries:
- numpy v1.19.5
- matplotlib v3.6.0
- shapely
- IPython
- scipy
------------------------------------------------------- <br />
git clone https://github.com/joabAM/arrayDesignOpt.git <br />
or git clone -b main http://intranet.igp.gob.pe:8082/arrayDesignOpt <br />
cd arrayDesignOpt  <br />
python setup.py bdist_wheel <br />
pip install /path/to../dist/arrayDesignOpt-0.1-py3-none-any.whl --force-reinstall (optional) <br />

## Usage
------------------------------------------------------- <br />
Simple script <br />
------------------------------------------------------- <br />
from arrayDesignOpt import main  <br />
arr = main.AntennaArray(N=64) <br />
norm_pattern, _domg, _cosxy, _radmin, _pat_int  = arr.getPattern() <br />
SLLdB, XBest, YBest, BeamOpt, Xlast1, Ylast1 = arr.minimizeSLL( _radmin,gain=10, plotPattern=True) <br />
