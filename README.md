# ArrayDesignOpt
Python files developed to simulate 2D random array patterns, and transform regular uniform arrays to pseudo random through side lobe minimization.



## Instalation
Requierement libraries:
- numpy v1.19.5
- matplotlib v3.6.0
- shapely
- IPython
- scipy

git clone https://github.com/joabAM/arrayDesignOpt.git <br />
cd arrayDesignOpt  <br />
python setup.py bdist_wheel <br />
pip install /path/to../dist/arrayDesignOpt-0.1-py3-none-any.whl --force-reinstall (optional) <br />

------------------------------------------------------- <br />
Simple script <br />
from arrayDesignOpt import main  <br />
arr = main.AntennaArray(N=64) <br />
norm_pattern, _domg, _cosxy, _radmin, _pat_int  = arr.getPattern() <br />
SLLdB, XBest, YBest, BeamOpt, Xlast1, Ylast1 = arr.minimizeSLL( _radmin,gain=10, plotPattern=True) <br />
