# ArrayDesignOpt
Python files developed to simulate 2D random array patterns, and transform regular uniform arrays to pseudo random through side lobe minimization.



## Instalation
Requierement libraries:
- numpy v1.19.5
- matplotlib v3.6.0
- shapely
- IPython
- scipy

git clone https://github.com/joabAM/arrayDesignOpt.git
cd arrayDesignOpt
python setup.py bdist_wheel
pip install /path/to../dist/arrayDesignOpt-0.1-py3-none-any.whl --force-reinstall (optional)

-------------------------------------------------------
Simple script
from arrayDesignOpt import main 
arr = main.AntennaArray(N=64)
norm_pattern, _domg, _cosxy, _radmin, _pat_int  = arr.getPattern()
SLLdB, XBest, YBest, BeamOpt, Xlast1, Ylast1 = arr.minimizeSLL( _radmin,gain=10, plotPattern=True)
