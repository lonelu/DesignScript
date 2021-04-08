# reload module in ipython
import importlib 

importlib.reload(sys.modules['ligand_database']) 
from ligand_database import * 