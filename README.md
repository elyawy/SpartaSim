# Sparta INDEL Simulator

A simple and fast simulator written in c++ with python bindings, Allowing advanced research in evolutionary biology.

## ⚡️ Quick start

- To compile the simulator c++ python package do the following:
  1. Clone this repo to a linux machine:
     ```bash
     git clone --recursive https://github.com/elyawy/SpartaSim.git
     ```
  3. Create a python 3.6+ environment.
  4. cd into the SpartaSim folder.
  5.  In the terminal run ```pip install .```.
  6.  For ease of use, it is recommended to install the module at: https://github.com/elyawy/SpartaPipeline/tree/main/utils
   
- Example use of the simulator:
  
```python
from elyawy.sparta import Simulator, Msa

tree_file = "path_to_newick_tree_file"
sim = Simulator(tree_file)

sim_params = {
    "root_length": 500,
    "insertion_rate": 0.03,
    "deletion_rate": 0.04,
    "length_distribution": "zipf",
    "length_parameter_insertion": 1.5,
    "length_parameter_deletion": 1.5
}
sim_params = list(sim_params.values())

sim.init_sim(*sim_params)

msa = sim()
```

- Generating summary statistics:
  
```python
msa_stats = msa.get_sum_stats()
# or with input file MSA:
msa_file = "path_to_fasta_msa_file"
msa = Msa(msa_file)
msa_stats = msa.get_sum_stats()
```

