import Sparta
from sim_creator import SimConfig
import numpy as np



simulator = Sparta.Sim("/home/elyawy/Development/Msc/sandbox/seqal_chrI_87389-87501_+_YAL030W.aln/RAxML_tree.tree")



sim_config = SimConfig(conf_file="/home/elyawy/Development/Msc/projects/Sparta/SpartaPipeline/test.conf", seq_lengths=[100,200])

SEED = 1#uuid.uuid4()
# SEED = hex(randint(0, 255)) if SEED == -1 else SEED
simulator.set_seed(SEED)

np.random.seed(SEED)

# minimum number of simulations for correction is 32, 200+ recommended.
sim_params = sim_config.get_random_sim(10000)

sum_stats_all = []
for params in sim_params:
    # print(params)
    numeric_params = [params[0],params[2][0], params[3][0], params[4], params[5]]

    simulator.init_sim(*params)

    sim_msa = simulator.run_sim()
    sum_stats_all.append(np.array(numeric_params + sim_msa.get_sum_stats()))


