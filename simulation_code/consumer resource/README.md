# Consumer-resource (MiCRM) generator

Generates the microbial consumer-resource model (MiCRM) simulations consumed by
[`simulations/consumer_resource`](../../simulations/consumer_resource).

## Upstream / attribution

The core MiCRM library in [`lib/`](lib/) — `run_mcrm.m`, `mcrm_params.m`,
`makePhyloConsumers.m`, `getConsumerPriors.m`, `getMetabolism.m`,
`getRandomConsumerMatrix.m`, `drchrnd.m`, `coarseGrainCommunityStructure.m`, etc. — is
**based on the `mcrm` (microbial consumer-resource model) library by J. Goldford**:

> https://github.com/jgoldford/mcrm

Please credit that repository when using this generator. `gen_data_mcrm.m` and
`example_mcrm.m` are the drivers written on top of that library to produce the
parameter/init files (`saved_sims/I*/params.mat`, `inits.mat`) for this project.

## Files
- `gen_data_mcrm.m` — main driver: samples a community and runs the MiCRM sweep.
- `example_mcrm.m` / `example.m` — minimal usage examples.
- `check_reproducibility.m` — reproducibility check.
- `lib/` — the MiCRM routines (see upstream above).
