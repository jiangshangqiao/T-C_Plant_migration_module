# T-C_Plant_migration_module
Source code for the plant migration module in the T&C model.

## Overview

The ** Plant migration module** extends the T&C model by simulating the spatial migration of plant species in response to environmental changes. The module represents plant dispersal, establishment, and colonization processes, enabling dynamic shifts in species distributions over time.

This repository contains the source code required to integrate plant migration processes into the T&C model.

## Feature

- Simulate plant migration across spatially explicit landscapes.
- Models seed dispersal and colonization of new grid cells.
- Incorporate environmental constraints on seedling establishment.
- Compatible with the T&C vegetation dynamic framework.
- Designed for climate change and ecosystem response studies.

---

## Repository Structure

```text
├── src/             # Source code
├── input_example/   # Example of input data and parameters
├── input_ori_example/         # Example of input data and parameters of the original version
├── docs/            # Documentation
└── README.md
```

## Requirments

- T&C model
- MATLAB
- Required T&C input datasets

---

## Installation

Clone the repository 

```bash
git clone https://github.com/jiangshangqiao/T-C_Plant_migration_module.git
```

Integrate the module into the T&C model according to your project workflow.

---

## Usage

1. Prepare climate and environmental input data.
2. Configure plant migration parameters.
3. Run the T&C model.
4. Analyze vegetation dynamics and migration results.

---

## Model Processes

The module includes:

- Seed production in the original T&C model as the fruit and flower carbon pool
- Seed dispersal
- Colonization of suitable habitats
- Seedling establishment
- Dynamic vegetation updates

## Inputs

Typical inputs include:

- Initial vegetation distribution
- Climate forcing
- Soil properties
- Topography
- Vegetation parameters
- Prescribed species priority
- Seed dispersal parameters
- Compared with the run file of the original T&C model, the extended model's run file includes an additional section for seed dispersal parameters. 

P.S. I have included both the original and the updated versions of the run and initialization files for comparison.

---

## Outputs

Typical outputs include:

- Hydrological fluxes
- Carbon fluxes
- Species distribution maps (variable: Ccrown_t)
- Time series of vegetation changes
- The current configuration outputs all variables throughout the entire simulation, which can generate a large volume of data. To reduce the output size, modify the OUTPUT_MANAGER_PAR1.m file in the original T&C model so that only the required variables and results are saved. 
---

## Citation

If you use this module in your research, please cite the relevant **T&C model publication** and any publication describing this plant migration module.

---

## Contributing

Contributions are welcome through GitHub Issues and Pull Requests.

---

## License

This project is licensed under the MIT License. See the LICENSE file for details.

This repository contains an extension module for the T&C model. Users should ensure they comply with the licensing terms of the original T&C model when integrating this module.

---

## Contact

For questions, suggestions, or bug reports, please open an Issue or contact the repository maintainer.

---

## Acknowledgments

This module builds upon the **T&C** ecohydrological model. We gratefully acknowledge the contributions of the T&C development team and the ecological modeling community.








