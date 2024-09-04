# Which arguments are required, optional and ignored for each simulation type?

## Marta's Continuous Model
### Lattice

## Required

- `time_step`
- `final_time`
- `no_possible_phenotypes`
- `tumour_selectivity`
- `CTL_selectivity`
- `binding_affinity`
- `affinity_range`
- `tumour_natural_prolif_rate`
- `tumour_natural_death_rate`
- `tumour_interaction_induced_rate`
- `CTL_natural_prolif_rate`
- `CTL_natural_death_rate`
- `CTL_interaction_induced_rate`
- `tumour_phenotypic_variation_probability`
- `A`
- `a`

## Optional

## Ignored

- `no_init_tumour_cells`
- `no_init_CTL_cells`
- `sequence_matrix_config`
- `affinity_matrix_config`
- `CTL_sequence_path`
- `tumour_sequence_path`
- `binding_scaling`

### Sequence

## Required

- `time_step`
- `final_time`
- `no_possible_phenotypes`
- `tumour_selectivity`
- `CTL_selectivity`
- `affinity_range`
- `tumour_natural_prolif_rate`
- `tumour_natural_death_rate`
- `tumour_interaction_induced_rate`
- `CTL_natural_prolif_rate`
- `CTL_natural_death_rate`
- `CTL_interaction_induced_rate`
- `tumour_phenotypic_variation_probability`
- `A`
- `a`
- `affinity_matrix_config`

## Optional


## Ignored

- `binding_affinity`
- `no_init_tumour_cells`
- `no_init_CTL_cells`
- `sequence_matrix_config`
- `CTL_sequence_path`
- `tumour_sequence_path`
- `binding_scaling`

## Discrete Model

### Lattice


## Required

- `time_step`
- `final_time`
- `no_possible_phenotypes`
- `tumour_selectivity`
- `CTL_selectivity`
- `binding_affinity`
- `affinity_range`
- `no_init_tumour_cells`
- `no_init_CTL_cells`
- `tumour_natural_prolif_rate`
- `tumour_natural_death_rate`
- `tumour_interaction_induced_rate`
- `CTL_natural_prolif_rate`
- `CTL_natural_death_rate`
- `CTL_interaction_induced_rate`
- `tumour_phenotypic_variation_probability`

## Optional


## Ignored

- `A`
- `a`
- `sequence_matrix_config`
- `affinity_matrix_config`
- `CTL_sequence_path`
- `tumour_sequence_path`
- `binding_scaling`

### Sequence

## Required

- `time_step`
- `final_time`
- `no_init_tumour_cells`
- `no_init_CTL_cells`
- `tumour_natural_prolif_rate`
- `tumour_natural_death_rate`
- `tumour_interaction_induced_rate`
- `CTL_natural_prolif_rate`
- `CTL_natural_death_rate`
- `CTL_interaction_induced_rate`
- `tumour_phenotypic_variation_probability`
- `sequence_matrix_config`
- `affinity_matrix_config`
- `CTL_sequence_path`
- `tumour_sequence_path`
- `binding_scaling`

## Optional

## Ignored

- `no_possible_phenotypes`
- `tumour_selectivity`
- `CTL_selectivity`
- `binding_affinity`
- `affinity_range`
- `A`
- `a`