# kalman_pll_testbench

This repository contains the Kalman PLL Testbench project for MATLAB R2024b. It is designed to simulate and analyze GNSS receiver dynamics using a Kalman filter-based Phase Lock Loop (PLL). The project integrates several components—including scintillation models and auxiliary libraries—some of which are managed as Git submodules.

## Getting Started

### Cloning the Repository

To ensure that all required submodules are properly cloned, you can use one of the following methods:

#### Option 1: Clone with Submodules Automatically

Clone the repository and its submodules in one step:

```bash
git clone --recurse-submodules <repository_url>
```

#### Option 2: Clone and Then Initialize Submodules
If you've already cloned the repository without submodules, navigate to the repository directory and run:

```bash
git submodule update --init --recursive
```

### Updating Submodules
To update the submodules to the latest commit on their tracked branches, run:

```bash
git submodule update --remote
```

## Repository Structure
Below is a simplified view of the repository structure:

```bash
kalman_pll_testbench/
├── .git/                      # Git repository and configuration files
├── libs/                      # Library folders (some managed as submodules)
│   ├── arfit/                 # ARFIT library (submodule)
│   ├── auxiliary_codes_old/   # Legacy auxiliary codes
│   ├── get_received_signal_functions/  # Signal generation functions
│   ├── kalman_pll/            # Kalman PLL implementation
│   ├── plots/                 # Plotting utilities
│   └── scintillation_models/  # Ionospheric scintillation models (submodules)
│       ├── cornell_scintillation_model/
│       └── refactored_tppsm/
├── scripts/                   # Additional scripts for simulations and tests
└── README.md                  # This file
```
*Note*: The actual file structure includes many internal Git folders (e.g., .git/modules/) and additional files. The above tree represents the primary directories you will interact with.

## Project Requirements
- MATLAB R2024b or newer: The project is developed and tested using MATLAB R2024b.
- Git: Used for version control and managing submodules.

## Contributing
Contributions are welcome! To contribute:

1. Fork the repository
2. Create a new branch for your feature or bufix
3. Submit a pull request with your changes.

For any issues or questions, please open an issue in the GitHub repository.

## License
This project is licensed under the [Your License Name Here] License.

## Additional Information
For further details on the project's design, usage, and configuration, please refer to the documentation in the /docs folder or contact the project maintainers.

## Authors information:

### Rodrigo de Lima Florindo
- Orcid: https://orcid.org/0000-0003-0412-5583
- Email: rdlfresearch@gmail.com


