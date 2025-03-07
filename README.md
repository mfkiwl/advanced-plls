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

## Git Hooks
The repository includes a folder named hook_examples that contains sample Git hook configurations. These hooks are designed to automate tasks such as deinitializing submodules when switching branches. Since Git does not track hook scripts by default (they reside in the .git/hooks directory), these examples are provided for your convenience.

To use a hook:

1. Copy the desired hook file (e.g., post-checkout) from the hook_examples folder into your local repository's .git/hooks directory.
2. Ensure the file is named exactly as required by Git (without any extension, e.g., rename post-checkout.sample to post-checkout).
3. Make the file executable. On Unix-like systems, you can run:
```bash
chmod +x .git/hooks/post-checkout
```
This setup can help streamline your workflow by automatically managing submodules and maintaining a clean working directory across branches.

## Contributing
Contributions are welcome! To contribute:

1. Fork the repository
2. Create a new branch for your feature or bufix
3. Submit a pull request with your changes.

For any issues or questions, please open an issue in the GitHub repository.

## Roadmap for Features Implementations

This section consolidates the current progress on the Adaptive AR/RBF modules with the planned future enhancements. The roadmap is organized into key areas covering smoothing, filtering, augmentation, adaptability, and Bayesian estimation.

*Note: The Stashed status signifies that this features will only be possibly implemented later, and it is not a priority for now.*

---

### 1. Completed Work

- **AR Parameter Estimation:**
  - **Aryule vs. ARfit:** Use the aryule method (or modify ARfit to remove the intercept vector).  
    **Status:** Done!
  - **Block AR Model:** Introduced the Block AR model parameter estimation algorithm.  
    **Status:** Done!
  - **Sliding AR Model:** Introduced the Sliding AR model parameter estimation algorithm.  
    **Status:** Done!

---

### 2. Planned Features

#### 2.1. Smoothing and Filtering
**Status:** Stashed.
- **Kalman Smoother:** Implement a Kalman Smoother.
  **Status:** Stashed.
- **Classical Batch Smoothing:** Implement a classical batch smoothing technique.
  **Status:** Stashed.

#### 2.3. New Augmentation Models to implement
- **Discrete Wiener Augmentation:** Implement an augmentation model using a second discrete Wiener model.
- **ARIMA Augmentation:** Implement an ARIMA augmentation model.

#### 2.2. Adaptive Augmentation Models
- **Dual Kalman Filter:** Introduce the Dual Kalman filter algorithm from Frederieke.  
  *Note:* Avoid using the error state (to prevent major refactoring). 
  **Status:** Stashed. (We've observed under the development of online learning AR models that it did not helped the single-frequency KF to not track the refractive + diffractive phase fluctuations. In addition, when the models are only trained with the diffractive phase, they present unmeaningful estimates of diffractive phase under refractive + diffractive phase scenario.)
- **RBF Network:** Introduce an RBF model that can be configured with any number of neurons.  
  **Status:** Stashed.

#### 2.4. Adaptability and Covariance Estimation
- **Robust Covariance Estimation:** Develop robust techniques for estimating state and measurement covariances.
  **Status:** Under Development.
- **Reinforcement Learning for Covariances:** Implement a reinforcement learning network technique to adapt state and measurement covariances.
  **Status:** Under Development.

#### 2.5. Bayesian Estimation Techniques
- **Extended Kalman Filter:** Implement the Extended Kalman Filter.
  **Status:** Under Development.
- **Unscented Kalman Filter:** Implement the Unscented Kalman Filter.
  **Status:** Under Development.
- **Cubature Kalman Filter:** Implement the Cubature Kalman Filter.
  **Status:** Under Development.
- **Particle Filters:** Explore Particle filters
  **Status:** Stashed.

## 2.6. Multi-frequency Models
---
- **Multi-Frequency Carrier Tracking:** Explore multi-frequency carrier phase tracking after addressing cycle-slip avoidance on single-frequency carrier phase tracking.
  **Status:** Stashed. (I think it'd be interesting to first stress out the single-frequency tracking scenario as much as possible before going into this.)

## License
This project is licensed under the [Your License Name Here] License.

## Additional Information
For further details on the project's design, usage, and configuration, please refer to the documentation in the /docs folder or contact the project maintainers.

## Authors' Information:

### Rodrigo de Lima Florindo
- Orcid: https://orcid.org/0000-0003-0412-5583
- Email: rdlfresearch@gmail.com


