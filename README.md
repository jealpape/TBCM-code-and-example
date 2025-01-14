# TBCM Code

## Description
This repository contains MATLAB code for the simulation of the Triangular Body-Cover Model (TBCM). The model is organized into folders according to its components, including subglottal tract models, muscle activation, and vocal tract.

## Project Structure
- `@BodyCoverModel`: Structure based on the BCM by Titze and Story
- `@MuscleActivation`: Muscle activation structure for the BCM (Titze's Rules)
- `@MuscleControlModel`: Structure for modeling muscle activation for the TBCM
- `@SubglottalTractModel`: Subglottal tract object introducing the WRA
- `@TriangularBodyCoverModel`: TBCM object (three masses)
- `@VocalTractModel`: Vocal tract object introducing the WRA
- `+IntrinsicMuscles`: Constants for the Kelvin model of each laryngeal muscle
- `+measure`: Set of functions to calculate aerodynamic features
- `vfsolver`: Previous version implemented as individual scripts, uses `solveFlow`
- `SimpleExample.m`: A basic example to understand how to use the code

## Example Codes

### `SimuBaseSignal.m`

This function-based script takes the following input values:
- Muscle activation
- Subglottal pressure
- Vowel
- Gender

Outputs include 50 ms of signals of interest:
- Geometry
- Subglottal pressure
- Collision pressure
- Output pressure (MIC)
- Glottal flow

### `code_gridSignal.m`

This script uses the `SimuBaseSignal.m` function within for-loops to compute signals and acoustic parameters for multiple phonation configurations.

## Acknowledgment
Special thanks to Dr. Gabriel Alzamendi for his hard work in building and implementing this code.

## References
- Galindo, G. E., Peterson, S. D., Erath, B. D., Castro, C., Hillman, R. E., & Zañartu, M. (2017). Modeling the pathophysiology of phonotraumatic vocal hyperfunction with a triangular glottal model of the vocal folds. Journal of Speech, Language, and Hearing Research, 60(9), 2452-2471.
- Alzamendi, G. A., Peterson, S. D., Erath, B. D., Hillman, R. E., & Zañartu, M. (2022). Triangular body-cover model of the vocal folds with coordinated activation of the five intrinsic laryngeal muscles. The Journal of the Acoustical Society of America, 151(1), 17-30.
