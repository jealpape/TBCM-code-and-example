# Body Cover Model

## Description
This folder contains implementations related to the body cover model used in vocal fold simulations.

## Content
- `AeroPressure2DrivingForces.m`: calculation of aerodynamic forces
- `BodyCoverModel.m`: BCM model object
- `CollisionForces.m`: calculation of collision forces
- `DampingForces.m`: calculation of damping forces
- `ElasticForces.m`: calculation of elastic forces
- `scalingVocalFold.m`: [under construction]
- `setSimulationParameter.m`: [under construction]
- `Simulate.m`: Numerical solution of the equations of motion of the masses

## References
- Story, B. H., & Titze, I. R. (1995). Voice simulation with a body‐cover model of the vocal folds. The Journal of the Acoustical Society of America, 97(2), 1249-1260.

---

# Muscle Activation

## Description
This folder contains the muscle activation model used in simulations.

## Content
- `MuscleActivation.m`: Muscle activation object
- `Rule2BodyCoverParameters.m`: Implementation of Titze2002 rules

## References
- Titze, I. R., & Story, B. H. (2002). Rules for controlling low-dimensional vocal fold models with muscle activation. The Journal of the Acoustical Society of America, 112(3), 1064-1076.
---

# Muscle Control Model

## Description
This folder includes the muscle control model for vocal fold simulations.

## Content
- `AuxMaterial/`: Folder with auxiliary materials.
- `OldVersion/`: Previous versions of the model.
- `ArithAngleCollision.m`: [under construction]
- `CalcBodyCoverParameters.m`: New implementation of Titze2002
- `ColReaction.m`: [under construction]
- `MuscleControlModel.m`: Motor control object
- `plotCAJointSpringConstants.m`: CAJ plotting code
- `plotCTJointSpringConstants.m`: CTJ plotting code
- `SimulateCAJ.m`: Numerical solution of the CAJ motion equation
- `SimulateCTJ.m`: Numerical solution of the CTJ motion equation
- `SimulatePosture.m`: Solution of the differential equations for posture variables
- `TorqueCol.m`: [under construction]

## References
- Titze, I. R. The Myoelastic Aerodynamic Theory of Phonation, 1st edition. National Center for Voice and Speech, 2006.
- Titze, I. R., & Hunter, E. J. (2007). A two-dimensional biomechanical model of vocal fold posturing. The Journal of the Acoustical Society of America, 121(4), 2254-2260.
- Titze, I. R., & Story, B. H. (2002). Rules for controlling low-dimensional vocal fold models with muscle activation. The Journal of the Acoustical Society of America, 112(3), 1064-1076.

---

# Subglottal Tract Model

## Description
This folder contains implementations related to the subglottal tract model. The object includes methods to compute total pressure and airflow for each segment of the tract, for example.

## Content
- `MatAux/`: Folder with auxiliary materials.
- `getStateSpaceWRAModelmodeA.m`: Obtain state-space representation of WRA, excluding PL
- `getStateSpaceWRAModelmodeB.m`: Obtain state-space representation of WRA, including PL
- `getSubglottalTract.m`: Functions for subglottal tract area
- `plotTract.m`: [under construction]
- `setSimulationParameter.m`: Initializer
- `Simulate.m`: [under construction]
- `SimulateWRAmodeA.m`: Solution of WRA equations, excluding PL
- `SimulateWRAmodeB.m`: Solution of WRA equations, including PL
- `SubglottalTractModel.m`: Subglottal tract object

## References
- Story, B. H. (1995). Physiologically-based speech simulation using an enhanced wave-reflection model of the vocal tract. the University of Iowa.

---

# Vocal Tract Model

## Description
This folder contains the vocal tract model, including specific configurations. The object includes methods to compute total pressure and airflow for each tract segment, for example.

## Content
- `MatAux/`: [Folder with auxiliary materials.]
- `getFemaleVocalTract_Story1998.m`: Functions for female tract area
- `getMaleVocalTract_Story1996.m`: Functions for male tract area
- `getMaleVocalTract_Story2008.m`: Functions for male tract area
- `getSimpleVocalTract.m`: Functions for base tract area
- `getStateSpaceWRAModel.m`: Obtain state-space representation of WRA
- `plotTract.m`: [under construction]
- `setSimulationParameter.m`: Initializer
- `SimulateSupraTractResponse.m`: [under construction]
- `SimulateWRA.m`: Solution
- `VocalTractModel.m`: Tract object

## References
- Story, B. H. (1995). Physiologically-based speech simulation using an enhanced wave-reflection model of the vocal tract. the University of Iowa.

---

# Triangular Body Cover Model

## Description
This folder includes implementations of the triangular body cover model.

## Content
- `Documentation/`: Folder with model documentation.
- `OldVersion/`: Previous versions of the model.
- `AeroPressure2DrivingForces.m`: Aerodynamic force
- `CollisionForces.m`: Collision force
- `DampingForces.m`: Damping force
- `ElasticForces.m`: Elastic force
- `scalingVocalFold.m`: [under construction]
- `setSimulationParameter.m`: [under construction]
- `Simulate.m`: Solution of the motion equations of the 3 masses in the model
- `TriangularBodyCoverModel.m`: TBCM mass object
- `zzCollisionForces.m`: [under construction]

## Usage Example
```matlab
% Example code to use the triangular body cover model
```

## References
- Galindo, G. E., Peterson, S. D., Erath, B. D., Castro, C., Hillman, R. E., & Zañartu, M. (2017). Modeling the pathophysiology of phonotraumatic vocal hyperfunction with a triangular glottal model of the vocal folds. Journal of Speech, Language, and Hearing Research, 60(9), 2452-2471.
- Alzamendi, G. A., Peterson, S. D., Erath, B. D., Hillman, R. E., & Zañartu, M. (2022). Triangular body-cover model of the vocal folds with coordinated activation of the five intrinsic laryngeal muscles. The Journal of the Acoustical Society of America, 151(1), 17-30.

---

# Intrinsic Muscles

## Description
This folder includes individual models of intrinsic muscles related to phonation, with constants for stress and strain calculations.

## Content
- `CricothyroidModel.m`: [under construction]
- `InterarytenoidModel.m`: [under construction]
- `LateralCricoarytenoidModel.m`: [under construction]
- `LigamentModel.m`: [under construction]
- `MucosaModel.m`: [under construction]
- `Muscle1DModel.m`: [under construction]
- `PosteriorCricoarytenoidModel.m`: [under construction]
- `ThyroarytenoidModel.m`: [under construction]

## Usage Example
```matlab
% Example code to use intrinsic muscle models
```

## References
- Titze, I. R. The Myoelastic Aerodynamic Theory of Phonation, 1st edition. National Center for Voice and Speech, 2006.