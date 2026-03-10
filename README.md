# IMKResiDisp – OpenSeesPy Material

This repository provides the C++ implementation of the IMKResiDisp material model, wrapped for use with **OpenSeesPy**. The model is developed based on the IMK model to more accurately predict the residual displacement of RC structures, and is implemented as a subclass of UniaxialMaterial in OpenSees.

## Files

- `IMKResiDisp.h` – Header file defining the IMKResiDisp class
- `IMKResiDisp.cpp` – Source file implementing the material model

## Installation

1. Copy `IMKResiDisp.h` and `IMKResiDisp.cpp` into the OpenSees source directory: OpenSees/SRC/material/uniaxial/
2. Recompile OpenSeesPy to include the IMKResiDisp material.

## Usage

After compiling, the material can be defined in Python as:

```python
# Example usage in OpenSeesPy
import openseespy.opensees as ops
ops.uniaxialMaterial('IMKResiDisp', MatTag, Ke, posUy_0, posUcap_0, posUu_0, posFy_0, posFcapFy_0, posFresFy_0, negUy_0, negUcap_0, negUu_0, negFy_0, negFcapFy_0, negFresFy_0, LAMBDA_S, LAMBDA_C, LAMBDA_A, c_S, c_C, c_A, alpha, Offset, beta_F, beta_U, theta_F, theta_U)
```
The material can then be used together with elements such as zeroLength to simulate concentrated plastic hinge behavior.

## Reference

If you use this code in your research, please cite the corresponding publication:
Wu, S., et al. (2026). Hysteretic behavior of RC columns under asymmetric decreasing-amplitude loading protocols and the improved model for residual displacement analysis. Journal of Building Engineering.
