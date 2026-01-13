# Plan - Aerosol Scattering Jacobians

## Phase 1: Analysis and Infrastructure
- [ ] Task: Audit `src/AtmScatter/AerosolScatter` to identify specific missing Jacobian code blocks.
- [ ] Task: Set up a baseline unit test environment for Aerosol Forward model outputs.
- [ ] Task: Conductor - User Manual Verification 'Phase 1: Analysis and Infrastructure' (Protocol in workflow.md)

## Phase 2: Implementation of TL and AD Models
- [ ] Task: Implement Tangent Linear (TL) updates for Aerosol scattering.
- [ ] Task: Implement Adjoint (AD) updates for Aerosol scattering.
- [ ] Task: Write tests to verify TL/AD consistency using the gradient test or finite differences.
- [ ] Task: Conductor - User Manual Verification 'Phase 2: Implementation of TL and AD Models' (Protocol in workflow.md)

## Phase 3: K-Matrix Integration and Validation
- [ ] Task: Develop the K-Matrix interface for the Aerosol modules.
- [ ] Task: Integrate the new Jacobians into the top-level `CRTM_K_Matrix_Module`.
- [ ] Task: Perform final validation against expected sensitivity profiles.
- [ ] Task: Conductor - User Manual Verification 'Phase 3: K-Matrix Integration and Validation' (Protocol in workflow.md)
