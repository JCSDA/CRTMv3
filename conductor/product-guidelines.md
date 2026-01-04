# Product Guidelines - CRTM v3.2.0

## Documentation & Code Comments
- **Minimalist Approach:** Code comments should prioritize explaining the "why" behind complex logic. The code itself should be written for clarity to minimize the need for descriptive "what" comments.

## Brand Messaging & Scientific Integrity
- **Standardized Terminology:** All communications and documentation must use clear, consistent terminology for physical quantities and model versions, adhering to the CRTM official documentation.
- **Authoritative & Collaborative:** CRTM is a community standard for operational forecasting. We emphasize its community-driven, open-source nature and actively encourage collaborative contributions.

## Coding Standards
- **JCSDA/JEDI Conventions:** Adhere strictly to the JCSDA/JEDI coding conventions to ensure seamless integration with the JEDI framework.
- **Strict Variable Declaration:** The use of `IMPLICIT NONE` is mandatory in all modules and subroutines to prevent implicit typing bugs.
- **Modern Fortran (Encouraged):** While Fortran 2008+ standards are encouraged for better maintainability, they are not a hard requirement to ensure compatibility for users with older compiler versions.

## Testing & Validation
- **Mandatory Regression Testing:** All modifications must pass the comprehensive regression test suite (e.g., `test_Simple`, `test_ClearSky`) to prevent scientific or performance regressions.
- **Unit Testing Requirement:** Each new module or major subroutine must be accompanied by a unit test located in `test/mains/unit`.