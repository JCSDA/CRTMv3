# Product Guidelines - CRTMv3 Missing Jacobians

## Development Principles
- **Clarity through Code:** Prioritize clean, self-explanatory variable naming and modular code structure. Rely on the code itself to convey meaning, keeping inline comments focused on "why" rather than "what" when necessary.
- **Mathematical Robustness:** Implement explicit checks and "safe" fallbacks for potential mathematical singularities, such as divisions by zero, ensuring the model remains stable across a wide range of input conditions.
- **Pragmatic Refactoring:** Balance the need for stability with the goal of improvement. Refactor existing code when it provides clear benefits to the clarity, maintainability, or performance of the new Jacobian implementations.

## Implementation Standards
- **Mathematical Integrity:** Ensure that every Tangent Linear and Adjoint implementation is a rigorous derivation of the corresponding Forward model.
- **Consistency:** Maintain a consistent style and structure across all new modules, following the established patterns in CRTMv3.
