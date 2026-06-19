# Contributing to ManyToMany

Thanks for your interest in improving **ManyToMany**! This guide covers how to build, test, and propose changes.

## Prerequisites

- [.NET 9.0 SDK](https://dotnet.microsoft.com/download/dotnet/9.0) (x64)
- A 64-bit OS (Windows, Linux, or macOS). The library is **x64-only by design**.
- Optional accelerators for full functionality: Intel MKL (PARDISO) and the CUDA Toolkit 11.0+. Neither is required — the library degrades gracefully to managed SIMD when they are absent.

## Building & running

```bash
# Restore and build the whole solution
dotnet build Numerical.sln -c Release

# Run the example suite (meshing + 2D/3D fracture mechanics)
dotnet run --project tests/Teste -c Release
```

The continuous-integration workflow (`.github/workflows/dotnet.yml`) restores, builds, and tests the solution on every push and pull request targeting `master`, so please make sure your branch builds cleanly before opening a PR.

## Project layout

ManyToMany is a layered stack of projects (`Relations` → `Matrices` → `Meshing` → `Nonlinear`/`Postprocess` → `Teste`). All public types share the `Numerical` namespace. See the [Project Structure](README.md#project-structure) section of the README for a file-by-file map.

When adding a feature, place it in the lowest layer that makes sense and avoid introducing upward dependencies (e.g. `Matrices` must not depend on `Meshing`).

## Coding guidelines

- **Target framework:** `net9.0`, C# `latest`, nullable reference types enabled.
- **Match the surrounding style.** Mirror the naming, formatting, and comment density of the file you are editing.
- **Performance matters.** This is numerical code: prefer `Span<T>`/`ReadOnlySpan<T>` in hot paths, reuse buffers via `ArrayPool<T>`, and avoid allocations inside inner loops. Follow the existing SIMD/parallel patterns rather than inventing new ones.
- **Keep the public surface tight.** Expose new API only when it is meant to be used externally; mark implementation helpers `internal`/`private`.
- **Document public APIs** with XML doc comments where the existing code does.

## Documentation changes

If your change adds, removes, or renames public API, please update:

1. The relevant section of the [README](README.md) — including its **Public API Reference** tables.
2. The matching deep-dive in [`Docs/`](Docs) when applicable.

The README is treated as the canonical, code-verified API reference; keep its signatures in sync with the source.

## Submitting changes

1. Create a feature branch off `master`.
2. Make focused commits with clear, descriptive messages.
3. Ensure `dotnet build Numerical.sln -c Release` succeeds.
4. Open a pull request describing **what** changed and **why**. For substantial or architectural changes, please open an issue to discuss the approach first.

## Reporting issues

When filing a bug, include:

- the module and method involved,
- a minimal code snippet that reproduces the problem,
- your OS, .NET version, and whether MKL/CUDA are installed,
- the expected vs. actual behavior (and any exception/stack trace).

## License

By contributing, you agree that your contributions will be licensed under the project's [GPLv3](LICENSE) license.
