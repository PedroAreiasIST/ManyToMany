// MeshConstants.cs - Shared constants for mesh library
// Eliminates duplicate EPSILON definitions across all files
// License: GPLv3

namespace Numerical;

/// <summary>
///     Shared constants and tolerance values for mesh operations.
///     All mesh library classes should use these values for consistency.
/// </summary>
public static class MeshConstants
{
    /// <summary>
    ///     Default numerical tolerance for geometric comparisons.
    ///     Used for: degenerate element detection, point coincidence, etc.
    /// </summary>
    public const double Epsilon = 1e-10;

    /// <summary>
    ///     Tolerance for detecting degenerate triangle areas.
    /// </summary>
    public const double DegenerateAreaTolerance = 1e-14;

    /// <summary>
    ///     Tolerance for detecting degenerate tetrahedron volumes.
    /// </summary>
    public const double DegenerateVolumeTolerance = 1e-15;

    /// <summary>
    ///     Default tolerance for node merging operations.
    /// </summary>
    public const double NodeMergeTolerance = 1e-12;

    /// <summary>
    ///     Small perturbation factor for grid point generation.
    /// </summary>
    public const double GridPerturbationFactor = 0.15;

    /// <summary>
    ///     Hexagonal grid row spacing factor (sqrt(3)/2).
    /// </summary>
    public const double HexRowSpacing = 0.866;
}
