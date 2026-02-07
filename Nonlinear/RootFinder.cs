using System.Buffers;
using System.Runtime.CompilerServices;

namespace Numerical;

/// <summary>
///     Advanced root finding algorithms with overloaded methods for different function types.
///     All methods are thread-safe as they contain no static mutable state.
/// </summary>
public static class RootFinder
{
    /// <summary>
    ///     Status codes returned by root finding algorithms
    /// </summary>
    public enum Status
    {
        OK = 0, // Converged to root within function tolerance
        Tolerance = 1, // Converged to bracket width tolerance (may not meet function tolerance)
        MaxIterations = 2, // Maximum iterations reached without convergence
        NoBracket = 3, // Function doesn't have opposite signs at endpoints
        BadInput = 4, // Invalid input (NaN or infinite values)
        NonFinite = 5, // Function returned non-finite value
        TooNarrow = 6 // Initial interval too narrow for computation
    }

    // ===== Tolerance and Convergence Parameters =====

    private const double FTOL = 1e-10; // Absolute function tolerance for convergence
    private const double RTOL = 1e-8; // Relative tolerance for interval width
    private const double ATOL = 1e-12; // Absolute tolerance for interval width
    private const int MAX_ITER = 100; // Maximum iterations before giving up
    private const double EPS_MACHINE = 2.2204460492503131e-16; // IEEE 754 double precision epsilon

    // ===== Algorithm-Specific Tuning Parameters =====

    // Inverse Quadratic Interpolation (IQI)
    private const double EPS_IQI = 1000.0 * EPS_MACHINE; // Denominator safety threshold for IQI
    private const double IQI_STEP_LIMIT_FACTOR = 3.0; // Maximum IQI step relative to interval

    // Newton-Raphson
    private const double NEWTON_STEP_SCALE_FACTOR = 0.5; // Newton step safety scaling
    private const double DERIV_STEP_CHECK_FACTOR = 0.5; // Factor for derivative quality check

    // Adaptive tolerance
    private const double FTOL_SCALE_FACTOR = 100.0; // Scale factor for adaptive tolerance
    private const double MAX_FTOL_SCALE = 100.0; // Cap for adaptive tolerance scaling

    // Recovery and robustness
    private const double GOLDEN_RATIO_COMPLEMENT = 0.3819660112501051; // (3 - √5) / 2 for golden section
    private const double DIVERGENCE_MULTIPLIER = 2.0; // Multiplier to detect divergence
    private const int MIN_ITER_DIVERGENCE_CHECK = 5; // Iterations before checking divergence
    internal const double BACKTRACK_FACTOR = 0.5; // Backtracking reduction factor

    // ITP Algorithm (Interpolate-Truncate-Project)
    private const int N0_ITP = 1; // ITP parameter: additional iterations over bisection
    private const double KTR_ITP = 0.2; // ITP parameter: truncation factor κ ∈ (0, ∞)
    private const double PEXP_ITP = 2.0; // ITP parameter: super-linear convergence exponent
    private const int N_SAMPLE_ITP = 7; // Number of interior points to sample for smart start

    // ===== Public API =====

    /// <summary>
    ///     Find root using derivative information (hybrid Newton-Raphson + IQI algorithm).
    ///     For functions that return both f(x) and f'(x).
    ///     Typical convergence: O(log ε) for smooth functions with Newton steps.
    /// </summary>
    /// <param name="xmin">Lower bound of search interval</param>
    /// <param name="xmax">Upper bound of search interval</param>
    /// <param name="func">Function returning (f(x), f'(x))</param>
    /// <returns>Tuple of (root, status)</returns>
    public static (double xsol, Status status) FindRoot(
        double xmin,
        double xmax,
        Func<double, (double f, double df)> func)
    {
        return RootG(xmin, xmax, func);
    }

    /// <summary>
    ///     Find root without derivative information (ITP algorithm).
    ///     For functions that only return f(x).
    ///     Optimal worst-case convergence: O(log log(1/ε)) iterations.
    /// </summary>
    /// <param name="xmin">Lower bound of search interval</param>
    /// <param name="xmax">Upper bound of search interval</param>
    /// <param name="func">Function returning f(x)</param>
    /// <returns>Tuple of (root, status)</returns>
    public static (double xsol, Status status) FindRoot(
        double xmin,
        double xmax,
        Func<double, double> func)
    {
        return RootITP(xmin, xmax, func);
    }

    // ===== Main Algorithms =====

    /// <summary>
    ///     Root finding with gradient (derivative) information.
    ///     Hybrid algorithm combining:
    ///     1. Newton-Raphson when derivative is reliable
    ///     2. Inverse Quadratic Interpolation (IQI) using three points
    ///     3. Secant method as fallback
    ///     4. Bisection as last resort
    ///     Always maintains a bracketing interval [a,b] with f(a)·f(b) &lt; 0.
    /// </summary>
    private static (double xsol, Status status) RootG(
        double xmin,
        double xmax,
        Func<double, (double f, double df)> func)
    {
        // ===== Input Validation =====

        if (!IsFinite(xmin) || !IsFinite(xmax))
            return (xmin, Status.BadInput);

        // Check if interval is wider than machine precision
        if (Math.Abs(xmax - xmin) <= Spacing(Math.Max(Math.Abs(xmin), Math.Abs(xmax))))
            return (xmin, Status.TooNarrow);

        // Normalize: ensure a < b
        double a = xmin, b = xmax;
        if (b < a) (a, b) = (b, a);

        // ===== Evaluate Endpoints =====

        var (fa, dfa, okA) = EvalFuncWithGradient(a, func);
        if (!okA)
        {
            // Try to recover by evaluating nearby point
            (fa, dfa, okA) = TryRecoverWithGradient(a, b, func);
            if (!okA) return (a, Status.NonFinite);
        }

        var (fb, dfb, okB) = EvalFuncWithGradient(b, func);
        if (!okB)
        {
            (fb, dfb, okB) = TryRecoverWithGradient(b, a, func);
            if (!okB) return (b, Status.NonFinite);
        }

        // Compute adaptive tolerance based on function scale
        // This prevents premature convergence for large-magnitude functions
        var fScaleInit = Math.Max(Math.Abs(fa), Math.Max(Math.Abs(fb), 1.0));
        var ftolAdaptive = FTOL * Math.Min(MAX_FTOL_SCALE, 1.0 + fScaleInit / FTOL_SCALE_FACTOR);

        // ===== Check if Endpoints are Already Roots =====

        if (Math.Abs(fa) <= ftolAdaptive)
            return (a, Status.OK);
        if (Math.Abs(fb) <= ftolAdaptive)
            return (b, Status.OK);

        // ===== Verify Bracketing =====

        if (!OppositeSign(fa, fb))
        {
            // No sign change - return point with smallest function value
            var root = Math.Abs(fa) <= Math.Abs(fb) ? a : b;
            return (root, Status.NoBracket);
        }

        // ===== Initialize Algorithm State =====

        // Maintain: |fb| >= |fc| (b is best point so far)
        if (Math.Abs(fa) < Math.Abs(fb))
            Swap3(ref a, ref b, ref fa, ref fb, ref dfa, ref dfb);

        // c stores the previous iterate or opposite bracket point
        var c = a;
        var fc = fa;
        var dfc = dfa;

        var prevStep = b - a; // Previous step size for step limiting
        var nonfinHits = 0; // Count consecutive non-finite evaluations

        // ===== Main Iteration Loop =====

        for (var iter = 1; iter <= MAX_ITER; iter++)
        {
            // Ensure b has smaller function magnitude than c
            if (Math.Abs(fc) < Math.Abs(fb))
                Swap3(ref b, ref c, ref fb, ref fc, ref dfb, ref dfc);

            var m = 0.5 * (c - b); // Midpoint displacement from b
            var width = Math.Abs(c - b); // Current bracket width
            var tolx = UlpScale(b); // Position-dependent tolerance

            // Update adaptive tolerance based on current function value
            ftolAdaptive = FTOL * Math.Min(MAX_FTOL_SCALE, Math.Max(1.0, Math.Abs(fb) / fScaleInit));

            // ===== Check Convergence =====

            // Converged if interval is too narrow or midpoint step is below tolerance
            if (Math.Abs(m) <= tolx || width <= Spacing(Math.Max(Math.Abs(b), Math.Abs(c))))
                return (b, Math.Abs(fb) <= ftolAdaptive ? Status.OK : Status.Tolerance);

            // ===== Check for Divergence =====

            // If function magnitude grows significantly, we may not have a valid bracket
            if (iter > MIN_ITER_DIVERGENCE_CHECK && Math.Abs(fb) > DIVERGENCE_MULTIPLIER * fScaleInit)
                return (b, Status.NoBracket);

            // ===== Choose Next Iterate =====

            var scale = 0.5 * Math.Abs(m); // Maximum allowed step (half interval)
            var useNewton = DerivGood(fb, dfb, scale);
            var newStep = 0.0;

            if (useNewton)
            {
                // ===== Newton-Raphson Step =====
                // Compute: x_new = x - f(x)/f'(x), limited to bracket

                newStep = -fb / dfb;
                newStep = Math.Max(-scale, Math.Min(scale, newStep));
            }
            else
            {
                // ===== Inverse Interpolation Steps =====

                // Safety thresholds to prevent division by near-zero differences
                var epsDenom = EPS_IQI * Math.Max(Math.Abs(fa), Math.Max(Math.Abs(fb), Math.Max(Math.Abs(fc), 1.0)));
                var epsValue = EPS_MACHINE *
                               Math.Max(Math.Abs(fa), Math.Max(Math.Abs(fb), Math.Max(Math.Abs(fc), 1.0)));

                // Check if we can use Inverse Quadratic Interpolation (IQI)
                // Requires three distinct function values and non-zero function magnitudes
                if (Math.Abs(fa - fb) > epsDenom &&
                    Math.Abs(fc - fb) > epsDenom &&
                    Math.Abs(fa - fc) > epsDenom &&
                    Math.Abs(fa) > epsValue &&
                    Math.Abs(fc) > epsValue)
                {
                    // ===== Inverse Quadratic Interpolation =====
                    // Interpolate inverse function: find x where p(f) = x passes through (fa,a), (fb,b), (fc,c)
                    // Then evaluate at f=0

                    var r = fb / fc;
                    var q = fb / fa;
                    var p = r * (r - q) * (b - a) + (1.0 - r) * q * (b - c);
                    var qDenom = (r - 1.0) * (q - 1.0) * (r - q);

                    if (Math.Abs(qDenom) > epsDenom)
                        newStep = p / qDenom;
                    else
                        newStep = m; // Fallback to bisection if denominator too small
                }
                else if (Math.Abs(fb - fa) > epsDenom && Math.Abs(fa) > epsValue)
                {
                    // ===== Secant Method =====
                    // Linear interpolation between two points

                    newStep = -fb * (b - a) / (fb - fa);
                }
                else
                {
                    // ===== Bisection =====
                    // Last resort when interpolation is unreliable

                    newStep = m;
                }

                // ===== Limit Interpolation Step Size =====
                // Prevent extrapolation too far outside bracket or taking steps larger than previous

                if ((newStep > 0 && newStep > IQI_STEP_LIMIT_FACTOR * m) ||
                    (newStep < 0 && newStep < IQI_STEP_LIMIT_FACTOR * m) ||
                    Math.Abs(newStep) > NEWTON_STEP_SCALE_FACTOR * Math.Abs(prevStep))
                    newStep = m; // Revert to bisection
            }

            // ===== Enforce Minimum Step Size =====
            // Prevent stalling due to roundoff by ensuring step is at least ULP

            EnforceMinStep(b, m >= 0 ? 1.0 : -1.0, ref newStep);

            // ===== Compute Trial Point =====

            var u = b + newStep;

            // Ensure u is at least tolx away from b
            if (Math.Abs(u - b) < tolx)
                u = b + Math.Sign(m) * tolx;

            // Ensure u is at least tolx away from c (prevent bracket collapse)
            if (Math.Abs(u - c) < tolx)
                u = c - Math.Sign(c - b) * tolx;

            // ===== Evaluate Function at Trial Point =====

            var (fu, dfu, okU) = EvalFuncWithGradient(u, func);
            if (!okU)
            {
                nonfinHits++;
                (fu, dfu, okU) = TryRecoverWithGradient(u, b, func);

                if (!okU)
                {
                    // After multiple failures, give up and return best known point
                    if (nonfinHits >= 3)
                    {
                        var bestPoint = Math.Abs(fb) <= Math.Abs(fc) ? b : c;
                        return (bestPoint, Status.NonFinite);
                    }

                    // Try bisection as recovery strategy
                    newStep = m;
                    EnforceMinStep(b, m >= 0 ? 1.0 : -1.0, ref newStep);
                    u = b + newStep;

                    if (Math.Abs(u - c) < tolx)
                        u = c - Math.Sign(c - b) * tolx;
                    if (Math.Abs(u - b) < tolx)
                        u = b + Math.Sign(m) * tolx;

                    (fu, dfu, okU) = EvalFuncWithGradient(u, func);
                    if (!okU)
                    {
                        var bestPoint = Math.Abs(fb) <= Math.Abs(fc) ? b : c;
                        return (bestPoint, Status.NonFinite);
                    }
                }
            }

            prevStep = newStep;

            // ===== Update Bracket =====
            // Maintain invariant: f(b) and f(c) have opposite signs

            if (OppositeSign(fb, fu))
            {
                // u and b bracket the root; move c to old b position
                a = b;
                fa = fb;
                dfa = dfb;
                c = u;
                fc = fu;
                dfc = dfu;
            }
            else
            {
                // u and c bracket the root; move a to u
                a = u;
                fa = fu;
                dfa = dfu;
            }

            // ===== Update Best Point =====
            // Move b to u if it has smaller function magnitude

            if (Math.Abs(fu) < Math.Abs(fb))
            {
                b = u;
                fb = fu;
                dfb = dfu;
            }

            // Re-establish invariant: |fc| >= |fb|
            if (Math.Abs(fc) < Math.Abs(fb))
                Swap3(ref b, ref c, ref fb, ref fc, ref dfb, ref dfc);
        }

        return (b, Status.MaxIterations);
    }

    /// <summary>
    ///     ITP (Interpolate-Truncate-Project) root finding algorithm with smart initial guess.
    ///     Combines optimal worst-case performance of bisection with fast average-case of interpolation.
    ///     Features:
    ///     1. Smart initialization: samples interior points to find tightest bracket
    ///     2. Regula falsi interpolation for initial refinement
    ///     3. Main ITP loop with provable O(log log(1/ε)) convergence
    /// </summary>
    private static (double xsol, Status status) RootITP(
        double xmin,
        double xmax,
        Func<double, double> func)
    {
        // ===== Normalize Interval =====

        var a = xmin <= xmax ? xmin : xmax;
        var b = xmin <= xmax ? xmax : xmin;

        // ===== Input Validation =====

        if (!IsFinite(a) || !IsFinite(b))
            return (0.5 * (a + b), Status.BadInput);

        var fa = func(a);
        var fb = func(b);

        if (!IsFinite(fa) || !IsFinite(fb))
            return (0.5 * (a + b), Status.NonFinite);

        // ===== Check Endpoints =====

        if (Math.Abs(fa) <= FTOL)
            return (a, Status.OK);
        if (Math.Abs(fb) <= FTOL)
            return (b, Status.OK);

        // ===== Require Bracketing =====

        if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
            return (Math.Abs(fa) <= Math.Abs(fb) ? a : b, Status.NoBracket);

        // ===== Smart Initial Guess: Sample Interior Points =====
        // Sample N_SAMPLE_ITP points to find the narrowest bracketing interval

        Span<double> xs = stackalloc double[N_SAMPLE_ITP];
        Span<double> fs = stackalloc double[N_SAMPLE_ITP];

        xs[0] = a;
        fs[0] = fa;
        xs[N_SAMPLE_ITP - 1] = b;
        fs[N_SAMPLE_ITP - 1] = fb;

        // Sample interior points uniformly
        for (var i = 1; i < N_SAMPLE_ITP - 1; i++)
        {
            xs[i] = a + (b - a) * i / (N_SAMPLE_ITP - 1.0);
            fs[i] = func(xs[i]);

            // Early termination if we hit a root during sampling
            if (IsFinite(fs[i]) && Math.Abs(fs[i]) <= FTOL)
                return (xs[i], Status.OK);
        }

        // ===== Find Tightest Bracket Among Sampled Points =====

        double bestA = a, bestB = b, bestFa = fa, bestFb = fb;
        var bestWidth = double.PositiveInfinity;
        var foundBracket = false;

        for (var i = 0; i < N_SAMPLE_ITP - 1; i++)
        {
            if (!IsFinite(fs[i]) || !IsFinite(fs[i + 1]))
                continue;

            // Check if consecutive points bracket a root (opposite signs)
            if ((fs[i] >= 0.0 && fs[i + 1] <= 0.0) || (fs[i] <= 0.0 && fs[i + 1] >= 0.0))
            {
                var widthSeg = Math.Abs(xs[i + 1] - xs[i]);

                // Prefer narrower brackets, or if equal width, prefer smaller function values
                if (!foundBracket ||
                    widthSeg < bestWidth ||
                    (widthSeg == bestWidth &&
                     Math.Max(Math.Abs(fs[i]), Math.Abs(fs[i + 1])) <
                     Math.Max(Math.Abs(bestFa), Math.Abs(bestFb))))
                {
                    bestA = xs[i];
                    bestB = xs[i + 1];
                    bestFa = fs[i];
                    bestFb = fs[i + 1];
                    bestWidth = widthSeg;
                    foundBracket = true;
                }
            }
        }

        // Use tightest bracket found, or fall back to original interval
        var ak = foundBracket ? bestA : a;
        var bk = foundBracket ? bestB : b;
        fa = foundBracket ? bestFa : fa;
        fb = foundBracket ? bestFb : fb;

        // ===== Try Regula Falsi Interpolation =====
        // Additional refinement before main ITP loop

        if (Math.Abs(fa - fb) > EPS_MACHINE * Math.Max(Math.Abs(fa), Math.Abs(fb)))
        {
            // Linear interpolation: find x where line through (ak,fa) and (bk,fb) crosses f=0
            var xInterp = (bk * fa - ak * fb) / (fa - fb);

            if (xInterp > ak && xInterp < bk)
            {
                var fc = func(xInterp);
                if (IsFinite(fc))
                {
                    if (Math.Abs(fc) <= FTOL)
                        return (xInterp, Status.OK);

                    // Update bracket based on which side the interpolated point falls
                    if ((fa >= 0.0 && fc <= 0.0) || (fa <= 0.0 && fc >= 0.0))
                    {
                        bk = xInterp;
                        fb = fc;
                    }
                    else
                    {
                        ak = xInterp;
                        fa = fc;
                    }
                }
            }
        }

        // ===== Standard ITP Algorithm =====

        // Compute tolerance based on interval magnitude
        var mabs = Math.Max(Math.Abs(ak), Math.Abs(bk));
        var tolx = Math.Max(ATOL, RTOL * mabs);

        // Early exit if already at tolerance
        if (Math.Abs(bk - ak) <= 2.0 * tolx)
            return (0.5 * (ak + bk), Status.Tolerance);

        // Compute maximum iterations: nMax = ceil(log2((b-a)/(2*tol))) + n0
        // This guarantees convergence even in worst case
        var nHalf = (int)Math.Ceiling(Math.Log((bk - ak) / (2.0 * tolx)) / Math.Log(2.0));
        if (nHalf < 0) nHalf = 0;
        var nMax = Math.Min(nHalf + N0_ITP, MAX_ITER);

        // ===== Main ITP Iteration Loop =====

        for (var j = 0; j <= nMax; j++)
        {
            var xh = 0.5 * (ak + bk); // Bisection point
            var width = bk - ak; // Current bracket width

            // ===== Compute Projection Radius =====
            // rProj ensures convergence by limiting how far from midpoint we can go
            // Formula: r_k = 2*tol*2^(nMax-k) - width/2

            var rProj = 2.0 * tolx * Math.Pow(2.0, nMax - j) - 0.5 * width;
            if (rProj < 0.0) rProj = 0.0;

            // ===== Compute Interpolation Point =====
            // Regula falsi between current bracket endpoints

            double xf;
            if (Math.Abs(fa - fb) > EPS_MACHINE * Math.Max(Math.Abs(fa), Math.Abs(fb)))
                xf = (bk * fa - ak * fb) / (fa - fb);
            else
                xf = xh; // Fallback to bisection if denominator too small

            // ===== Compute Truncation Step =====
            // Adds controlled perturbation to avoid getting stuck
            // σ = sign(xh - xf): direction toward midpoint
            // δ = min(κ * |xh - xf|, width^2): adaptive perturbation

            double sigma = Math.Sign(xh - xf);
            var delta = width > EPS_MACHINE
                ? Math.Min(KTR_ITP * Math.Pow(width, PEXP_ITP), Math.Abs(xh - xf))
                : 0.0;

            var xt = xf + sigma * delta; // Truncated point

            // ===== Project onto Valid Region =====
            // Ensures xt stays within [xh - rProj, xh + rProj]

            double xp;
            if (xt <= xh - rProj)
                xp = xh - rProj;
            else if (xt >= xh + rProj)
                xp = xh + rProj;
            else
                xp = xt;

            // ===== Evaluate Function at Projected Point =====

            var fc = func(xp);
            if (!IsFinite(fc))
            {
                // Fallback to midpoint on non-finite value
                xp = xh;
                fc = func(xh);
                if (!IsFinite(fc))
                    return (xh, Status.NonFinite);
            }

            // ===== Check Convergence =====

            if (Math.Abs(fc) <= FTOL)
                return (xp, Status.OK);
            if (Math.Abs(width) <= 2.0 * tolx)
                return (xp, Status.Tolerance);

            // ===== Update Bracket =====
            // Maintain opposite signs at endpoints

            if ((fa >= 0.0 && fc <= 0.0) || (fa <= 0.0 && fc >= 0.0))
            {
                bk = xp;
                fb = fc;
            }
            else
            {
                ak = xp;
                fa = fc;
            }
        }

        return (Math.Abs(fa) <= Math.Abs(fb) ? ak : bk, Status.MaxIterations);
    }

    // ===== Helper Methods =====

    /// <summary>
    ///     Check if a value is finite (not NaN or infinity)
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal static bool IsFinite(double x)
    {
        return !double.IsNaN(x) && !double.IsInfinity(x);
    }

    /// <summary>
    ///     Check if two values have opposite signs (strict definition for bracketing).
    ///     Returns false if either value is zero, since zero indicates we've found a root.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool OppositeSign(double x, double y)
    {
        return (x < 0 && y > 0) || (x > 0 && y < 0);
    }

    /// <summary>
    ///     Unit roundoff relative to x (machine epsilon scaled by magnitude)
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double UnitRoundoff(double x)
    {
        return EPS_MACHINE * Math.Max(1.0, Math.Abs(x));
    }

    /// <summary>
    ///     Spacing function: distance to next representable double (ULP-aware)
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Spacing(double x)
    {
        return UnitRoundoff(x);
    }

    /// <summary>
    ///     ULP-based scale for convergence checking.
    ///     Returns max(ATOL, RTOL * max(1, |x|))
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double UlpScale(double x)
    {
        return Math.Max(ATOL, RTOL * Math.Max(1.0, Math.Abs(x)));
    }

    /// <summary>
    ///     Check if derivative is good enough to use for Newton step.
    ///     Criteria:
    ///     1. Derivative must be finite
    ///     2. Derivative must not be too small (avoid division by near-zero)
    ///     3. Newton step must be within allowed bounds
    ///     4. Function must not be too flat (|f| vs |f'| * step check)
    /// </summary>
    private static bool DerivGood(double fx, double dfx, double stepBound)
    {
        // Check 1: Derivative must be finite
        if (!IsFinite(dfx))
            return false;

        // Check 2: Derivative must be significant relative to its magnitude
        if (Math.Abs(dfx) <= UnitRoundoff(dfx))
            return false;

        // Check 3: Newton step must be within bounds
        var stepN = -fx / dfx;
        if (Math.Abs(stepN) > stepBound)
            return false;

        // Check 4: Function must not be too flat (ensures sufficient decrease)
        if (Math.Abs(fx) > Math.Abs(dfx) * stepBound * DERIV_STEP_CHECK_FACTOR)
            return false;

        return true;
    }

    /// <summary>
    ///     Evaluate function and gradient, catching exceptions and checking for finite results.
    ///     Sets derivative to zero if function is non-finite or if derivative is non-finite.
    /// </summary>
    /// <returns>Tuple of (f, df, ok) where ok indicates if function value is finite</returns>
    private static (double f, double df, bool ok) EvalFuncWithGradient(
        double x,
        Func<double, (double f, double df)> func)
    {
        try
        {
            var (f, df) = func(x);
            var ok = IsFinite(f);

            // Set derivative to zero if function is non-finite (safety)
            if (!ok) df = 0.0;

            // Set derivative to zero if it's non-finite (but keep trying if f is finite)
            if (!IsFinite(df)) df = 0.0;

            return (f, df, ok);
        }
        catch
        {
            // Return safe values on exception
            return (0.0, 0.0, false);
        }
    }

    /// <summary>
    ///     Attempt to recover from non-finite function value by trying nearby points.
    ///     Strategy:
    ///     1. Try midpoint between x0 and x1
    ///     2. Try golden ratio point (approximately 0.382 of the way from x0 to x1)
    /// </summary>
    private static (double f, double df, bool ok) TryRecoverWithGradient(
        double x0,
        double x1,
        Func<double, (double f, double df)> func)
    {
        // Strategy 1: Try midpoint
        var xm = 0.5 * (x0 + x1);
        var result = EvalFuncWithGradient(xm, func);
        if (result.ok)
            return result;

        // Strategy 2: Try golden ratio point (avoids rational fractions)
        var xg = x0 + GOLDEN_RATIO_COMPLEMENT * (x1 - x0);
        return EvalFuncWithGradient(xg, func);
    }

    /// <summary>
    ///     Ensure step size is at least the unit roundoff to avoid stalling.
    ///     Modifies 'step' in place if it's too small.
    /// </summary>
    private static void EnforceMinStep(double xfrom, double direction, ref double step)
    {
        var tmin = UnitRoundoff(xfrom);
        if (Math.Abs(step) < tmin)
            step = Math.Sign(direction) * tmin;
    }

    /// <summary>
    ///     Swap three pairs of values simultaneously (for algorithm state updates)
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Swap3(
        ref double x1, ref double x2,
        ref double f1, ref double f2,
        ref double d1, ref double d2)
    {
        (x1, x2) = (x2, x1);
        (f1, f2) = (f2, f1);
        (d1, d2) = (d2, d1);
    }
}

// ===== Trust Region Newton-Dogleg Solver =====

/// <summary>
///     Delegate for computing residuals: r = f(x) where f: R^n -> R^n
/// </summary>
public delegate void ResidualFunc(ReadOnlySpan<double> x, Span<double> r);

/// <summary>
///     Delegate for solving linear system with Jacobian: solve J(x) * p = rhs for p.
///     The Jacobian is implicitly evaluated at the current x (captured by closure).
/// </summary>
public delegate void LinearSolveFunc(ReadOnlySpan<double> rhs, Span<double> p);

/// <summary>
///     Delegate for Jacobian-vector products at current x (captured by closure).
///     If transpose==false: y = J(x) * v (forward mode)
///     If transpose==true:  y = J(x)^T * v (adjoint mode, strongly recommended for performance)
///     If your implementation doesn't support transpose, ignore the flag; solver will fall back to numerical J^T*r.
/// </summary>
public delegate void JacobianVectorFunc(ReadOnlySpan<double> v, Span<double> y, bool transpose);

/// <summary>
///     Configuration options for Trust Region Newton-Dogleg solver
/// </summary>
public readonly record struct TRNOptions(
    int MaxIterations = 50, // Maximum number of iterations
    double ErrorTolerance = 1e-8, // Convergence tolerance for scaled step size
    double AcceptanceEta = 1e-8, // Minimum trust ratio ρ for accepting step (actual/predicted reduction)
    int MaxBacktracks = 5, // Maximum number of backtracking steps along dogleg direction
    double ModelDenomFloor = 1e-30, // Floor value to prevent division by zero in trust ratio
    bool HasJTMultiply = true, // Set false if Jv doesn't support transpose; solver uses numerical J^T*r
    bool EnforcePerCoordinateBox = false // If true, additionally enforce |p_i| <= scale_i * Delta per coordinate
);

/// <summary>
///     Result returned by Trust Region Newton-Dogleg solver
/// </summary>
public readonly record struct TRNResult(
    int Iterations, // Number of iterations performed
    double FinalPhi, // Final objective value: 0.5 * ||r(x)||^2
    double LastStepScaledInf, // Scaled infinity norm of last step: max_i |p_i| / scale_i
    bool Converged, // True if converged within error tolerance
    string Message // Diagnostic message
);

/// <summary>
///     Trust Region Newton-Dogleg solver for nonlinear systems f(x) = 0.
///     Uses a fixed-radius trust region (Delta=1) with per-coordinate scaling and dogleg step strategy.
///     Convergence criterion: scaled infinity-norm of accepted step falls below ErrorTolerance.
/// </summary>
public static class TrustRegionNewtonDogleg
{
    /// <summary>
    ///     Solve nonlinear system f(x) = 0 using trust region method with dogleg step.
    ///     Trust Region: The step p is constrained to satisfy ||D^{-1} p||_2 <= Delta = 1,
    ///     where D = diag(scale) and scale[i] represents the maximum allowable displacement for x[i].
    ///     Dogleg Strategy: Combines Cauchy point (steepest descent) and Newton point:
    ///     - If Newton step fits in trust region: use it (full Newton)
    ///     - Otherwise: dogleg path from origin -> Cauchy point -> Newton point
    ///     Stopping Criterion: max_i |p_i| / scale_i <= ErrorTolerance
    ///     (Based on step size, not residual norm, for robust convergence detection)
    /// </summary>
    /// <param name="x0">Initial guess (length n)</param>
    /// <param name="scale">Per-coordinate maximum displacement (length n, all positive)</param>
    /// <param name="residual">Function to compute residual vector r(x)</param>
    /// <param name="solveLinear">Function to solve J(x) * p = rhs for p</param>
    /// <param name="jv">Function to compute J*v and optionally J^T*v</param>
    /// <param name="xOut">Output: final solution (length n, must be pre-allocated)</param>
    /// <param name="options">Solver configuration options</param>
    /// <returns>Result structure with convergence information</returns>
    public static TRNResult Solve(
        ReadOnlySpan<double> x0,
        ReadOnlySpan<double> scale,
        ResidualFunc residual,
        LinearSolveFunc solveLinear,
        JacobianVectorFunc jv,
        Span<double> xOut,
        TRNOptions options = default
    )
    {
        var n = x0.Length;

        // Handle trivial case
        if (n == 0)
            return new TRNResult(0, 0.0, 0.0, true, "Trivial problem: n=0");

        // ===== Validate Inputs =====

        if (xOut.Length != n)
            throw new ArgumentException($"xOut length ({xOut.Length}) must equal x0 length ({n}).");

        if (scale.Length != n)
            throw new ArgumentException($"scale length ({scale.Length}) must equal x0 length ({n}).");

        // Validate array contents for finiteness and positivity
        for (var i = 0; i < n; ++i)
        {
            if (!RootFinder.IsFinite(x0[i]))
                throw new ArgumentException($"x0[{i}] = {x0[i]} is not finite (NaN or Inf).");

            if (!(scale[i] > 0.0) || !RootFinder.IsFinite(scale[i]))
                throw new ArgumentException($"scale[{i}] = {scale[i]} must be finite and positive.");
        }

        // Fixed trust radius in scaled coordinates (Delta = 1)
        const double Delta = 1.0;

        // ===== Allocate Workspace from Pool =====
        // Using ArrayPool reduces GC pressure for repeated calls

        var pool = ArrayPool<double>.Shared;
        var xArr = pool.Rent(n); // Current iterate
        var rArr = pool.Rent(n); // Residual at current x
        var rTrialArr = pool.Rent(n); // Residual at trial point
        var gArr = pool.Rent(n); // Gradient g = J^T * r
        var JgArr = pool.Rent(n); // J * g (for Cauchy point)
        var pUArr = pool.Rent(n); // Cauchy (steepest descent) point
        var pBArr = pool.Rent(n); // Newton (Gauss-Newton) point
        var pArr = pool.Rent(n); // Final dogleg step
        var tmpArr = pool.Rent(n); // Temporary workspace (rhs, Jp, etc.)
        var xTrialArr = pool.Rent(n); // Trial point x + step*p
        var bArr = pool.Rent(n); // Vector pB - pU for dogleg computation

        try
        {
            // ===== Initialize Iterate and Residual =====

            x0.CopyTo(xArr);
            var x = xArr.AsSpan(0, n);
            var r = rArr.AsSpan(0, n);
            residual(x, r);

            // Objective function: φ(x) = 0.5 * ||r(x)||^2
            var phi = 0.5 * Dot(r, r);

            var lastStepScaledInf = double.PositiveInfinity;

            // ===== Main Iteration Loop =====

            for (var iter = 0; iter < options.MaxIterations; ++iter)
            {
                // ===== Compute Gradient g = J^T * r =====
                // Use user-provided transpose if available, otherwise numerical differentiation

                var g = gArr.AsSpan(0, n);
                if (options.HasJTMultiply)
                    jv(r, g, true); // Adjoint mode: g = J^T * r
                else
                    NumericalJTtimesVector(x, scale, residual, r, g);

                // ===== Compute Cauchy Point (Steepest Descent) =====
                // pU = -α * g, where α = (g·g) / (Jg·Jg) minimizes φ(x - t*g) in quadratic model

                var Jg = JgArr.AsSpan(0, n);
                jv(g, Jg, false); // Forward mode: Jg = J * g

                var gg = Dot(g, g); // ||g||^2
                var JgJg = Dot(Jg, Jg); // ||J*g||^2
                var alpha = gg / Math.Max(JgJg, options.ModelDenomFloor); // Optimal Cauchy step length

                var pU = pUArr.AsSpan(0, n);
                ScaleTo(-alpha, g, pU); // pU = -α * g

                // ===== Compute Newton Point (Gauss-Newton) =====
                // Solve J * pB = -r for pB

                var pB = pBArr.AsSpan(0, n);
                var rhs = tmpArr.AsSpan(0, n);
                ScaleTo(-1.0, r, rhs); // rhs = -r
                solveLinear(rhs, pB); // pB = J^{-1} * (-r)

                // ===== Compute Dogleg Step Inside Trust Region =====
                // Trust region: ||D^{-1} * p||_2 <= Delta, where D = diag(scale)

                var p = pArr.AsSpan(0, n);
                var nB = Scaled2Norm(pB, scale); // ||pB||_D (scaled 2-norm of Newton point)
                var nU = Scaled2Norm(pU, scale); // ||pU||_D (scaled 2-norm of Cauchy point)

                if (nB <= Delta)
                {
                    // Case 1: Full Newton step fits in trust region
                    pB.CopyTo(p);
                }
                else if (nU >= Delta)
                {
                    // Case 2: Even Cauchy point is outside trust region
                    // Take radial projection of Cauchy point onto trust region boundary
                    var s = Delta / Math.Max(nU, options.ModelDenomFloor);
                    AxpyTo(s, pU, default, p); // p = s * pU
                }
                else
                {
                    // Case 3: Cauchy inside, Newton outside
                    // Take dogleg path: p = pU + τ * (pB - pU) with τ ∈ (0,1) such that ||p||_D = Delta

                    var b = bArr.AsSpan(0, n);
                    for (var i = 0; i < n; ++i)
                        b[i] = pB[i] - pU[i];

                    // Work in scaled coordinates: A = pU/scale, B = b/scale
                    // Solve quadratic: ||A + τ*B||^2 = Delta^2
                    // Expanded: (B·B) τ^2 + 2(A·B) τ + (A·A - Delta^2) = 0

                    double a2 = 0.0, b2 = 0.0, ab = 0.0;
                    for (var i = 0; i < n; ++i)
                    {
                        var A = pU[i] / scale[i];
                        var B = b[i] / scale[i];
                        a2 += A * A; // ||A||^2
                        b2 += B * B; // ||B||^2
                        ab += A * B; // A·B
                    }

                    var Delta2 = Delta * Delta;

                    // Discriminant: Δ = (A·B)^2 + (B·B) * (Delta^2 - A·A)
                    var disc = ab * ab + b2 * (Delta2 - a2);
                    var tau = 0.0;

                    if (disc <= 0.0 || b2 <= options.ModelDenomFloor)
                    {
                        // Degenerate case: no valid intersection or b too small
                        // Fallback: take scaled projection along b direction
                        var normB = Math.Sqrt(b2);
                        tau = (Delta - Math.Sqrt(a2)) / Math.Max(normB, options.ModelDenomFloor);
                        tau = Math.Clamp(tau, 0.0, 1.0);
                    }
                    else
                    {
                        // Standard case: τ = (-A·B + √Δ) / (B·B)
                        // Take positive root to move toward pB
                        tau = (-ab + Math.Sqrt(disc)) / Math.Max(b2, options.ModelDenomFloor);
                        tau = Math.Clamp(tau, 0.0, 1.0);
                    }

                    // Construct dogleg step: p = pU + τ * (pB - pU)
                    for (var i = 0; i < n; ++i)
                        p[i] = pU[i] + tau * (pB[i] - pU[i]);
                }

                // ===== Optional Per-Coordinate Box Constraint =====
                // Additionally enforce |p_i| <= Delta * scale_i for each coordinate

                if (options.EnforcePerCoordinateBox)
                    for (var i = 0; i < n; ++i)
                    {
                        var lim = Delta * scale[i];
                        if (p[i] > lim) p[i] = lim;
                        if (p[i] < -lim) p[i] = -lim;
                    }

                // ===== Compute Predicted Reduction =====
                // For Gauss-Newton model: pred(s) = -s*(g·p) - 0.5*s^2*||Jp||^2
                // where s is the step scaling factor

                var Jp = tmpArr.AsSpan(0, n);
                jv(p, Jp, false); // Jp = J * p

                var gp = Dot(g, p); // Directional derivative
                var Jp2 = Dot(Jp, Jp); // ||J*p||^2

                // ===== Backtracking Line Search Along Dogleg Direction =====

                var xTrial = xTrialArr.AsSpan(0, n);
                var rTrial = rTrialArr.AsSpan(0, n);

                var stepScale = 1.0; // Initial step scaling (full step)
                var phiTrial = double.PositiveInfinity;
                var rho = double.NegativeInfinity; // Trust ratio: actual/predicted reduction

                for (var bt = 0; bt <= options.MaxBacktracks; ++bt)
                {
                    // Try scaled step: x_trial = x + stepScale * p
                    for (var i = 0; i < n; ++i)
                        xTrial[i] = x[i] + stepScale * p[i];

                    residual(xTrial, rTrial);
                    phiTrial = 0.5 * Dot(rTrial, rTrial);

                    // Compute trust ratio ρ = (actual reduction) / (predicted reduction)
                    var predS = -stepScale * gp - 0.5 * stepScale * stepScale * Jp2;

                    if (predS > 0.0)
                        rho = (phi - phiTrial) / Math.Max(predS, options.ModelDenomFloor);
                    else
                        // Predicted reduction is non-positive; accept if actual decrease
                        rho = phi - phiTrial > 0 ? 1.0 : double.NegativeInfinity;

                    // Accept if: (1) objective decreased, AND (2) trust ratio sufficient
                    if (phiTrial < phi && rho >= options.AcceptanceEta)
                        break;

                    // Give up after max backtracks (accept best attempt even if not good)
                    if (bt == options.MaxBacktracks)
                        break;

                    // Reduce step size and try again
                    stepScale *= RootFinder.BACKTRACK_FACTOR;
                }

                // ===== Decide Step Acceptance =====

                var accept = phiTrial < phi && rho >= options.AcceptanceEta;

                if (!accept)
                    // Reject step: keep current x, r, phi; trust region may need adjustment
                    // With fixed trust region, we just try again next iteration
                    continue;

                // ===== Accept Step =====

                for (var i = 0; i < n; ++i)
                    x[i] = xTrial[i];

                rTrial.CopyTo(r);
                phi = phiTrial;

                // ===== Check Convergence =====
                // Compute scaled infinity norm of accepted step: max_i |stepScale * p_i| / scale_i

                lastStepScaledInf = 0.0;
                for (var i = 0; i < n; ++i)
                    lastStepScaledInf = Math.Max(lastStepScaledInf, Math.Abs(stepScale * p[i]) / scale[i]);

                if (lastStepScaledInf <= options.ErrorTolerance)
                    return new TRNResult(
                        iter + 1,
                        phi,
                        lastStepScaledInf,
                        true,
                        "Converged: scaled step size below tolerance.");
            }

            // Maximum iterations reached without convergence
            return new TRNResult(
                options.MaxIterations,
                phi,
                lastStepScaledInf,
                false,
                "Maximum iterations reached without convergence.");
        }
        finally
        {
            // ===== Clean Up =====
            // CRITICAL: Copy final solution to output BEFORE returning arrays to pool

            new ReadOnlySpan<double>(xArr, 0, n).CopyTo(xOut);

            // Return all rented arrays to pool
            pool.Return(xArr);
            pool.Return(rArr);
            pool.Return(rTrialArr);
            pool.Return(gArr);
            pool.Return(JgArr);
            pool.Return(pUArr);
            pool.Return(pBArr);
            pool.Return(pArr);
            pool.Return(tmpArr);
            pool.Return(xTrialArr);
            pool.Return(bArr);
        }
    }

    // ===== Vector Operation Helpers =====

    /// <summary>
    ///     Dot product: a · b
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Dot(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        var s = 0.0;
        for (var i = 0; i < a.Length; ++i)
            s += a[i] * b[i];
        return s;
    }

    /// <summary>
    ///     Scale vector: out = α * a
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void ScaleTo(double alpha, ReadOnlySpan<double> a, Span<double> outVec)
    {
        for (var i = 0; i < a.Length; ++i)
            outVec[i] = alpha * a[i];
    }

    /// <summary>
    ///     AXPY operation: out = α * x + y.
    ///     If y is empty, performs out = α * x.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void AxpyTo(double alpha, ReadOnlySpan<double> x, ReadOnlySpan<double> y, Span<double> outVec)
    {
        if (y.Length == 0)
            for (var i = 0; i < x.Length; ++i)
                outVec[i] = alpha * x[i];
        else
            for (var i = 0; i < x.Length; ++i)
                outVec[i] = alpha * x[i] + y[i];
    }

    /// <summary>
    ///     Compute scaled 2-norm: ||v||_D where D = diag(scale).
    ///     Returns: sqrt(sum((v_i / scale_i)^2))
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Scaled2Norm(ReadOnlySpan<double> v, ReadOnlySpan<double> scale)
    {
        var s = 0.0;
        for (var i = 0; i < v.Length; ++i)
        {
            var t = v[i] / scale[i];
            s += t * t;
        }

        return Math.Sqrt(s);
    }

    /// <summary>
    ///     Fallback computation of g = J^T * r via numerical differentiation of φ(x) = 0.5 * ||r(x)||^2.
    ///     Uses centered finite differences: g_i ≈ (φ(x + h*e_i) - φ(x - h*e_i)) / (2h)
    ///     where h = 1e-8 * max(scale[i], 1.0) for adaptive step sizing.
    /// </summary>
    private static void NumericalJTtimesVector(
        ReadOnlySpan<double> x,
        ReadOnlySpan<double> scale,
        ResidualFunc residual,
        ReadOnlySpan<double> rAtX,
        Span<double> gOut)
    {
        var n = x.Length;
        var pool = ArrayPool<double>.Shared;
        var xWorkArr = pool.Rent(n);
        var rWorkArr = pool.Rent(n);

        try
        {
            x.CopyTo(xWorkArr);
            var xw = xWorkArr.AsSpan(0, n);
            var rw = rWorkArr.AsSpan(0, n);

            for (var i = 0; i < n; ++i)
            {
                // Adaptive step size based on coordinate scale
                var hi = 1e-8 * Math.Max(1.0, scale[i]);
                var xi = xw[i];

                // Forward perturbation: φ(x + h*e_i)
                xw[i] = xi + hi;
                residual(xw, rw);
                var phiPlus = 0.5 * Dot(rw, rw);

                // Backward perturbation: φ(x - h*e_i)
                xw[i] = xi - hi;
                residual(xw, rw);
                var phiMinus = 0.5 * Dot(rw, rw);

                // Restore original value
                xw[i] = xi;

                // Centered difference: g_i = (φ+ - φ-) / (2h)
                gOut[i] = (phiPlus - phiMinus) / (2.0 * hi);
            }
        }
        finally
        {
            pool.Return(xWorkArr);
            pool.Return(rWorkArr);
        }
    }
}