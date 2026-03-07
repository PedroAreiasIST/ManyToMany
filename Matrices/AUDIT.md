# Matrices Directory Audit

Date: 2026-03-07
Scope: `Matrices/` (source, project configuration, and repository hygiene)

## Executive summary

I found three high-priority issues and one medium-priority issue:

1. **Repository hygiene issue:** build artifacts (`bin/`, `obj/`, IDE user files) are tracked in git.
2. **Potentially dangerous runtime side effect:** native dependency auto-install is enabled by default.
3. **Diagnostics gap:** many `catch` blocks swallow exceptions without logging/telemetry.
4. **API contract caveat:** `Matrix.Equals` uses tolerance while `GetHashCode` does not.

## Findings

### 1) Build artifacts are committed (High)

**What I observed**
- The `Matrices` project currently tracks generated outputs under `Matrices/bin/` and `Matrices/obj/`, plus a user-specific settings file.

**Evidence**
- `git ls-files Matrices | rg '^(Matrices/(bin|obj)/|Matrices/Folder\.DotSettings\.user)'`

**Impact**
- Frequent merge conflicts and noisy diffs.
- Large repository size and slower clone/fetch.
- Higher risk of stale/binary drift being mistaken for source changes.

**Recommendation**
- Add/verify `.gitignore` rules for `Matrices/bin/`, `Matrices/obj/`, and `*.DotSettings.user`.
- Remove already tracked generated artifacts from git history/index.

---

### 2) Auto-install of MKL is enabled by default (High)

**What I observed**
- `NativeLibraryConfig` enables package-manager-based MKL installation by default (`EnableAutoInstall = true`) and interactive prompts by default.

**Evidence**
- `EnableAutoInstall` default is `true`.
- `GetMKLLibraries()` may call `TryAutoInstallMKL()` when discovery fails/verification fails.
- Installation paths include package manager commands (`choco`, `apt-get`, `dnf`, `yum`, `brew`) and may rely on `sudo`.

**Impact**
- Surprising side effects for consumers of a numerical library.
- Security and compliance concerns in production/CI/air-gapped environments.
- Potential hangs/failures in non-interactive execution contexts.

**Recommendation**
- Flip default to opt-in (`EnableAutoInstall = false`).
- Require explicit host-application consent/configuration before any install attempt.
- Treat auto-install as a separate helper tool, not default library behavior.

---

### 3) Exception swallowing in native library discovery/install paths (High)

**What I observed**
- Numerous bare `catch` blocks suppress exceptions silently, including in CUDA version checks and path discovery.

**Evidence**
- `GetCudaVersion()` swallows exceptions while probing library loading.
- Search path building has multiple empty `catch` blocks.
- Similar patterns appear in MKL verification loops.

**Impact**
- Hard to diagnose deployment/runtime failures.
- Operators receive false negatives (“not available”) with no root-cause context.
- Increases MTTR and support burden.

**Recommendation**
- At minimum, log suppressed exceptions at debug level with context.
- Aggregate probe errors and expose them via a diagnostics API for troubleshooting.
- Reserve fully silent failure only for intentionally best-effort probes with compensating diagnostics elsewhere.

---

### 4) Equality/hash semantics mismatch is documented but risky (Medium)

**What I observed**
- `Matrix.Equals(Matrix?)` is tolerance-based, while `GetHashCode()` hashes raw sampled values.
- The code comments explicitly warn against hash-based key usage under tolerance semantics.

**Impact**
- Easy misuse by consumers (dictionary/set keys), causing hard-to-reproduce bugs.

**Recommendation**
- Consider one or more of:
  - Provide a dedicated `IEqualityComparer<Matrix>` for tolerance semantics.
  - Make default `Equals` exact and expose tolerance-based comparison as explicit methods.
  - Keep current behavior but add analyzer/docs/examples to prevent key misuse.

## Validation notes

- Attempted to run project build validation, but `dotnet` is not installed in this environment.
- Audit findings are based on source/config review and repository state inspection.
