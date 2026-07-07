# Dead-Code Sweep (conservative)

Branch: `overnight-todo-sweep` — TODO low "Clean up unused / dead code".
Date: 2026-07-07.

## Summary

- **Truly-dead items removed (category a): 0**
- **Stale `#[allow(dead_code)]` attributes dropped (category b): 0** (none exist to drop)
- **Kept as experiment / deferred / live-but-pub (category c): all candidates**

No source files were changed. The systematic mechanisms this sweep relies on both
came back empty, and the only theoretical remaining dead code (unreferenced `pub`
items) is either an actively-iterated experiment path or reachable from an
example/bin target that the reference tooling does not index.

## What was checked

1. **`#[allow(dead_code)]` / `#[allow(unused*)]` / `#[expect(dead_code)]` attributes**
   across `rust/flashlfq-core/src` and `rust/flashlfq-py`.
   → **None found.** The only `#[allow(...)]` attributes in the tree are
   `#[allow(clippy::too_many_arguments)]` (14 sites), which are lint-noise
   suppressions, not dead-code markers, and are out of scope.

2. **Compiler dead-code / unused warnings.** Forced a full recompile of
   `flashlfq-core` (touched `lib.rs`), then built the **entire workspace with
   `--all-targets`** (lib, bins, examples, tests, benches) in release.
   → **0 warnings.** No `dead_code`, `never used`, `never read`, or
   `never constructed` diagnostics anywhere. Every private item is reachable.

3. **`cargo clippy --workspace --all-targets`.**
   → **0 dead/unused findings.**

4. **Unreferenced `pub` items** (which never trigger `dead_code`). Spot-checked the
   most recently-added public surface, the msalign exporter in
   `feature_export.rs` (`write_ms1_align_file`, `resolved_to_ms1align_features`,
   `write_ms1_align`, `Ms1AlignFeature`). Serena `find_referencing_symbols`
   returned **zero** references for these — but a repo-wide grep shows they ARE
   used by `examples/detect_features_tsv.rs`, gated behind the `MSALIGN_OUT` env
   var (commit e212d14). **Lesson: an empty Serena reference result is NOT proof a
   `pub` item is dead** — Serena does not index example/bin targets, so pub-item
   removal cannot be justified from Serena alone.

## Classification of candidates

| Candidate | Class | Reason |
|---|---|---|
| (any `#[allow(dead_code)]`) | — | None exist in the codebase |
| All private fns/structs in `flashlfq-core` | (c) live | 0 compiler/clippy dead-code warnings ⇒ all reachable |
| `feature_export::{write_ms1_align, write_ms1_align_file, resolved_to_ms1align_features, Ms1AlignFeature}` | (c) keep | msalign exporter (e212d14); live via `examples/detect_features_tsv.rs`, gated by `MSALIGN_OUT` |
| `joint_fit.rs` / `JOINT_FIT` path | (c) keep | opt-in joint multi-envelope decon experiment (6305f9c) |
| `feature_refinement::refine_feature_shift_neighbor` + `REFINE_METHOD` gate | (c) keep | opt-in score-ordered neighbor refinement (34bc15a, 6a45ee5, fdac6b9) |
| `shift_composite` / spectral-averaging path | (c) keep | averaging composite; gated param (46c1af4, 805017a) |
| `isotope_shift_decon.rs` | (c) keep | FlashLFQ-style −1/0/+1 shift-decon primitive |
| `trace_kernel.rs` experiment gates (`DETECT_PROFILE`, etc.) | (c) keep | env-gated detector experiments |

## Deliberately KEPT experiment paths

Per the sweep's safety constraints, every opt-in / env- or param-gated experiment
was preserved untouched:

- Joint multi-envelope decon — `joint_fit.rs`, `JOINT_FIT`
- Score-ordered neighbor refinement — `refine_feature_shift_neighbor`, `REFINE_METHOD`
- Spectral-averaging composite — `shift_composite`
- Isotope shift-decon primitive — `isotope_shift_decon.rs`
- msalign exporter — `feature_export.rs`, `MSALIGN_OUT`
- Detector profiling gates — `trace_kernel.rs`, `DETECT_PROFILE`

## Verification

- `cargo build --release --workspace --all-targets` → green, **0 warnings**.
- `cargo clippy --release --workspace --all-targets` → green, **0 dead/unused**.
- No source files modified, so the pre-existing green test baseline is unchanged
  (known unrelated failure `raw_quant_matches_mzml_quant` excepted).

## Outcome

Nothing was clearly safe to remove under the conservative mandate, so nothing was
removed. This is the documented, intended outcome: the tree is already free of
dead-code markers and compiler/clippy-detectable dead code, and the remaining
low-reference `pub` items are live experiment paths that must be preserved.
