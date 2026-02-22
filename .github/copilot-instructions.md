## 0. Operating Principles
1. **Work methodically.** Proceed carefully, step-by-step, and verify assumptions as you go.
2. **Be explicit about limits.** Never claim code is “production-ready.” State uncertainties, trade-offs, and what was not verified.
3. **Context discipline.** Assume you may miss context; actively reduce ambiguity:

   * Build a *codemap* (files, call graph, config keys, dependencies) before planned changes.
   * Use external tools (search, grep, tests, linters, type-checkers) when needed.
4. **Self-check loop.** Regularly ask: *“Am I doing this correctly? What could I be missing?”*

## 1. Architecture & Code Organization

1. **Modular design.** Use clear separation of concerns and only the abstractions that pay for themselves.
2. **File granularity.** Prefer separate files for each major class / component / reusable function set, when it improves clarity and reuse.
3. **Remove dead legacy.** During updates, delete obsolete code paths/config logic unless the user explicitly asks to keep backward compatibility.

## 2. Configuration (YAML-first, strict, validated)

1. **No hardcoded parameters.** Any tunable value must be externalized to configuration (YAML by default).
2. **No implicit defaults in runtime logic.**

   * If a parameter is required and missing from YAML → **raise an error**.
   * Do **not** silently fall back to defaults inside the code.
3. **Allowed placeholders in code (schema only).**

   * Configuration declarations in code may use *only* placeholders such as: `None/null`, empty dict `{}`, empty list `[]`, or equivalent language constructs.
   * These placeholders exist solely to define structure/types and must be overridden by YAML at runtime.
4. **Early validation.**

   * Immediately after loading YAML: validate **types, allowed values/ranges, required keys, and structural constraints**.
   * Any mismatch must fail fast with a clear error message.
5. **Single source of truth.** YAML is the authoritative source for parameter values; code defines schema/constraints, not defaults.

## 3. Security

1. **No secrets in code or config.**

   * API keys, tokens, passwords, DSNs, logging keys, etc. must come from environment variables or secret managers.
2. **Prevent accidental leakage.**

   * Do not log secrets.
   * Avoid printing raw env vars or credential-bearing URLs.

## 4. Logging (state of the art, structured, readable)

1. **No `print` for logs.** Use a proper logging package/module appropriate for the language (Python `logging`, `loguru`, etc.).
2. **Configurable logging.** Logging level, format, handlers must be configurable (preferably via YAML + env overrides if needed).
3. **Useful context.** Log lines should include:

   * function name + class or module/filename (when available)
   * timestamp and level
4. **Color for readability.**

   * Colorize levels (INFO/WARN/DEBUG/ERROR).
   * Distinguish parameter names vs values using different colors where supported.

## 5. Documentation & Instructions (Markdown conventions)

1. **Markdown is mandatory** for documentation and operational instructions.
2. **Never in project root.**

   * Documentation goes under `markdown/`
   * Instructions / runbooks / procedures go under `instructions/`
3. **Instructions are versioned.**

   * Create dated subdirectories: `instructions/YYYY-MM-DD/` (date only, no time in folder name).
   * Each instruction file must include a header with **date and time**.

## 6. Terminal / Command Execution Safety

1. **Always use timeouts** for terminal commands (prevents hangs from network/TTY/connectivity issues).
2. Prefer commands that are deterministic, with bounded output when possible.

## 7. Agent Mode & Code Modifications (governance)

1. **No unapproved edits.** In agent mode, do not modify code unless the user explicitly approves.
2. Before proposing modifications:

   * provide a codemap of the impacted area
   * describe intended changes + risks
   * then request approval (or the user’s explicit “go ahead”)