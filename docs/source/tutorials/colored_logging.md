# Colored Logging

PolyzyMD uses per-module colored log output so you can quickly identify
which subsystem produced each message. Colors are applied to both
structured log lines (via Python's `logging` module) and CLI status
messages (via `colored_echo`).

## Color Scheme

Each module group is assigned a distinct color for **INFO** and **DEBUG**
messages. **WARNING** is always amber/yellow and **ERROR** is always red,
regardless of module.

| Module Group | Logger prefix | Color description |
|---|---|---|
| CLI / Workflow | `polyzymd.cli`, `polyzymd.workflow` | Lavender (bright blue on 16-color terminals) |
| Building | `polyzymd.builders` | Sage green (bright green) |
| Simulation | `polyzymd.simulation.runner`, `.continuation` | Warm peach (bright magenta) |
| Progress | `polyzymd.simulation.progress`, `.signals` | Steel blue (bright cyan) |
| Core | `polyzymd.core` | Parchment (dark yellow) |
| Data | `polyzymd.data` | Aqua/mint (dark cyan) |
| Exporters | `polyzymd.exporters` | Lilac (dark magenta) |
| Utils | `polyzymd.utils` | Neutral gray (white) |

The parenthesized colors are the 16-color ANSI fallbacks used on
terminals that do not support 256-color or truecolor mode (common on HPC
login nodes).

## Terminal Support Levels

PolyzyMD auto-detects your terminal's color capability:

| Level | Detection | What you see |
|---|---|---|
| **Truecolor** (24-bit RGB) | `COLORTERM=truecolor` or `COLORTERM=24bit` | Full RGB tints from the table above |
| **256-color** | `TERM` contains `256color` | Closest xterm-256 palette match |
| **16-color** (BASIC) | Any other TTY (e.g. `TERM=xterm-16color`) | Bright/dark ANSI color fallbacks |
| **None** | Pipe, `TERM=dumb`, `NO_COLOR` set, or `--no-color` flag | No color codes emitted |

Most personal terminals (iTerm2, Windows Terminal, GNOME Terminal) support
truecolor. HPC login nodes typically report `TERM=xterm-16color` and use
the 16-color fallbacks.

## Disabling Colors

There are three ways to disable colored output:

### 1. `--no-color` CLI flag

Pass `--no-color` **before** the subcommand name:

```bash
polyzymd --no-color check-progress -c config.yaml
polyzymd --no-color build -c config.yaml
polyzymd --no-color recover -c config.yaml -r 1
```

> **Important:** `--no-color` is a *global* option and must appear before
> the subcommand. Placing it after (e.g.
> `polyzymd check-progress --no-color`) will produce an error.

### 2. `NO_COLOR` environment variable

Set `NO_COLOR` to any non-empty value. This follows the
[no-color.org](https://no-color.org) convention and is useful for
shell scripts or CI pipelines:

```bash
export NO_COLOR=1
polyzymd check-progress -c config.yaml   # no colors
```

Or for a single command:

```bash
NO_COLOR=1 polyzymd check-progress -c config.yaml
```

### 3. Automatic detection

Colors are automatically disabled when:

- Standard error is not a TTY (e.g. output is piped or redirected)
- `TERM=dumb`

## Preserved Click Styling

A small number of CLI messages use Click's built-in styling
(`click.style`) for **green** (success) and **red** (error) indicators.
These are not affected by the per-module color scheme and remain as-is
regardless of the `--no-color` flag. Examples:

- "Build complete!" (green)
- "Simulation is complete -- nothing to recover." (green)
- Error exit messages (red)

## For Developers

The color infrastructure lives in `src/polyzymd/cli/colors.py`:

- **`ColoredFormatter`** -- a `logging.Formatter` subclass that wraps log
  output in ANSI escape codes based on logger name and level.
- **`colored_echo(message, phase=..., level=...)`** -- drop-in replacement
  for `click.echo()` that applies phase-aware coloring. The `phase`
  argument maps to a module prefix (e.g. `phase="build"` uses the
  builders color).
- **`setup_colored_logging(verbose=..., no_color=...)`** -- one-call
  initialization that replaces the root logger handler with a
  `ColoredFormatter`.
- **`ColorScheme`** -- registry with longest-prefix matching. You can call
  `scheme.register(prefix, rgb, xterm256, basic_ansi)` to add custom
  module colors.
