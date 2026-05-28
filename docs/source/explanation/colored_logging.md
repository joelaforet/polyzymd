# Why PolyzyMD uses colored logging

PolyzyMD workflows can produce output from several subsystems in one run:
configuration checks, system building, simulation setup, progress reporting,
analysis, and export steps. Colored logging gives those messages a visual
structure so users and new contributors can quickly tell what kind of work is
happening and how urgently they need to respond.

The colors are meant as a readability aid, not as the only source of meaning.
The message text and severity level still carry the essential information, so
logs remain useful when color is disabled in scripts, continuous integration,
or redirected output.

## Color Scheme

For routine **INFO** and **DEBUG** messages, color mainly answers "which part of
PolyzyMD is speaking?" **WARNING** is always amber/yellow and **ERROR** is
always red, regardless of subsystem, because severity should be visible before
module identity.

| Output area | What it usually means | Color description |
|---|---|---|
| CLI / workflow | Commands, orchestration, and workflow-level decisions | Lavender (bright blue on 16-color terminals) |
| Building | System preparation and molecular construction | Sage green (bright green) |
| Simulation | OpenMM run setup, continuation, and recovery messages | Warm peach (bright magenta) |
| Progress | Runtime progress and signal-handling updates | Steel blue (bright cyan) |
| Core | Shared PolyzyMD data structures and core services | Parchment (dark yellow) |
| Data | Bundled data and resource-loading messages | Aqua/mint (dark cyan) |
| Exporters | Format conversion and export messages | Lilac (dark magenta) |
| Utilities | General supporting utilities | Neutral gray (white) |

The parenthesized colors are the 16-color ANSI fallbacks used on
terminals that do not support 256-color or truecolor mode (common on HPC
login nodes).

## Terminal Support Levels

PolyzyMD adapts its own colorized log output to the terminal it detects:

| Level | Detection | What you see |
|---|---|---|
| **Truecolor** (24-bit RGB) | `COLORTERM=truecolor` or `COLORTERM=24bit` | Full RGB tints from the table above |
| **256-color** | `TERM` contains `256color` | Closest xterm-256 palette match |
| **16-color** (BASIC) | Any other TTY (e.g. `TERM=xterm-16color`) | Bright/dark ANSI color fallbacks |
| **None** | Pipe, `TERM=dumb`, `NO_COLOR` set, or `--no-color` flag | No PolyzyMD color codes emitted |

Most personal terminals (iTerm2, Windows Terminal, GNOME Terminal) support
truecolor. HPC login nodes typically report `TERM=xterm-16color` and use
the 16-color fallbacks.

## Disabling Colors

There are three ways to disable PolyzyMD's color infrastructure:

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

The flag disables PolyzyMD-managed color handling. Some command-line messages
may still use Click's own styling behavior; those messages are not governed by
PolyzyMD's subsystem color scheme.

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

## Click-styled messages

A small number of CLI messages use Click's built-in success or error styling
instead of the subsystem colors described above. Treat those messages as
plain success/error indicators rather than as signs that a particular PolyzyMD
module produced the output. Their behavior follows Click's terminal handling,
not PolyzyMD's color scheme.

Examples include:

- "Build complete!" (green)
- "Simulation is complete -- nothing to recover." (green)
- Error exit messages (red)

## Design rationale for contributors

Colored output should make PolyzyMD easier to scan without making logs harder
to capture, search, or compare. Contributors should therefore use the shared
logging and color helpers rather than adding ad hoc ANSI escape codes or
module-specific terminal checks.

Using the shared helpers keeps behavior consistent across interactive shells,
HPC login nodes, redirected output, and automation. It also ensures output
respects common controls such as `NO_COLOR`, pipes, `TERM=dumb`, and terminal
color capability. If a new workflow area needs distinct presentation, prefer
extending the shared conventions over creating a separate color system.
