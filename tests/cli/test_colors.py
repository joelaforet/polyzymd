"""Tests for progress bar rendering utilities."""

from polyzymd.cli import colors


class TestRenderProgressBar:
    """Unit tests for render_progress_bar()."""

    def setup_method(self):
        """Disable color so we can inspect raw characters."""
        colors.set_color_support(colors.TerminalColorSupport.NONE)

    def test_default_width_is_40(self):
        bar = colors.render_progress_bar(0.5, "running")
        assert len(bar) == 40

    def test_custom_width(self):
        bar = colors.render_progress_bar(0.5, "running", width=20)
        assert len(bar) == 20

    def test_zero_fraction_all_empty(self):
        bar = colors.render_progress_bar(0.0, "not_started")
        assert bar == "\u2591" * 40

    def test_full_fraction_all_filled(self):
        bar = colors.render_progress_bar(1.0, "completed")
        assert bar == "\u2588" * 40

    def test_half_fraction(self):
        bar = colors.render_progress_bar(0.5, "running", width=10)
        assert bar == "\u2588" * 5 + "\u2591" * 5

    def test_fraction_clamped_above_one(self):
        bar = colors.render_progress_bar(1.5, "completed", width=10)
        assert bar == "\u2588" * 10

    def test_fraction_clamped_below_zero(self):
        bar = colors.render_progress_bar(-0.5, "not_started", width=10)
        assert bar == "\u2591" * 10

    def test_truecolor_wraps_ansi(self):
        colors.set_color_support(colors.TerminalColorSupport.TRUECOLOR)
        bar = colors.render_progress_bar(1.0, "completed", width=5)
        assert "\033[38;2;" in bar
        assert "\033[0m" in bar
        assert "\u2588" * 5 in bar

    def test_basic_color_wraps_ansi(self):
        colors.set_color_support(colors.TerminalColorSupport.BASIC)
        bar = colors.render_progress_bar(0.5, "failed", width=10)
        assert "\033[91m" in bar
        assert "\033[0m" in bar

    def test_extended_color_wraps_ansi(self):
        colors.set_color_support(colors.TerminalColorSupport.EXTENDED)
        bar = colors.render_progress_bar(0.5, "interrupted", width=10)
        assert "\033[38;5;" in bar
        assert "\033[0m" in bar

    def test_unknown_status_uses_not_found_color(self):
        colors.set_color_support(colors.TerminalColorSupport.TRUECOLOR)
        bar = colors.render_progress_bar(0.0, "some_unknown_status", width=5)
        assert "\033[38;2;128;128;128m" in bar
