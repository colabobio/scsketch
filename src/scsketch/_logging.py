from __future__ import annotations

import logging
from typing import Literal

LogLevel = Literal["debug", "info", "warning", "error", "critical"]

_LEVELS: dict[LogLevel, int] = {
    "debug": logging.DEBUG,
    "info": logging.INFO,
    "warning": logging.WARNING,
    "error": logging.ERROR,
    "critical": logging.CRITICAL,
}


def configure_logging(verbosity: LogLevel = "warning") -> logging.Logger:
    """
    Configure scsketch logging without touching the root logger.

    Mirrors the pattern used in jupyter-scatter: a module-specific logger with its own
    handler, and a `verbosity` knob to control emitted detail.
    """
    logger = logging.getLogger("scsketch")

    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("[%(name)s] %(levelname)s: %(message)s"))
        logger.addHandler(handler)

    logger.setLevel(_LEVELS[verbosity])
    logger.propagate = False
    return logger

