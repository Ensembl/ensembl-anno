"""Logging configuration for the GMB pipeline."""

from __future__ import annotations

import logging
import os
import sys
from pathlib import Path
from typing import IO


_LOG_FILE_HANDLE: IO[str] | None = None
_ORIGINAL_STDOUT = sys.stdout
_ORIGINAL_STDERR = sys.stderr


class TeeStream:
    """Write text to both the original stream and the log file."""

    def __init__(self, stream: IO[str], log_file: IO[str]):
        self._stream = stream
        self._log_file = log_file

    def write(self, data: str) -> int:
        self._stream.write(data)
        self._log_file.write(data)
        return len(data)

    def flush(self) -> None:
        self._stream.flush()
        self._log_file.flush()

    def isatty(self) -> bool:
        return self._stream.isatty()

    def fileno(self) -> int:
        return self._stream.fileno()

    @property
    def encoding(self) -> str | None:
        return self._stream.encoding


def _install_stdio_tee(log_file: str) -> None:
    """Mirror stdout/stderr into *log_file*.

    Most of the current pipeline still reports progress with ``print``.  Tee
    streams make file logging useful immediately while allowing modules to move
    to structured logging incrementally.
    """
    global _LOG_FILE_HANDLE

    Path(log_file).parent.mkdir(parents=True, exist_ok=True)
    if _LOG_FILE_HANDLE is not None and not _LOG_FILE_HANDLE.closed:
        _LOG_FILE_HANDLE.close()

    _LOG_FILE_HANDLE = open(log_file, "a", buffering=1)
    sys.stdout = TeeStream(_ORIGINAL_STDOUT, _LOG_FILE_HANDLE)  # type: ignore[assignment]
    sys.stderr = TeeStream(_ORIGINAL_STDERR, _LOG_FILE_HANDLE)  # type: ignore[assignment]


def restore_stdio() -> None:
    """Restore stdout/stderr after installing a tee stream."""
    global _LOG_FILE_HANDLE

    sys.stdout = _ORIGINAL_STDOUT
    sys.stderr = _ORIGINAL_STDERR
    if _LOG_FILE_HANDLE is not None and not _LOG_FILE_HANDLE.closed:
        _LOG_FILE_HANDLE.close()
    _LOG_FILE_HANDLE = None


def resolve_log_file(output_dir: str, log_file: str | None, no_log_file: bool = False) -> str | None:
    """Return the effective log-file path for a pipeline run."""
    if no_log_file:
        return None
    if log_file:
        return log_file
    return os.path.join(output_dir, "gmb.log")


def setup_logging(
    level: int = logging.INFO,
    log_file: str | None = None,
    capture_stdio: bool = False,
) -> logging.Logger:
    """Configure and return the gmb root logger.

    When ``capture_stdio`` is true, stdout/stderr are mirrored to ``log_file``.
    The logger itself writes to stderr, so logging records are captured by the
    same tee without needing a separate file handler.
    """
    if log_file and capture_stdio:
        _install_stdio_tee(log_file)

    logger = logging.getLogger("gmb")
    logger.handlers.clear()
    logger.propagate = False

    formatter = logging.Formatter(
        "%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    stream_handler = logging.StreamHandler(sys.stderr)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    if log_file and not capture_stdio:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    logger.setLevel(level)
    return logger
