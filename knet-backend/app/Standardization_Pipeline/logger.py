"""This module contains the logger configuration for the application."""

import logging
from datetime import datetime

from pytz import timezone


def get_logger(name: str) -> logging.Logger:
    """Get a logger with the specified name."""
    logging.Formatter.converter = lambda *args: datetime.now(tz=timezone("US/Eastern")).timetuple()
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)

    formatter = logging.Formatter("%(asctime)s - %(name)s - %(funcName)s() - %(levelname)s - %(message)s")

    stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    return logger
