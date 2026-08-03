"""Provide CLI for application."""

import logging

import click

from gene_jar import __version__
from gene_jar.config import get_config
from gene_jar.utils import initialize_logs

@click.group()
@click.version_option(__version__)
def cli() -> None:
    """Short description of CLI.

    \b
        $ echo "provide a multiline description with a leading \\b"
        $ echo "more commands here"
        $ echo "otherwise, a single indent will pick up proper formatting"

    Conclude by summarizing additional commands
    """  # noqa: D301
    log_level = logging.DEBUG if get_config().debug else logging.INFO
    initialize_logs(log_level)
