"""Models used by gene-jar."""

from enum import StrEnum


class ServiceEnvironment(StrEnum):
    """Supported service environments."""

    DEV = "dev"
    PROD = "prod"
    TEST = "test"
