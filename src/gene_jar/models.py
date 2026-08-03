from enum import Enum


class ServiceEnvironment(StrEnum):
    DEV = "dev"
    PROD = "prod"
    TEST = "test"