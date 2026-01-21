from enum import Enum


class ServiceEnvironment(str, Enum):
    DEV = "dev"
    PROD = "prod"
    TEST = "test"