"""Read and provide runtime configuration."""

from functools import cache

from pydantic_settings import BaseSettings, SettingsConfigDict

from gene_harmony.models import ServiceEnvironment


class Settings(BaseSettings):
    """Create app settings

    This is not a singleton, so every new call to this class will re-compute
    configuration settings, defaults, etc.
    """

    model_config = SettingsConfigDict(
        env_prefix="GENE_HARMONY_",
        env_file=".env",
        env_file_encoding="utf-8",
        extra="ignore"
    )

    env: ServiceEnvironment = ServiceEnvironment.DEV
    debug: bool = False
    test: bool = False


@cache
def get_config() -> Settings:
    """Get runtime configuration.

    This function is cached, so the config object only gets created/calculated once.

    :return: Settings instance
    """
    return Settings()
