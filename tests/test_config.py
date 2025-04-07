import pytest

from gbtools.utils import EntrezConfig, load_config


class TestEntrezConfig:
    def test_with_api_key(self):
        """Expect config initialization without API key to succeed
        with a default value of `None`"""
        config = EntrezConfig(email="jegsamson.dev@gmail.com", api_key="dummy-api-key")
        assert hasattr(config, "email")
        assert hasattr(config, "api_key")
        assert hasattr(config, "max_tries")
        assert hasattr(config, "sleep_between_tries")
        assert hasattr(config, "path")
        assert config.allowed_queries_per_second == 10

    def test_without_email(self):
        """Expect config initialization without user email key to fail"""
        with pytest.raises(TypeError) as excinfo:
            EntrezConfig(api_key="dummy-api-key")

    def test_without_api_key(self):
        """Initialize config without NCBI API key"""
        config = EntrezConfig(email="jegsamson.dev@gmail.com")
        assert hasattr(config, "email")
        assert hasattr(config, "api_key")
        assert hasattr(config, "max_tries")
        assert hasattr(config, "sleep_between_tries")
        assert hasattr(config, "path")
        assert config.api_key is None
        assert config.allowed_queries_per_second == 3

    def test_config_save(self):
        """Expect to save or overwrite config file located at /tmp/gbtools_auth.json"""
        config = EntrezConfig(email="jegsamson.dev@gmail.com", api_key="dummy-api-key")
        config.save()
        assert config.path.exists() == True
        assert config.path.is_file() == True


def test_load_config():
    """Expect to load JSON file as an EntrezConfig object"""
    config = load_config()
    assert hasattr(config, "email")
    assert hasattr(config, "api_key")
    assert hasattr(config, "max_tries")
    assert hasattr(config, "sleep_between_tries")
    assert hasattr(config, "path")
