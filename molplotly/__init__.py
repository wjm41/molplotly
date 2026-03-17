from importlib.metadata import PackageNotFoundError, version

from .main import add_molecules

__all__ = ["add_molecules"]

try:
    __version__ = version("molplotly")
except PackageNotFoundError:
    __version__ = "unknown"
