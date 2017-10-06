from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


import shi7.shi7 as shi7
import shi7.shi7_learning as shi7_learning

__all__ = ["shi7", "shi7_learning"]