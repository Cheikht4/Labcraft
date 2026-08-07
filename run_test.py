import pytest
from tests.test_lamp.test_complex_enumeration import test_complex_enumeration_counts
try:
    test_complex_enumeration_counts()
except Exception as e:
    import traceback
    traceback.print_exc()
