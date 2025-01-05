import pytest
from variant_caller.sv_detector import SVDetector

def test_sv_detector_initialization():
    detector = SVDetector(min_size=50, max_size=1000)
    assert detector.min_size == 50
    assert detector.max_size == 1000

# Add more tests...
