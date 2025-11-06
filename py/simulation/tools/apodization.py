from enum import Enum, unique


@unique
class apo_types(Enum):
    NONE = 0
    HAMMING = 1
    HANNING = 2
    GAUSSIAN = 3
