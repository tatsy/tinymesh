try:
    import unittest2 as unittest
except:
    import unittest

import numpy as np
from tinymesh import Mesh

class TestConstruct(unittest.TestCase):
    def test_tetrahedron(self):
        # Construct a simple tetrahedron
        a = np.asarray([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype='double')
        b = np.asarray([0, 1, 2, 0, 2, 3, 0, 3, 1, 3, 2, 1], dtype='uint32')
        try:
            mesh = Mesh(a, b)
        except Exception as e:
            self.fail('Mesh construction from vertices/indices was failed!')
