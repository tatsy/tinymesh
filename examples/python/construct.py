import os
import sys
from plyfile import PlyData, PlyElement
import numpy as np
from pytinymesh import Mesh

def main(filename):
    # Construct a simple tetrahedron
    a = np.asarray([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype='double')
    b = np.asarray([0, 1, 2, 0, 2, 3, 0, 3, 1, 3, 2, 1], dtype='uint32')
    mesh = Mesh(a, b)
    mesh.save('tetra.ply')

    # Load .PLY with "plyfile" and construct a mesh by "tinymesh"
    plydata = PlyData.read(filename)
    vx = plydata['vertex']['x']
    vy = plydata['vertex']['y']
    vz = plydata['vertex']['z']
    verts = np.stack([vx, vy, vz], axis=1).astype('double')
    elems = np.concatenate(plydata['face']['vertex_indices']).astype('uint32')
    mesh = Mesh(verts, elems)
    mesh.save('output.ply')


if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print("[USAGE] python %s %s" % (os.path.basename(__file__), "INPUT_FILE"))
    else:
        main(sys.argv[1])
