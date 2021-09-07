import os
import sys

import numpy as np

from tinymesh import Mesh, hole_fill, implicit_fairing, remesh_triangular


def main(filename):
    mesh = Mesh(filename)

    # Freeze faces
    for i in range(mesh.num_vertices()):
        v = mesh.vertex(i)
        v.set_is_static(True)

    # Fill holes
    hole_fill(mesh, np.pi / 6.0)
    print('[ OK ] hole filling')

    # Then, remesh
    remesh_triangular(mesh, short_length=0.5, long_length=2.0)
    print('[ OK ] remeshing')

    # Finally, fair
    implicit_fairing(mesh)
    print('[ OK ] fairing')

    # Save
    base, ext = os.path.splitext(filename)
    outfile = base + "_remesh" + ext
    mesh.save(outfile)
    print('[ OK ] saved to %s' % (outfile))


if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print("[USAGE] python %s %s" % (os.path.basename(__file__), "INPUT_FILE"))
    else:
        main(sys.argv[1])
