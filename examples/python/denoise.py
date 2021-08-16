import os
import sys

import numpy as np

from tinymesh import Mesh, hole_fill, smooth_taubin, denoise_normal_gaussian


def main(filename):
    mesh = Mesh(filename)

    # Fill holes
    hole_fill(mesh, np.pi / 6.0)
    print('[ OK ] hole filling')

    # Add noise
    for i in range(mesh.num_vertices()):
        v = mesh.vertex(i)
        pos = v.pos() + np.random.uniform(-0.01, 0.01, size=(3))
        v.set_pos(pos)

    base, ext = os.path.splitext(filename)
    outfile = base + "_noise" + ext
    mesh.save(outfile)

    # Then, denoise
    denoise_normal_gaussian(mesh, sigma=0.2, iterations=10)
    print('[ OK ] denoising')

    # Save
    base, ext = os.path.splitext(filename)
    outfile = base + "_denoise" + ext
    mesh.save(outfile)
    print('[ OK ] saved to %s' % (outfile))


if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print("[USAGE] python %s %s" % (os.path.basename(__file__), "INPUT_FILE"))
    else:
        main(sys.argv[1])
