import os
import sys
from tinymesh import Mesh, hole_fill, remesh_triangular


def main(filename):
    mesh = Mesh(filename)

    # Freeze faces
    for i in range(mesh.num_faces()):
        f = mesh.face(i)
        f.set_is_static(True)

    # Fill holes
    hole_fill(mesh)

    # Then, remesh
    remesh_triangular(mesh)

    # Save
    base, ext = os.path.splitext(filename)
    outfile = base + "_remesh" + ext
    mesh.save(outfile)


if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print("[USAGE] python %s %s" % (os.path.basename(__file__), "INPUT_FILE"))
    else:
        main(sys.argv[1])
