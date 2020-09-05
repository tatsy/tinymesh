import os
import sys
from plyfile import PlyData, PlyElement
import numpy as np
from tinymesh import Mesh


def main(filename):
    # Load .PLY with "plyfile"
    plydata = PlyData.read(filename)
    vx = plydata['vertex']['x']
    vy = plydata['vertex']['y']
    vz = plydata['vertex']['z']
    verts = np.stack([vx, vy, vz], axis=1).astype('double')
    faces = np.concatenate(plydata['face']['vertex_indices']).astype('uint32')

    # Construct halfedge structure by "tinymesh"
    mesh = Mesh(verts, faces)

    # Get vertices and indices back to Python
    verts = mesh.get_vertices()
    faces = mesh.get_vertex_indices()

    verts = np.asarray(verts, dtype='float32')
    faces = np.asarray(faces, dtype='int32').reshape((-1, 3))

    vert_el = np.array([tuple(v) for v in verts], dtype=[('x', 'f4'), ('y', 'f4'), ('z', 'f4')])
    face_el = np.array([(tuple(f),) for f in faces], dtype=[('vertex_indices', 'i4', (3,))])
    plydata = PlyData([PlyElement.describe(vert_el, 'vertex'),
                       PlyElement.describe(face_el, 'face')], text=False)

    base, ext = os.path.splitext(filename)
    outfile = base + "_copy" + ext
    plydata.write(outfile)


if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print("[USAGE] python %s %s" % (os.path.basename(__file__), "INPUT_FILE"))
    else:
        main(sys.argv[1])
