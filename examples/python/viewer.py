import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection


def normalize(x):
    return x / np.linalg.norm(x)


def frustum(left, right, bottom, top, znear, zfar):
    M = np.zeros((4, 4))
    M[0, 0] = +2.0 * znear / (right - left)
    M[1, 1] = +2.0 * znear / (top - bottom)
    M[2, 2] = -(zfar + znear) / (zfar - znear)
    M[0, 2] = (right + left) / (right - left)
    M[2, 1] = (top + bottom) / (top - bottom)
    M[2, 3] = -2.0 * znear * zfar / (zfar - znear)
    M[3, 2] = -1.0
    return M


class Viewer(object):
    def __init__(self, width=500, height=500, dpi=100):
        self.width = width
        self.height = height
        self.aspect = width / height
        self.dpi = 100
        self.set_identity()

    def set_identity(self):
        self.M = np.eye(4)
        self.V = np.eye(4)
        self.P = np.eye(4)

    def perspective(self, fovy, aspect, znear, zfar):
        h = np.tan(0.5 * np.radians(fovy)) * znear
        w = h * aspect
        P = frustum(-w, w, -h, h, znear, zfar)
        self.P = self.P @ P

    def lookat(self, ex, ey, ez, cx, cy, cz, ux, uy, uz):
        e = np.array((ex, ey, ez))
        c = np.array((cx, cy, cz))
        u = np.array((ux, uy, uz))
        z = normalize(e - c)
        x = normalize(np.cross(u, z))
        y = normalize(np.cross(z, x))
        rot = np.eye(4)
        rot[0, :3] = x
        rot[1, :3] = y
        rot[2, :3] = z
        tran = np.eye(4)
        tran[:3, 3] = -e
        V = rot @ tran
        self.V = self.V @ V

    def translate(self, x, y, z):
        M = np.array([
            [1.0, 0.0, 0.0, x],
            [0.0, 1.0, 0.0, y],
            [0.0, 0.0, 1.0, z],
            [0.0, 0.0, 0.0, 1.0],
        ])
        self.M = self.M @ M

    def rotate(self, theta, wx, wy, wz):
        rad = np.radians(theta)
        st = np.sin(rad)
        ct = np.cos(rad)
        K = np.array([[0.0, -wz, wy], [wz, 0.0, -wx], [-wy, wx, 0.0]])
        R = np.eye(3) + st * K + (1.0 - ct) * (K @ K)
        M = np.eye(4)
        M[:3, :3] = R
        self.M = self.M @ M

    def scale(self, sx, sy, sz):
        M = np.eye(4)
        M[0, 0] = sx
        M[1, 1] = sy
        M[2, 2] = sz
        self.M = self.M @ M

    def visualization(self,
                      verts,
                      faces,
                      colors=None,
                      wireframe=False,
                      shade=True,
                      title="",
                      save=False,
                      filename="output.png"):
        V = verts.copy()
        F = faces.copy()
        if colors is not None:
            C = colors[:, :3].copy()
        else:
            C = 0.7 * np.ones_like(verts)

        ## coordinate transformation
        MV = self.V @ self.M
        MVP = self.P @ MV
        V = V[F]
        C = C[F]

        N = np.cross(V[:, 1, :] - V[:, 0, :], V[:, 2, :] - V[:, 0, :])
        N /= np.linalg.norm(N, axis=-1, keepdims=True)
        N = N @ MV[:3, :3].T
        N /= np.linalg.norm(N, axis=-1, keepdims=True)
        L = np.c_[np.zeros(N.shape), np.ones((*N.shape[:-1], 1))] @ MV.T
        L = -L[:, :3] / L[:, 3:4]
        L /= np.linalg.norm(L, axis=-1, keepdims=True)
        NdotL = np.maximum(0.0, np.sum(N * L, axis=-1, keepdims=True))

        V = np.c_[V, np.ones((*V.shape[:-1], 1))] @ MVP.T
        V /= V[:, :, 3:4]
        T = V[:, :, :2]

        ## z-sort
        Z = -V[:, :, 2].mean(axis=1)
        zmin, zmax = Z.min(), Z.max()
        Z = (Z - zmin) / (zmax - zmin)
        l = np.argsort(Z)

        C = C.mean(axis=1)
        if shade:
            C *= NdotL
        T, C = T[l, :], C[l, :]

        if wireframe:
            linewidth = 0.1
            edgecolor = "black"
        else:
            linewidth = 1.0
            edgecolor = C

        xsiz = int(self.width / self.dpi)
        ysiz = int(self.height / self.dpi)
        fig = plt.figure(figsize=(xsiz, ysiz), dpi=self.dpi)
        ax = fig.add_axes([0, 0, 1, 1], xlim=[-1, 1], ylim=[-1, 1], aspect=1.0 / self.aspect, frameon=False)
        ax.axis("off")
        collection = PolyCollection(T,
                                    closed=True,
                                    linewidth=linewidth,
                                    facecolor=C,
                                    edgecolor=edgecolor,
                                    antialiaseds=True)
        ax.add_collection(collection)

        if len(title) != 0:
            plt.title(title)

        if save:
            plt.savefig(filename, bbox_inches="tight")

        plt.show()

    def mesh_visualization(self, mesh, *args, **kwargs):
        verts = mesh.get_vertices()
        faces = mesh.get_vertex_indices()
        verts = np.asarray(verts, dtype="float32")
        faces = np.asarray(faces, dtype="uint32").reshape((-1, 3))
        self.visualization(verts, faces, *args, **kwargs)
