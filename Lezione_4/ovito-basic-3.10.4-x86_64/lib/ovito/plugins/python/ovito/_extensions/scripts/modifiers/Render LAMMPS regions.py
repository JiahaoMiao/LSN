##### Render LAMMPS regions #####
#
# Generates a mesh representation from LAMMPS *region* commands.
# [See documentation...](manual:modifiers.render_lammps_regions)
#
# Supported styles: *block*, *cone*, *cylinder*, *ellipsoid*, *plane*, *prism*, *sphere*.
#
# Supported keywords: *rotate*, *move*, *side*, *open*.

import ovito
from ovito.data import *
import numpy as np
from ovito.vis import SurfaceMeshVis
from ovito.pipeline import ModifierInterface
from ovito.modifiers import ColorByTypeModifier
from traits.api import Bool, Code, Instance
from ovito.traits import *
import numpy.typing as npt
from typing import Any, List
from collections import namedtuple
from abc import abstractmethod

RegionRotation = namedtuple("rotation", ["theta", "px", "py", "pz", "rx", "ry", "rz"])

class RegionStyle:
    def __init__(self):
        self.move = None
        self.rotation = None
        self.open = []
        self.side = None
        self.units = None
        self.region_id = None
        self.style = None

    def parse_general_parameters(self, tokens: List[str]):
        if tokens[0] != "region":
            raise SyntaxError("LAMMPS command does not start with expected keyword 'region'.")
        if len(tokens) <= 2:
            raise IndexError("Incomplete region command: region <ID> <style> [style-parameters] [keyword] [keyword-parameters] ...")
        # Second and third arg must be region id and region style
        self.region_id = tokens[1]
        self.style = tokens[2]
        return tokens[3:]

    def parse_kw_open(self, tokens: List[str]):
        try:
            face = int(tokens[1])
            self.open.append(face)
        except ValueError as exc:
            raise SyntaxError("Expected face index after keyword 'open'") from exc
        return tokens[2:]

    def parse_kw_units(self, tokens: List[str]):
        if tokens[1] not in ["box", "lattice"]:
            raise SyntaxError("Illegal units command")
        self.units = tokens[1]
        print('Sorry, units keyword is not supported by this modifier. Always using units "box".')
        return tokens[2:]

    def parse_kw_side(self, tokens: List[str]):
        if tokens[1] not in ["in", "out"]:
            raise SyntaxError("Expected side keyword args {in|out}")
        if tokens[1] == "out":
            self.side = 1
        return tokens[2:]

    def parse_kw_move(self, tokens: List[str]):
        try:
            self.move = list(map(float,tokens[1:4]))
        except ValueError as exc:
            raise SyntaxError("Expected keyword args <x> <y> <z> after keyword 'move'.\nNote that LAMMPS equal-style variables v_ are not supported.") from exc
        return tokens[4:]

    def parse_kw_rotate(self, tokens: List[str]):
        try:
            self.rotation = RegionRotation(*map(float,tokens[1:8]))
        except Exception as exc:
            raise SyntaxError("Expected keyword args <theta>, <px>, <py>, <pz>, <rx>, <ry>, <rz> after keyword 'rotate'.\nNote that LAMMPS equal-style variables v_ are not supported.") from exc
        return tokens[8:]

    def parse_keyword_parameters(self, keywords: List[str]):
        while keywords:
            if keywords[0] == 'open':
                keywords = self.parse_kw_open(keywords)
            elif keywords[0] == 'side':
                keywords = self.parse_kw_side(keywords)
            elif keywords[0] == 'units':
                keywords = self.parse_kw_units(keywords)
            elif keywords[0] == 'move':
                keywords = self.parse_kw_move(keywords)
            elif keywords[0] == 'rotate':
                keywords = self.parse_kw_rotate(keywords)
            elif keywords[0] == 'region':
                raise SyntaxError("Invalid keyword 'region'. Please put each LAMMPS region command on a new line.")
            else:
                raise SyntaxError("Invalid keyword argument:\n" + " ".join(keywords) + "\n^")

        return keywords

    def rotate_verts(self, verts: npt.ArrayLike):
        u = np.array((self.rotation.rx, self.rotation.ry, self.rotation.rz))
        u_x, u_y, u_z = u/np.linalg.norm(u)
        sin, cos = np.sin(self.rotation.theta), np.cos(self.rotation.theta)

        R = np.array(((cos + u_x**2*(1-cos), u_y*u_x*(1-cos)+u_z*sin, u_z*u_x*(1-cos) - u_y*sin),
                    (u_x*u_y*(1-cos) - u_z*sin, cos + u_y**2*(1-cos), u_z*u_y*(1-cos)+u_x*sin),
                    (u_x*u_z*(1-cos) + u_y*sin, u_y*u_z*(1-cos) - u_x*sin, cos + u_z**2*(1-cos))))

        offset = np.array((self.rotation.px, self.rotation.py, self.rotation.pz))
        verts -= offset
        verts = (R.T @ verts.T).T
        verts += offset
        return verts

    def inv_rotate_verts(self, verts: npt.ArrayLike):
        u = np.array((-self.rotation.rx, -self.rotation.ry, -self.rotation.rz))
        u_x, u_y, u_z = u/np.linalg.norm(u)
        sin, cos = np.sin(self.rotation.theta), np.cos(self.rotation.theta)

        R = np.array(((cos + u_x**2*(1-cos), u_y*u_x*(1-cos)+u_z*sin, u_z*u_x*(1-cos) - u_y*sin),
                    (u_x*u_y*(1-cos) - u_z*sin, cos + u_y**2*(1-cos), u_z*u_y*(1-cos)+u_x*sin),
                    (u_x*u_z*(1-cos) + u_y*sin, u_y*u_z*(1-cos) - u_x*sin, cos + u_z**2*(1-cos))))

        offset = np.array((self.rotation.px, self.rotation.py, self.rotation.pz))
        verts = verts - offset
        verts = (R.T @ verts.T).T
        verts = verts + offset
        return verts

    def generate_mesh(self, mesh: SurfaceMesh):
        # Add element type to 'Region ID' property of SurfaceMesh
        mesh.regions_.count += 1
        mesh.regions_['Region ID_'][-1] = mesh.regions.count
        # LAMMPS region IDs should be unique
        existing_region_id = mesh.regions['Region ID'].type_by_name(self.region_id, raise_error=False)
        if existing_region_id:
            print(f"Warning: Region ID '{self.region_id}' already exists in the output mesh (numeric id {existing_region_id.id}). Consider renaming it.")
        region_type = mesh.regions_['Region ID_'].add_type_id(mesh.regions.count, mesh.regions, self.region_id)
        # Generate vertices and faces
        verts, faces = self.generate_verts_and_faces(mesh.domain)
        # Apply 'move'/'rotate' options to verts
        if self.rotation is not None:
            verts = self.rotate_verts(verts)
        if self.move is not None:
            verts += self.move
        # Add vertices to surface mesh
        new_vertex_index = mesh.create_vertices(verts)
        # Apply 'side' option to faces
        if self.side == 1:
            if not isinstance(faces, np.ndarray):
                faces = [face[::-1] for face in faces]
            else:
                faces = np.flip(faces, axis = 1)
        # Add faces to surface mesh
        if not isinstance(faces, np.ndarray):
            for _, face in enumerate(faces):
                for j, _ in enumerate(face):
                    face[j] += new_vertex_index
        else:
            faces += new_vertex_index
        new_face_index = mesh.create_faces(faces)
        # Set 'Region' property of new faces
        mesh.faces["Region_"][new_face_index:] = mesh.regions.count - 1

    @abstractmethod
    def generate_verts_and_faces(self, domain: SimulationCell):
        ...

    @abstractmethod
    def is_inside_implementation(self, coords):
        ...

    def is_inside(self, coords):
        if self.rotation is not None:
            coords = self.inv_rotate_verts(coords)
        if self.move is not None:
            coords = coords - self.move
        inside = self.is_inside_implementation(coords)
        if self.side == 1:
            return np.invert(inside)
        return inside

class Block(RegionStyle):
    def __init__(self, tokens: List[str]):
        super().__init__()
        self.xhi, self.yhi, self.zhi, self.xlo, self.ylo, self.zlo = None, None, None, None, None, None
        tokens = self.parse_general_parameters(tokens)
        tokens = self.parse_style_parameters(tokens)
        tokens = self.parse_keyword_parameters(tokens)

    def parse_style_parameters(self, tokens: List[str]):
        '''
        block style args = xlo xhi ylo yhi zlo zhi
        '''
        try:
            self.xlo, self.xhi, self.ylo, self.yhi, self.zlo, self.zhi  = map(lambda v: float(v) if v not in ('INF', 'EDGE') else v, tokens[:6])
        except Exception as exc:
            raise SyntaxError("Expected style args <xlo>, <xhi>, <ylo>, <yhi>, <zlo>, <zhi> after style 'block'.\nNote that LAMMPS equal-style variables <v_> are not supported.") from exc
        return tokens[6:]

    def generate_verts_and_faces(self, domain: SimulationCell):
        cell_diameter_vec = np.sum(domain[:,0:3], axis = 0)
        cell_diameter = np.linalg.norm(cell_diameter_vec)
        cell_center = domain[0:3,3] + cell_diameter_vec/2

        if self.xhi == "EDGE":
            self.xhi = domain[0,3] + np.linalg.norm(domain[0,:3])
        if self.yhi == "EDGE":
            self.yhi = domain[1,3] + np.linalg.norm(domain[1,:3])
        if self.zhi == "EDGE":
            self.zhi = domain[2,3] + np.linalg.norm(domain[2,:3])
        if self.xlo == "EDGE":
            self.xlo = domain[0,3]
        if self.ylo == "EDGE":
            self.ylo = domain[1,3]
        if self.zlo == "EDGE":
            self.zlo = domain[2,3]

        if self.xhi == "INF":
            self.xhi = cell_center[0] + cell_diameter/2
        if self.yhi == "INF":
            self.yhi = cell_center[1] + cell_diameter/2
        if self.zhi == "INF":
            self.zhi = cell_center[2] + cell_diameter/2
        if self.xlo == "INF":
            self.xlo = cell_center[0] - cell_diameter/2
        if self.ylo == "INF":
            self.ylo = cell_center[1] - cell_diameter/2
        if self.zlo == "INF":
            self.zlo = cell_center[2] - cell_diameter/2

        assert self.xhi > self.xlo, "xhi must be > xlo"
        assert self.yhi > self.ylo, "yhi must be > ylo"
        assert self.zhi > self.zlo, "zhi must be > zlo"

        verts = np.array([[0.,0.,0.], [0.,0.,1.], [0.,1.,1.], [0.,1.,0.], [1.,0.,0.], [1.,0.,1.], [1.,1.,1.], [1.,1.,0.]])
        verts[:,0] *= (self.xhi-self.xlo)
        verts[:,0] += self.xlo
        verts[:,1] *= (self.yhi-self.ylo)
        verts[:,1] += self.ylo
        verts[:,2] *= (self.zhi-self.zlo)
        verts[:,2] += self.zlo

        faces = []
        if 1 not in self.open:
            faces.append([0,1,2,3])
        if 2 not in self.open:
            faces.append([7,6,5,4])
        if 3 not in self.open:
            faces.append([0,4,5,1])
        if 4 not in self.open:
            faces.append([3,2,6,7])
        if 5 not in self.open:
            faces.append([0,3,7,4])
        if 6 not in self.open:
            faces.append([1,5,6,2])

        return verts, np.array(faces)

    def is_inside_implementation(self, coords):
        hi = np.all(np.less_equal(coords, [self.xhi, self.yhi, self.zhi]), axis = 1)
        lo = np.all(np.greater_equal(coords, [self.xlo, self.ylo, self.zlo]), axis = 1)
        return np.logical_and(hi,lo)

class Prism(RegionStyle):
    def __init__(self, tokens: List[str]):
        super().__init__()
        self.xlo, self.xhi, self.ylo, self.yhi, self.zlo, self.zhi, self.xy, self.xz, self.yz = None, None, None, None, None, None, None, None, None
        tokens = self.parse_general_parameters(tokens)
        tokens = self.parse_style_parameters(tokens)
        tokens = self.parse_keyword_parameters(tokens)

    def parse_style_parameters(self, tokens: List[str]):
        '''
        prism args = xlo xhi ylo yhi zlo zhi xy xz yz
        xlo,xhi,ylo,yhi,zlo,zhi = bounds of untilted prism (distance units)
        xy = distance to tilt y in x direction (distance units)
        xz = distance to tilt z in x direction (distance units)
        yz = distance to tilt z in y direction (distance units)

        For a prism region, a non-zero tilt factor in any pair of dimensions cannot
        be used if both the lo/hi values in either of those dimensions are INF. E.g. if the xy tilt is non-zero,
        then xlo and xhi cannot both be INF, nor can ylo and yhi.
        '''
        try:
            self.xlo, self.xhi, self.ylo, self.yhi, self.zlo, self.zhi, self.xy, self.xz, self.yz = map(lambda v: float(v) if v not in ('INF', 'EDGE') else v, tokens[:9])
        except Exception as exc:
            raise Exception("Expected style args <xlo>, <xhi>, <ylo>, <yhi>, <zlo>, <zhi>, <xy>, <xz>, <yz> after style 'prism'.\nNote that LAMMPS equal-style variables <v_> are not supported.") from exc

        prism_tilt_error = "A non-zero tilt factor in any pair of dimensions cannot be used if both the lo/hi values in either of those dimensions are INF"
        if self.xy != 0. and ((self.xhi == "INF" and self.xlo == "INF") or (self.yhi == "INF" and self.ylo == "INF")):
            raise Exception(prism_tilt_error)
        if self.xz != 0. and ((self.xhi == "INF" and self.xlo == "INF") or (self.zhi == "INF" and self.zlo == "INF")):
            raise Exception(prism_tilt_error)
        if self.yz != 0. and ((self.zhi == "INF" and self.zlo == "INF") or (self.yhi == "INF" and self.ylo == "INF")):
            raise Exception(prism_tilt_error)
        return tokens[9:]

    def generate_verts_and_faces(self, domain: SimulationCell):
        cell_diameter_vec = np.sum(domain[:,0:3], axis = 0)
        cell_diameter = np.linalg.norm(cell_diameter_vec)
        cell_center = domain[0:3,3] + np.sum(domain[:,0:3], axis = 0)/2

        if self.xhi == "EDGE":
            self.xhi = domain[0,3]+np.linalg.norm(domain[0,:3])
        if self.yhi == "EDGE":
            self.yhi = domain[1,3]+np.linalg.norm(domain[1,:3])
        if self.zhi == "EDGE":
            self.zhi = domain[2,3]+np.linalg.norm(domain[2,:3])
        if self.xlo == "EDGE":
            self.xlo = domain[0,3]
        if self.ylo == "EDGE":
            self.ylo = domain[1,3]
        if self.zlo == "EDGE":
            self.zlo =  domain[2,3]

        if self.xhi == "INF":
            self.xhi = cell_center[0] + cell_diameter
        if self.yhi == "INF":
            self.yhi = cell_center[1] + cell_diameter
        if self.zhi == "INF":
            self.zhi = cell_center[2] + cell_diameter
        if self.xlo == "INF":
            self.xlo = cell_center[0] - cell_diameter
        if self.ylo == "INF":
            self.ylo = cell_center[1] - cell_diameter
        if self.zlo == "INF":
            self.zlo = cell_center[2] - cell_diameter

        assert self.xhi > self.xlo, "xhi must be > xlo"
        assert self.yhi > self.ylo, "yhi must be > ylo"
        assert self.zhi > self.zlo, "zhi must be > zlo"

        verts = np.array([[0.,0.,0.], [0.,0.,1.], [0.,1.,1.], [0.,1.,0.], [1.,0.,0.], [1.,0.,1.],  [1.,1.,1.], [1.,1.,0.]])
        verts[:,0] *= (self.xhi-self.xlo)
        verts[:,1] *= (self.yhi-self.ylo)
        verts[:,2] *= (self.zhi-self.zlo)

        M = np.zeros((3,3))
        M[0][0] = 1.
        M[1][1] = 1.
        M[2][2] = 1.
        M[2][1] = self.yz/(self.zhi-self.zlo)
        M[1][0] = self.xy/(self.yhi-self.ylo)
        M[2][0] = self.xz/(self.zhi-self.zlo)
        verts = (M.T @ verts.T).T
        verts[:,0] += self.xlo
        verts[:,1] += self.ylo
        verts[:,2] += self.zlo

        faces = []
        if 1 not in self.open:
            faces.append([0,1,2,3])
        if 2 not in self.open:
            faces.append([7,6,5,4])
        if 3 not in self.open:
            faces.append([0,4,5,1])
        if 4 not in self.open:
            faces.append([3,2,6,7])
        if 5 not in self.open:
            faces.append([0,3,7,4])
        if 6 not in self.open:
            faces.append([1,5,6,2])

        return verts, np.array(faces)

    def is_inside_implementation(self, coords):
        h = np.array(((self.xhi - self.xlo, self.xy, self.xz),
                    (0., self.yhi - self.ylo, self.yz),
                    (0., 0., self.zhi - self.zlo)
                    ))
        x = np.linalg.inv(h)*(coords - [self.xlo, self.ylo, self.zlo])
        return np.all(np.logical_and(x >=0.0, x <= 1.), axis = 1)

class Cylinder(RegionStyle):
    def __init__(self, tokens: List[str]):
        super().__init__()
        self.dim, self.c1, self.c2, self.radius, self.lo, self.hi = None, None, None, None, None, None
        tokens = self.parse_general_parameters(tokens)
        tokens = self.parse_style_parameters(tokens)
        tokens = self.parse_keyword_parameters(tokens)
        self.resolution = 50

    def parse_style_parameters(self, tokens: List[str]):
        '''
        cylinder args = dim c1 c2 radius lo hi
        dim = x or y or z = axis of cylinder
        c1,c2 = coords of cylinder axis in other 2 dimensions (distance units)
        radius = cylinder radius (distance units)
        c1,c2, and radius can be a variable (see below)
        lo,hi = bounds of cylinder in dim (distance units)
        '''
        self.dim = tokens[0]
        if self.dim not in ["x","y","z"]:
            raise ValueError("Expected dimension {x|y|z}}")
        try:
            self.c1, self.c2, self.radius, self.lo, self.hi = list(map(lambda v: float(v) if v not in ('INF', 'EDGE') else v, tokens[1:6]))
        except Exception as exc:
            raise SyntaxError("Expected style args <dim>, <c1>, <c2>, <radius>, <lo>, <hi> after style 'cylinder'.\nNote that LAMMPS equal-style variables v_ are not supported.") from exc
        return tokens[6:]

    @staticmethod
    def points_on_circle(r, n, dim):
        v = np.zeros(shape=(n, 3))
        i,j = [[1,2], [0,2], [0,1]][["x", "y", "z"].index(dim)]
        v[:,i] = np.cos(np.linspace(0, 2*np.pi, n, endpoint=False)) * r
        v[:,j] = np.sin(np.linspace(0, 2*np.pi, n, endpoint=False)) * r
        return v

    def generate_verts_and_faces(self, domain: SimulationCell):
        cell_diameter_vec = np.sum(domain[:,0:3], axis = 0)
        cell_diameter = np.linalg.norm(cell_diameter_vec)
        cell_center = domain[0:3,3] + np.sum(domain[:,0:3], axis = 0)/2

        i = ["x", "y", "z"].index(self.dim)
        if self.lo == "EDGE" :
            self.lo = domain[i,3]
        if self.hi == "EDGE":
            self.hi = domain[i,3]+np.linalg.norm(domain[i,:3])
        if self.lo == "INF":
            self.lo = cell_center[i] - cell_diameter/2
        if self.hi == "INF":
            self.hi = cell_center[i] + cell_diameter/2
        assert self.hi > self.lo, "hi must be > lo"

        bottom_plane = self.points_on_circle(self.radius, self.resolution, self.dim)
        top_plane = self.points_on_circle(self.radius, self.resolution, self.dim)

        bottom_plane[:,i] += self.lo
        top_plane[:,i] += self.hi

        verts = np.concatenate((top_plane, bottom_plane))
        verts = np.concatenate((top_plane, bottom_plane))
        if self.dim == "x":
            verts[:,1] += self.c1
            verts[:,2] += self.c2
        if self.dim == "y":
            verts[:,0] += self.c1
            verts[:,2] += self.c2
        if self.dim == "z":
            verts[:,0] += self.c1
            verts[:,1] += self.c2

        b_face = np.arange(bottom_plane.shape[0])
        t_face = np.arange(bottom_plane.shape[0]+top_plane.shape[0]-1, bottom_plane.shape[0]-1, -1)
        faces = []
        if not 1 in self.open:
            faces.append(b_face)
        if not 2 in self.open:
            faces.append(t_face)
        if not 3 in self.open:
            for i in range(self.resolution):
                wall_fragment = [(i+1)%self.resolution, i, t_face[::-1][i], t_face[::-1][(i+1)%self.resolution]]
                faces.append(wall_fragment)
        return verts, faces

    def is_inside_implementation(self, coords):
        if self.dim == "x":
            ax = np.logical_and(coords[:,0] >= self.lo, coords[:,0] <= self.hi)
            c = np.less_equal(np.linalg.norm(coords[:,1:] - [self.c1, self.c2], axis = 1), self.radius)
        if self.dim == "y":
            ax = np.logical_and(coords[:,1] >= self.lo, coords[:,1] <= self.hi)
            c = np.less_equal(np.linalg.norm(coords[:,0::2] - [self.c1, self.c2], axis = 1), self.radius)
        if self.dim == "z":
            ax = np.logical_and(coords[:,2] >= self.lo, coords[:,2] <= self.hi)
            c = np.less_equal(np.linalg.norm(coords[:,:2] - [self.c1, self.c2], axis = 1), self.radius)
        return np.logical_and(ax, c)

class Cone(RegionStyle):
    def __init__(self, tokens: List[str]):
        super().__init__()
        self.dim, self.c1, self.c2, self.radlo, self.radhi, self.lo, self.hi = None, None, None, None, None, None, None
        tokens = self.parse_general_parameters(tokens)
        tokens = self.parse_style_parameters(tokens)
        tokens = self.parse_keyword_parameters(tokens)
        self.resolution = 50

    def parse_style_parameters(self, tokens: List[str]):
        '''
        cone style args = dim c1 c2 radlo radhi lo hi
        dim = x or y or z = axis of cone
        c1,c2 = coords of cone axis in other 2 dimensions (distance units)
        radlo,radhi = cone radii at lo and hi end (distance units)
        lo,hi = bounds of cone in dim (distance units)
        '''
        self.dim = tokens[0]
        if self.dim not in ["x","y","z"]:
            raise ValueError("Cone style parameter 'dimension' must be x, y, or z")
        try:
            self.c1, self.c2, self.radlo, self.radhi, self.lo, self.hi = list(map(lambda v: float(v) if v not in ('INF', 'EDGE') else v, tokens[1:7]))
        except Exception as exc:
            raise Exception(f"Expected style args <dim>, <c1>, <c2>, <radius1>, <radius2>, <lo>, <hi> after style 'cone'.\nNote that LAMMPS equal-style variables v_ are not supported.") from exc
        return tokens[7:]

    def generate_verts_and_faces(self, domain: SimulationCell):
        cell_diameter_vec = np.sum(domain[:,0:3], axis = 0)
        cell_diameter = np.linalg.norm(cell_diameter_vec)
        cell_center = domain[0:3,3] + np.sum(domain[:,0:3], axis = 0)/2
        i = ["x", "y", "z"].index(self.dim)
        if self.lo == "EDGE" :
            self.lo = domain[i,3]
        if self.hi == "EDGE":
            self.hi = domain[i,3]+np.linalg.norm(domain[i,:3])
        if self.lo == "INF":
            self.lo = cell_center[i] - cell_diameter/2
        if self.hi == "INF":
            self.hi = cell_center[i] + cell_diameter/2
        assert self.hi > self.lo, "hi must be > lo"

        bottom_plane = Cylinder.points_on_circle(self.radlo, self.resolution, self.dim)
        top_plane = Cylinder.points_on_circle(self.radhi, self.resolution, self.dim)

        i = ["x", "y", "z"].index(self.dim)
        bottom_plane[:,i] += self.lo
        top_plane[:,i] += self.hi

        verts = np.concatenate((bottom_plane, top_plane))
        if self.dim == "x":
            verts[:,1] += self.c1
            verts[:,2] += self.c2
        if self.dim == "y":
            verts[:,0] += self.c1
            verts[:,2] += self.c2
        if self.dim == "z":
            verts[:,0] += self.c1
            verts[:,1] += self.c2

        b_face = np.arange(bottom_plane.shape[0])
        t_face = np.arange(bottom_plane.shape[0]+top_plane.shape[0]-1, bottom_plane.shape[0]-1, -1)
        faces = []
        if not 1 in self.open:
            faces.append(b_face)
        if not 2 in self.open:
            faces.append(t_face)
        if not 3 in self.open:
            for i in range(self.resolution):
                wall_fragment = [(i+1)%self.resolution, i, t_face[::-1][i], t_face[::-1][(i+1)%self.resolution]]
                faces.append(wall_fragment)
        return verts, faces

    def is_inside_implementation(self, coords):
        if self.dim == "x":
            ax = np.logical_and(coords[:,0] >= self.lo, coords[:,0] <= self.hi)
            r = self.radlo + (coords[:,0]-self.lo)*(self.radhi - self.radlo)/(self.hi - self.lo)
            c = np.less_equal(np.linalg.norm(coords[:,1:] - [self.c1, self.c2], axis = 1), r)
        if self.dim == "y":
            ax = np.logical_and(coords[:,1] >= self.lo, coords[:,1] <= self.hi)
            r = self.radlo + (coords[:,1]-self.lo)*(self.radhi - self.radlo)/(self.hi - self.lo)
            c = np.less_equal(np.linalg.norm(coords[:,0::2] - [self.c1, self.c2], axis = 1), r)
        if self.dim == "z":
            ax = np.logical_and(coords[:,2] >= self.lo, coords[:,2] <= self.hi)
            r = self.radlo + (coords[:,2]-self.lo)*(self.radhi - self.radlo)/(self.hi - self.lo)
            c = np.less_equal(np.linalg.norm(coords[:,:2] - [self.c1, self.c2], axis = 1), r)
        return  np.logical_and(ax, c)

class Sphere(RegionStyle):
    def __init__(self, tokens: List[str]):
        super().__init__()
        self.x, self.y, self.z, self.radius = None, None, None, None
        tokens = self.parse_general_parameters(tokens)
        tokens = self.parse_style_parameters(tokens)
        tokens = self.parse_keyword_parameters(tokens)
        self.resolution = 4

    def parse_style_parameters(self, tokens: List[str]):
        '''
        sphere style args = x y z radius
        x,y,z = center of sphere (distance units)
        radius = radius of sphere (distance units)
        '''
        try:
            self.x, self.y, self.z, self.radius  = list(map(float, tokens[:4]))
        except:
            raise Exception(f"Expected <x>, <y>, <z>, <radius> after style 'sphere'.\nNote that LAMMPS equal-style variables v_ are not supported.")
        return tokens[4:]

    @staticmethod
    def createIcoSphere(resolution: int):

        X = 0.525731112119133606
        Z = 0.850650808352039932
        N = 0.0

        vertices = [
            (-X,N,Z), (X,N,Z), (-X,N,-Z), (X,N,-Z),
            (N,Z,X), (N,Z,-X), (N,-Z,X), (N,-Z,-X),
            (Z,X,N), (-Z,X, N), (Z,-X,N), (-Z,-X, N)
        ]

        triangles = [(1, 4, 0), (4, 9, 0), (4, 5, 9), (8, 5, 4), (1, 8, 4),
                    (1, 10, 8), (10, 3, 8), (8, 3, 5), (3, 2, 5), (3, 7, 2),
                    (3, 10, 7), (10, 6, 7), (6, 11, 7), (6, 0, 11), (6, 1, 0),
                    (10, 1, 6), (11, 0, 9), (2, 11, 9), (5, 2, 9), (11, 2, 7)]

        new_vertex_cache = {}
        for i in range(resolution):
            new_faces = []
            for face in triangles:
                mid = [0,0,0]
                for edge in range(3):
                    first = face[edge]
                    second = face[(edge+1)%3]
                    #Has edge been cut before?
                    key =  (first, second) if first < second else (second, first)
                    index = new_vertex_cache.get(key)
                    if not index:
                        edge0 = vertices[first]
                        edge1 = vertices[second]
                        point = np.asarray(edge1) + np.asarray(edge0)
                        point /= np.linalg.norm(point)
                        vertices.append(point)
                        new_vertex_cache[key] = len(vertices) - 1
                    mid[edge] = new_vertex_cache[key]
                new_faces.append((face[0], mid[0], mid[2]))
                new_faces.append((face[1], mid[1], mid[0]))
                new_faces.append((face[2], mid[2], mid[1]))
                new_faces.append((mid[0], mid[1], mid[2]))
            triangles = new_faces

        faces = np.reshape(triangles, (-1, 3))
        verts = np.reshape(vertices, (-1, 3))
        return verts, faces

    def generate_verts_and_faces(self, domain: SimulationCell):
        verts, faces = self.createIcoSphere(self.resolution)
        verts *= self.radius
        verts += [self.x,self.y,self.z]
        return verts, faces

    def is_inside_implementation(self, coords):
        return np.less_equal(np.linalg.norm(coords - [self.x, self.y, self.z], axis = 1), self.radius)

class Ellipsoid(RegionStyle):
    def __init__(self, tokens: List[str]):
        super().__init__()
        self.x, self.y, self.z, self.a, self.b, self.c = None, None, None, None, None, None
        tokens = self.parse_general_parameters(tokens)
        tokens = self.parse_style_parameters(tokens)
        tokens = self.parse_keyword_parameters(tokens)
        self.resolution = 5

    def parse_style_parameters(self, tokens: List[str]):
        '''
        ellipsoid style args = x y z a b c
        x,y,z = center of ellipsoid (distance units)
        a,b,c = half the length of the principal axes of the ellipsoid (distance units)
        x,y,z,a,b and c can be a variable (see below)
        '''
        try:
            self.x, self.y, self.z, self.a, self.b, self.c = list(map(float, tokens[:6]))
        except:
            raise Exception(f"Expected <x>, <y>, <z>, <a>, <b>, <c> after style 'ellipsoid'.\nNote that LAMMPS equal-style variables v_ are not supported.")
        return tokens[6:]

    def generate_verts_and_faces(self, domain: SimulationCell):
        verts, faces = Sphere.createIcoSphere(self.resolution)
        M = np.array(((self.a, 0, 0),(0, self.b, 0),(0, 0, self.c)))
        verts = (verts @ M) + (self.x, self.y, self.z)
        return verts, faces

    def is_inside_implementation(self, coords):
        rc = np.sqrt(self.a * self.a * self.b * self.b * self.c * self.c)
        return np.less_equal(np.linalg.norm((coords - [self.x, self.y, self.z])*[self.b*self.c,self.a*self.c,self.a*self.b], axis = 1), rc)

class Plane(RegionStyle):
    def __init__(self, tokens: List[str]):
        super().__init__()
        self.px, self.py, self.pz, self.nx, self.ny, self.nz  = None, None, None, None, None, None
        tokens = self.parse_general_parameters(tokens)
        tokens = self.parse_style_parameters(tokens)
        tokens = self.parse_keyword_parameters(tokens)
        if self.move is not None:
            raise SyntaxError("Illegal keyword 'move' after region style 'plane'.")
        if self.rotation is not None:
            raise SyntaxError("Illegal keyword 'rotate' after region style 'plane'.")

    def parse_style_parameters(self, tokens: List[str]):
        """
        plane style args = px py pz nx ny nz
        px,py,pz = point on the plane (distance units)
        nx,ny,nz = direction normal to plane (distance units)
        """
        try:
            self.px, self.py, self.pz, self.nx, self.ny, self.nz = list(map(float, tokens[:6]))
        except:
            raise Exception(f"Expected <px>, <py>, <pz>, <nx>, <ny>, <nz> after style 'plane'.\nNote that LAMMPS equal-style variables v_ are not supported.")
        return tokens[6:]

    # intersection function
    @staticmethod
    def intersect_line_plane(p0: npt.ArrayLike, p1: npt.ArrayLike, p: npt.ArrayLike, n: npt.ArrayLike, epsilon=1e-8):
        """
        p0 and p1 are points on the line
        p is a point on the plane (plane coordinate).
        n is the normal vector of the plane
        Return a Vector or None (in case of no intersection)).
        """
        n /= np.linalg.norm(n)
        eij = np.array(p1) - np.array(p0)
        p = np.array(p)
        dot = np.dot(n, eij)
        if abs(dot) > epsilon:
            # The factor of the point between p0 -> p1 (0 - 1)
            # if 'fac' is between (0 - 1) the point intersects with the segment.
            fac = (- np.dot(n, (p0-p))) / dot
            if (fac >= 0.0) and (fac <= 1.0):
                return list(p0 + eij*fac)
        # The segment is parallel to plane.
        return None

    @staticmethod
    def angle_between(v1: npt.ArrayLike, v2: npt.ArrayLike, n: npt.ArrayLike):
        """
        Signed angle between v1 and v2, both lying on the plane with plane normal n
        """
        x = np.cross(v1, v2)
        c = np.sign(np.dot(x,n)) * np.linalg.norm(x)
        return np.arctan2(c,np.dot(v1,v2))

    def generate_verts_and_faces(self, domain: SimulationCell):
        xlo = domain[0,3]
        xhi = xlo + domain[0,0]
        ylo = domain[1,3]
        yhi = ylo + domain[1,1]
        zlo = domain[2,3]
        zhi = zlo + domain[2,2]

        cell_coords = np.array([[xhi, yhi, zhi], [xlo, yhi, zhi], [xhi, yhi, zlo], [xhi, ylo, zhi],\
                                [xlo, ylo, zhi], [xlo, yhi, zlo],[xhi, ylo, zlo], [xlo, ylo, zlo]])
        edges = [[0,1], [1,4], [4, 7], [1,5], [0,2], [2,5], [5,7], [2,6] , [0,3], [3,6], [6,7], [3,4]]

        coords = []
        for edge in edges:
            p = self.intersect_line_plane(cell_coords[edge[0]], cell_coords[edge[1]], (self.px,self.py,self.pz), (self.nx,self.ny,self.nz))
            if (p is not None) and (p not in coords):
                coords.append(p)

        centroid = np.mean(coords, axis = 0)
        angles = np.array([self.angle_between(coords[0]-centroid, coord-centroid, (self.nx,self.ny,self.nz)) for coord in coords])

        if len(coords) > 2:
            return np.array(coords), [np.argsort(angles)]
        else:
            raise RuntimeError("Plane is located outside of simulation cell.")

    def is_inside_implementation(self, coords):
        dot = np.dot(coords - [self.px, self.py, self.pz], [self.nx, self.ny, self.nz])
        return np.greater_equal(dot, 0)

class RenderLAMMPSRegionsModifier(ModifierInterface):
    """
    Base: :py:class:`ovito.pipeline.ModifierInterface`

    This :ref:`Python-based modifier <writing_custom_modifiers>` allows visualizing
    spatial regions with different 3d geometries as defined by the `region <https://docs.lammps.org/region.html>`__ command
    of the `LAMMPS <https://docs.lammps.org/>`__ simulation code.
    See also the corresponding :ref:`user manual page <manual:modifiers.render_lammps_regions>` for more information on this OVITO modifier.

    The :py:class:`~ovito.data.SurfaceMesh` generated by the modifier is accessible
    in the pipeline output :py:class:`~ovito.data.DataCollection` under the unique :py:attr:`~ovito.data.DataObject.identifier` *lammps-regions*::

        data = pipeline.compute()
        mesh = data.surfaces['lammps-regions']

    """

    commands = Code(value='region 1 sphere -20 0 0 15\nregion 2 cylinder z 10 0 5 INF INF', label='', ovito_group='Region commands')
    """
    The `LAMMPS region command(s) <https://docs.lammps.org/region.html>`__ to be interpreted by the modifier.
    If multiple commands are specified, they must be separated by new-line characters.

    :Default: ``'region 1 sphere -20 0 0 15'``
    """

    color_by_type = Bool(default_value=False, label='Color by region', ovito_group='Options')
    """
    Controls whether the modifier should color the parts of the surface mesh by region.
    The effect of this option is the same as if you would additionally insert the :py:class:`~ovito.modifiers.ColorByTypeModifier` into the pipeline.

    If turned off, the mesh will have a uniform color by default and you can make use of the color mapping options of the :py:attr:`vis` element
    to configure a different coloring scheme if necessary.

    If turned on, the modifier will generate the :ref:`mesh region property <surface-mesh-region-properties-list>` ``Color``,
    and the mesh will be rendered using the values stored in this property. The modifier assigns a prescribed set of colors to the regions
    in a cyclic manner.
    Controls whether the modifier should color the parts of the surface mesh by region.
    The effect of this option is the same as if you would additionally insert the :py:class:`~ovito.modifiers.ColorByTypeModifier` into the pipeline.

    If turned off, the mesh will have a uniform color by default and you can make use of the color mapping options of the :py:attr:`vis` element
    to configure a different coloring scheme if necessary.

    If turned on, the modifier will generate the :ref:`mesh region property <surface-mesh-region-properties-list>` ``Color``,
    and the mesh will be rendered using the values stored in this property. The modifier assigns a prescribed set of colors to the regions
    in a cyclic manner.

    :Default: ``False``

    .. versionadded:: 3.10.0
    """

    select_atoms = Bool(default_value=False, label='Select interior particles', ovito_group='Options')
    """
    Turns on the selection of all particles that are located within any of the generated regions.

    :Default: ``False``
    """

    vis = OvitoObject(SurfaceMeshVis, title='LAMMPS Regions', smooth_shading=False, clip_at_domain_boundaries=True)
    """
    The :py:class:`~ovito.vis.SurfaceMeshVis` element that is responsible for rendering the surface mesh generated by the modifier.
    It controls the visual appearance of the region surfaces.
    """

    def modify(self, data: DataCollection, **kwargs):
        if data.cell is None:
            raise RuntimeError('This modifier requires a simulation cell to be defined.')

        # Create SurfaceMesh object or use existing one from another modifier in the pipeline.
        if 'lammps-regions' not in data.surfaces:

            # Disable PBCs for the surface mesh domain.
            # Note: For the time being, we have to clone the entire DataCollection, because the Python API has no way of copying just the SimulationCell object.
            copy = data.clone()
            copy.cell_.pbc = (False, False, False)
            copy.cell_.vis = None

            surface_mesh = data.surfaces.create('lammps-regions', domain=copy.cell, title='LAMMPS Regions', vis=self.vis)
        else:
            surface_mesh = data.surfaces['lammps-regions']
        if 'Region ID' not in surface_mesh.regions:
            surface_mesh.regions_.create_property('Region ID', components=1, dtype=int)
        if "Region" not in surface_mesh.faces:
            surface_mesh.faces_.create_property('Region')

        if self.select_atoms:
            sel = data.particles_.create_property('Selection', data=0)

        for command in self.cleanup_input(self.commands):
            region = self.parse_region_command(str(command))
            region.generate_mesh(surface_mesh)
            if self.select_atoms:
                np.logical_or(sel, region.is_inside(data.particles.positions), out=sel[...])
        surface_mesh.connect_opposite_halfedges()

        # Color mesh faces by region type
        if self.color_by_type:
            data.apply(ColorByTypeModifier(property = 'Region ID', operate_on = 'regions:lammps-regions/regions'))

    def parse_region_command(self, command: str):
        tokens = command.split()
        supported_region_styles = {
            'block': Block,
            'cylinder': Cylinder,
            'sphere': Sphere,
            'cone': Cone,
            'ellipsoid': Ellipsoid,
            'plane': Plane,
            'prism': Prism,
        }
        style_name = tokens[2]
        if not style_name in supported_region_styles:
            raise SyntaxError(f"Unknown or unsupported LAMMPS region style: {style_name}.\nSupported styles are: {list(supported_region_styles.keys())}.")
        return supported_region_styles[style_name](tokens)

    def cleanup_input(self, code: str):
        code = code.strip()
        lines = []
        input = iter(code.split("\n"))
        line = next(input).strip()
        while line:
            l = line
            # If the last printable character on the line is a “&” character, the command is assumed to continue on the next line
            while line.endswith('&'):
                line=next(input).strip()
                l+=line
            # All characters from the first “#” character onward are treated as comment and discarded
            if not l.startswith('#'):
                lines.append(l.replace("&", '').split("#")[0])
            try:
                line = next(input).strip()
            except:
                return lines
        return lines