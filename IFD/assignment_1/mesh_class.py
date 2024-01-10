import numpy as np
import copy

class mesh:
    nPoints = 0
    nFaces = 0
    nCells = 0
    points_arr = []
    faces_arr = []
    cells_arr = []
    nBoundaryPatches = 0
    boundary_patches_arr = []


    def __init__(self, points_list, faces_list, cells_list, moving_list = [], fixed_list = [], FAB_list = []) -> None:
        self.nPoints = len(points_list)
        self.nFaces = len(faces_list)
        self.nCells = len(cells_list)
        cell_label = 0
        for p in points_list:
            individual_point = point(p[0], p[1], p[2])
            self.points_arr.append(individual_point)
        for f in faces_list:
            inter_face = [self.points_arr[i] for i in f]
            individual_face = face(inter_face)
            self.faces_arr.append(individual_face)
        for c in cells_list:
            inter_cell = [self.faces_arr[i] for i in c]
            individual_cell = cell(inter_cell)
            individual_cell.label = cell_label
            cell_label += 1
            self.cells_arr.append(individual_cell)
        
        self.boundary_patches_arr.append(boundary_patches([self.faces_arr[i] for i in moving_list], [self.faces_arr[i] for i in fixed_list], [self.faces_arr[i] for i in FAB_list]))


    def numOfPoints(self):
        return self.nPoints
    def numOfFaces(self):
        return self.nFaces
    def numOfCells(self):
        return self.nCells
    
    def cellVolumes(self):
        vol_list = []
        for c in self.cells_arr:
            vol_list.append(c.cellVolume())
        return vol_list

    def addBoundaryPatch(self):
        None

    def faceInWhichCell(self, f):
        cell_list = []
        for c in self.cells_arr:
            if f in c.faces:
                cell_list.append(c)
        return cell_list


    def boundary_face_attach_Cell(self):
        cell_list = []
        for patch in self.boundary_patches_arr:
            for f in patch.all_faces:
                cell_list.append(self.faceInWhichCell(f)[0])
        return cell_list
    
    def neighbourCell(self, c):
        neighbour_list = []
        faces = c.faces
        for f in faces:
            face_neighbour = self.faceInWhichCell(f)
            for c_i in face_neighbour:
                if c != c_i:
                    neighbour_list.append(c_i)
        return neighbour_list

    def faceOwnerNeighbour(self, f):
        face_neighbours = self.faceInWhichCell(f)
        if len(face_neighbours) == 1:
            return (face_neighbours[0], None)
        else:
            normal = f.faceAreaVector()
            center_1 = face_neighbours[0].cellCenter()
            center_2 = face_neighbours[1].cellCenter()
            c12 = center_2 - center_1
            if np.dot(normal, c12) > 0:
                return (face_neighbours[0], face_neighbours[1])
            else:
                return (face_neighbours[1], face_neighbours[0])
            


class point:
    x = 0
    y = 0
    z = 0
    arr = np.array([x,y,z])
    def __init__(self, x, y, z) -> None:
        self.x = x
        self.y = y
        self.z = z
        self.arr = np.array([x,y,z])
    def move(self, mv_arr):
        self.arr = self.arr + mv_arr
        self.x = self.arr[0]
        self.y = self.arr[1]
        self.z = self.arr[2]
        return None


class face:
    nVertices = 0
    vertices = []
    def __init__(self, vertices) -> None:
        self.nVertices = len(vertices)
        self.vertices = []
        for v in vertices:
            self.vertices.append(v)

    def faceCenter(self):
        center_arr = np.array([0,0,0])
        for v in self.vertices:
            center_arr = center_arr + v.arr
        center_arr = center_arr / self.nVertices
        return center_arr
    
    def setCenterTo(self, mv_arr):
        for v in self.vertices:
            v.move(mv_arr)
        return None

    def setCenterOrigin(self):
        center_arr = self.faceCenter()
        self.setCenterTo(-center_arr)
        return center_arr

    def faceAreaVector(self):
        area_vec = np.array([0,0,0])
        v0 = self.vertices[0].arr
        center_arr = self.setCenterTo(-v0)
        for i in range(1,self.nVertices-1):
            v1 = self.vertices[i].arr
            v2 = self.vertices[i+1].arr
            area_vec = area_vec + np.cross(v1, v2)

        self.setCenterTo(v0)
        return area_vec/2
    
    def faceCentorid(self):
        v0 = self.vertices[0]
        mass_list = []
        for i in range(1,self.nVertices - 1):
            v1 = self.vertices[i]
            v2 = self.vertices[i+1]
            triangle_face = face([v0, v1, v2])
            area = np.linalg.norm(triangle_face.faceAreaVector())
            mass_center = triangle_face.faceCenter()
            mass_list.append(mass(mass_center[0], mass_center[1], mass_center[2], area))
        centroid = np.array([0,0,0])
        total_mass = 0
        for m in mass_list:
            total_mass += m.mass
            centroid = centroid + m.arr*m.mass
        centroid = centroid/total_mass
        return centroid
        

        

class cell:
    nFaces = 0
    faces = []
    label = -1
    def __init__(self, faces) -> None:
        self.nFaces = len(faces)
        self.faces = []
        for f in faces:
            self.faces.append(f)

    def findVertices(self):
        vertices = []
        for f in self.faces:
            for v in f.vertices:
                if (v not in vertices):
                    vertices.append(v)
        return vertices

    def cellCenter(self):
        center_arr = np.array([0,0,0])
        vertices = self.findVertices()
        for v in vertices:
            center_arr = center_arr + v.arr
        center_arr = center_arr/len(vertices)
        return center_arr



    def setCenterOrigin(self):
        center_arr = np.array([0,0,0])
        vertices = self.findVertices()
        for v in vertices:
            center_arr = center_arr + v.arr
        center_arr = center_arr/len(vertices)
        mv_arr = -center_arr
        for v in vertices:
            v.move(mv_arr)
        return center_arr
        
    def setCenterTo(self, mv_arr):
        vertices = self.findVertices()
        for v in vertices:
            v.move(mv_arr)

    def cellVolume(self):
        vol = 0
        center_arr = self.setCenterOrigin()
        for f in self.faces:
            n = f.nVertices
            v2 = f.vertices[0].arr
            for i in range(1,n-1):
                v0 = f.vertices[i].arr
                v1 = f.vertices[i+1].arr
                vol += np.abs(np.dot(v2, np.cross(v0, v1)))
        vol /= 6
        self.setCenterTo(center_arr)
        return vol
    
    def cellCentroid(self):
        centroid = np.array([0,0,0])
        mass_list = []
        total_mass = 0
        center_arr = self.setCenterOrigin()
        for f in self.faces:
            n = f.nVertices
            v2 = f.vertices[0].arr
            for i in range(1,n-1):
                v0 = f.vertices[i].arr
                v1 = f.vertices[i+1].arr
                vol = np.abs(np.dot(v2, np.cross(v0, v1)))
                if vol > 1e-3:
                #Face f is on the opposite, form a tetrahedral with the origin
                    tetrahedral_center = (v0 + v1 + v2)/4
                    mass_list.append(mass(tetrahedral_center[0], tetrahedral_center[1], tetrahedral_center[2], vol))
        for m in mass_list:
            total_mass += m.mass
            centroid = centroid + m.arr*m.mass
        centroid = centroid/total_mass
        self.setCenterTo(center_arr)
        centroid = centroid + center_arr
        return centroid


class boundary_patches:
    nMovingWall = 0
    nFixedWalls = 0
    nFrontAndBack = 0
    nFaces = 0
    moving_faces = []
    fixed_faces = []
    FAB_faces = []
    all_faces = []
    def __init__(self, moving_f_arr, fixed_f_arr, FAB_f_arr) -> None:
        self.nMovingWall = len(moving_f_arr)
        self.nFixedWalls = len(fixed_f_arr)
        self.nFrontAndBack = len(FAB_f_arr)
        self.nFaces = self.nMovingWall + self.nFixedWalls + self.nFrontAndBack
        self.moving_faces = []
        self.fixed_faces = []
        self.FAB_faces = []
        for f in moving_f_arr:
            self.moving_faces.append(f)
        for f in fixed_f_arr:
            self.fixed_faces.append(f)
        for f in FAB_f_arr:
            self.FAB_faces.append(f)
        self.all_faces.extend(self.moving_faces)
        self.all_faces.extend(self.fixed_faces)
        self.all_faces.extend(self.FAB_faces)
        return None
        
class mass:
    x = 0
    y = 0
    z = 0
    arr = np.array([x,y,z])
    mass = 0
    def __init__(self, x, y, z, m) -> None:
        self.x = x
        self.y = y
        self.z = z
        self.arr = np.array([x,y,z])
        self.mass = m

    def move(self, mv_arr):
        self.arr = self.arr + mv_arr
        self.x = self.arr[0]
        self.y = self.arr[1]
        self.z = self.arr[2]
        return None




    

points_list = [[0,0,0], [2,0,0], [1,0,1], [0,0,1], [0,1,0], [1,1,0], [1,1,1], [0,1,1], [0,2,0], [1,2,0], [1,2,1], [0,2,1]]
faces_list = [[0,1,2,3], [0,3,7,4], [5,4,7,6], [5,6,1], [6,2,1], [6,7,3,2], [5,1,0,4], [8,4,7,11], [8,11,10,9], [9,10,6,5], [10,11,7,6],[5,4,8,9]]
cells_list = [[0,1,2,3,4,5,6], [2,7,8,9,10,11]]
moving_list = [0,5]


M = mesh(points_list, faces_list, cells_list,moving_list)
cell = M.cells_arr[0]

print(cell.cellVolume())

f = M.faces_arr[1]
print(f.faceAreaVector())

print(M.neighbourCell(cell)[0].label)

f = M.faces_arr[2]
(ownerCell, neighbourCell) = M.faceOwnerNeighbour(f)
print(f.faceCentorid())
print(M.cells_arr[0].cellCentroid())


